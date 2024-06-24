import pandas as pd
import pm4py
from pandas import DataFrame
from pm4py.objects.log.obj import EventLog, Trace

from special.estimation.metrics import get_singletons, get_doubletons, completeness, coverage, \
    sampling_effort_abundance, sampling_effort_incidence, hill_number_asymptotic, entropy_exp, simpson_diversity


class SpeciesEstimator:
    """
    A class for the estimation of diversity and completeness profiles of trace-based species definitions
    """

    def __init__(self, species_retrieval_function, d0: bool = True, d1: bool =False, d2: bool =False, c0: bool =True, c1: bool =True,
                 l_n: list =[.9, .95, .99], step_size: int|None=None):
        """
        :param species_retrieval_function: a function mapping a trace to a list of corresponding species
        :param d0: flag indicating if D0(=species richness) should be included
        :param d1: flag indicating if D1(=exponential Shannon entropy) should be included
        :param d2: flag indicating if D2(=Simpson diversity index) should be included
        :param c0: flag indicating if C0(=completeness) should be included
        :param c1: flag indicating if C1(=coverage) should be included
        :param l_n: list of desired completeness values for estimation additional sampling effort
        :param step_size: the number of added traces after which the profiles are updated. Use None if
        """
        # TODO add differentiation between abundance and incidence based data
        self.include_abundance = True
        self.include_incidence = True

        self.include_d0 = d0
        self.include_d1 = d1
        self.include_d2 = d2

        self.include_c0 = c0
        self.include_c1 = c1

        self.l_n = l_n

        self.step_size = step_size

        self.species_retrieval = species_retrieval_function

        # reference sample stats
        self.reference_sample_abundance = {}
        self.reference_sample_incidence = {}

        self.incidence_current_total_species_count = 0
        self.abundance_current_total_species_count = 0
        self.incidence_sample_size = 0
        self.abundance_sample_size = 0
        self.current_spatial_aggregation = 0

        self.metrics = {"abundance_no_observations": [], "incidence_no_observations": [],
                        "abundance_sum_species_counts": [], "incidence_sum_species_counts": [],
                        "degree_of_aggregation": [], "abundance_singletons": [], "incidence_singletons": [],
                        "abundance_doubletons": [], "incidence_doubletons": []}

        if self.include_d0:
            self.metrics["abundance_sample_d0"] = []
            self.metrics["incidence_sample_d0"] = []
            self.metrics["abundance_estimate_d0"] = []
            self.metrics["incidence_estimate_d0"] = []

        if self.include_d1:
            self.metrics["abundance_sample_d1"] = []
            self.metrics["incidence_sample_d1"] = []
            self.metrics["abundance_estimate_d1"] = []
            self.metrics["incidence_estimate_d1"] = []

        if self.include_d2:
            self.metrics["abundance_sample_d2"] = []
            self.metrics["incidence_sample_d2"] = []
            self.metrics["abundance_estimate_d2"] = []
            self.metrics["incidence_estimate_d2"] = []

        if self.include_c0:
            self.metrics["abundance_c0"] = []
            self.metrics["incidence_c0"] = []

        if self.include_c1:
            self.metrics["abundance_c1"] = []
            self.metrics["incidence_c1"] = []

        for l in self.l_n:
            self.metrics["abundance_l_" + str(l)] = []
            self.metrics["incidence_l_" + str(l)] = []

    def apply(self, data: pd.DataFrame | EventLog | Trace) -> None:
        """
        add all observations of an event log and update diversity and completeness profiles once afterward.
        If parameter step_size is set to an int, profiles are additionally updated along the way according to
        the step size
        :param data: the event log containing the trace observations
        """
        if isinstance(data, pd.DataFrame):
            return self.apply(pm4py.convert_to_event_log(data))
        if isinstance(data, EventLog):
            for tr in data:
                return self.apply(tr)
        if isinstance(data, Trace):
            self.add_observation(data)
            # if step size is set, update metrics after <step_size> many traces
            if self.step_size is None:
                return
            elif self.incidence_sample_size % self.step_size == 0:
                self.update_metrics()


    def add_observation(self, observation: Trace) -> None:
        """
        adds a single observation
        :param observation: the trace observation
        """
        # retrieve species from current observation
        trace_retrieved_species_abundance = self.species_retrieval(observation)
        trace_retrieved_species_incidence = set(trace_retrieved_species_abundance)

        # update species abundances/incidences
        for s in trace_retrieved_species_abundance:
            self.reference_sample_abundance[s] = self.reference_sample_abundance.get(s, 0) + 1

        for s in trace_retrieved_species_incidence:
            self.reference_sample_incidence[s] = self.reference_sample_incidence.get(s, 0) + 1

        # update current number of observation for each model
        self.abundance_sample_size = self.abundance_sample_size + len(
            trace_retrieved_species_abundance)
        self.incidence_sample_size = self.incidence_sample_size + 1

        # update current sum of all observed species for each model
        self.abundance_current_total_species_count = self.abundance_current_total_species_count + len(
            trace_retrieved_species_abundance)
        self.incidence_current_total_species_count = self.incidence_current_total_species_count + len(
            trace_retrieved_species_incidence)

        #update current degree of spatial aggregation
        self.current_spatial_aggregation = 1 - self.incidence_current_total_species_count / self.abundance_current_total_species_count



    def update_metrics(self) -> None:
        """
        updates the diversity and completeness profiles based on the current observations
        """
        # update number of observations so far
        self.metrics["abundance_no_observations"].append(self.abundance_sample_size)
        self.metrics["incidence_no_observations"].append(self.incidence_sample_size)

        #update number of species seen so far
        self.metrics["abundance_sum_species_counts"].append(self.abundance_current_total_species_count)
        self.metrics["incidence_sum_species_counts"].append(self.incidence_current_total_species_count)

        #update degree of spatial aggregation
        self.metrics["degree_of_aggregation"].append(self.current_spatial_aggregation)

        #update singleton and doubleton counts
        self.metrics["abundance_singletons"].append(get_singletons(self.reference_sample_abundance))
        self.metrics["incidence_singletons"].append(get_singletons(self.reference_sample_incidence))

        self.metrics["abundance_doubletons"].append(get_doubletons(self.reference_sample_abundance))
        self.metrics["incidence_doubletons"].append(get_doubletons(self.reference_sample_incidence))

        #update diversity profile
        if self.include_d0:
            self.__update_d0()
        if self.include_d1:
            self.__update_d1()
        if self.include_d2:
            self.__update_d2()

        #update completeness profile
        if self.include_c0:
            self.__update_c0()
        if self.include_c1:
            self.__update_c1()

        #update estimated sampling effort for target completeness
        for l in self.l_n:
            self.__update_l(l)

    def __update_d0(self) -> None:
        """
        updates D0 (=species richness) based on the current observations
        """
        #update sample metrics
        self.metrics["abundance_sample_d0"].append(len(self.reference_sample_abundance))
        self.metrics["incidence_sample_d0"].append(len(self.reference_sample_incidence))

        #update estimated metrics
        self.metrics["abundance_estimate_d0"].append(
            hill_number_asymptotic(0, self.reference_sample_abundance, self.abundance_sample_size))
        self.metrics["incidence_estimate_d0"].append(
            hill_number_asymptotic(0, self.reference_sample_incidence, self.incidence_sample_size, abundance=False))

    def __update_d1(self) -> None:
        """
        updates D1 (=exponential of Shannon entropy) based on the current observations
        """
        #update sample metrics
        self.metrics["abundance_sample_d1"].append(entropy_exp(self.reference_sample_abundance))
        self.metrics["incidence_sample_d1"].append(entropy_exp(self.reference_sample_incidence))

        #update estimated metrics
        self.metrics["abundance_estimate_d1"].append(
            hill_number_asymptotic(1, self.reference_sample_abundance, self.abundance_sample_size))
        self.metrics["incidence_estimate_d1"].append(
            hill_number_asymptotic(1, self.reference_sample_incidence, self.incidence_sample_size, abundance=False))

    def __update_d2(self) -> None:
        """
        updates D2 (=Simpson Diversity Index) based on the current observations
        """
        #update sample metrics
        self.metrics["abundance_sample_d2"].append(simpson_diversity(self.reference_sample_abundance))
        self.metrics["incidence_sample_d2"].append(simpson_diversity(self.reference_sample_incidence))

        #update estimated metrics
        self.metrics["abundance_estimate_d2"].append(
            hill_number_asymptotic(2, self.reference_sample_abundance, self.abundance_sample_size))
        self.metrics["incidence_estimate_d2"].append(
            hill_number_asymptotic(2, self.reference_sample_incidence, self.incidence_sample_size, abundance=False))

    def __update_c0(self) -> None:
        """
        updates C0 (=completeness) based on the current observations
        """
        self.metrics["abundance_c0"].append(completeness(self.reference_sample_abundance))
        self.metrics["incidence_c0"].append(completeness(self.reference_sample_incidence))

    def __update_c1(self) -> None:
        """
        updates C1 (=coverage) based on the current observations
        """
        self.metrics["abundance_c1"].append(
            coverage(self.reference_sample_abundance, self.abundance_sample_size))
        self.metrics["incidence_c1"].append(
            coverage(self.reference_sample_incidence, self.incidence_sample_size))

    def __update_l(self, g: float) -> None:
        """
        updates l_g (=expected number additional observations for reaching completeness g) based on the current
        observations
        :param g: desired  completeness
        """
        self.metrics["abundance_l_" + str(g)].append(
            sampling_effort_abundance(g, self.reference_sample_abundance, self.abundance_sample_size))
        self.metrics["incidence_l_" + str(g)].append(
            sampling_effort_incidence(g, self.reference_sample_incidence, self.incidence_sample_size))

    def print_metrics(self) -> None:
        """
        prints the Diversity and Completeness Profile of the current observations
        """
        print("### SAMPLE STATS ###")
        print("Abundance")
        print("%-30s %s" % ("     No Observations:", str(self.metrics["abundance_no_observations"])))
        print("%-30s %s" % ("     Total Species Count:", str(self.metrics["abundance_sum_species_counts"])))
        print("%-30s %s" % ("     Singletons:", str(self.metrics["abundance_singletons"])))
        print("%-30s %s" % ("     Doubletons:", str(self.metrics["abundance_doubletons"])))
        print("Incidence")
        print("%-30s %s" % ("     No Observations:", str(self.metrics["incidence_no_observations"])))
        print("%-30s %s" % ("     Total Species Count:", str(self.metrics["incidence_sum_species_counts"])))
        print("%-30s %s" % ("     Singletons:", str(self.metrics["incidence_singletons"])))
        print("%-30s %s" % ("     Doubletons:", str(self.metrics["incidence_doubletons"])))
        print("%-30s %s" % ("Degree of Aggregation:", str(self.metrics["degree_of_aggregation"])))
        print()
        print("### DIVERSITY AND COMPLETENESS PROFILE ###")
        print("Abundance")
        if self.include_d0:
            print("%-30s %s" % ("     D0 - sample:", (self.metrics["abundance_sample_d0"])))
            print("%-30s %s" % ("     D0 - estimate:", str(self.metrics["abundance_estimate_d0"])))
        if self.include_d1:
            print("%-30s %s" % ("     D1 - sample:", str(self.metrics["abundance_sample_d1"])))
            print("%-30s %s" % ("     D1 - estimate:", str(self.metrics["abundance_estimate_d1"])))
        if self.include_d2:
            print("%-30s %s" % ("     D2 - sample:", str(self.metrics["abundance_sample_d2"])))
            print("%-30s %s" % ("     D2 - estimate:", str(self.metrics["abundance_estimate_d2"])))
        if self.include_c0:
            print("%-30s %s" % ("     C0:", str(self.metrics["abundance_c0"])))
        if self.include_c1:
            print("%-30s %s" % ("     C1:", str(self.metrics["abundance_c1"])))
        for l in self.l_n:
            print("%-30s %s" % ("     l_" + str(l) + ":", str(self.metrics["abundance_l_" + str(l)])))
        print("Incidence")
        if self.include_d0:
            print("%-30s %s" % ("     D0 - sample:", (self.metrics["incidence_sample_d0"])))
            print("%-30s %s" % ("     D0 - estimate:", str(self.metrics["incidence_estimate_d0"])))
        if self.include_d1:
            print("%-30s %s" % ("     D1 - sample:", str(self.metrics["incidence_sample_d1"])))
            print("%-30s %s" % ("     D1 - estimate:", str(self.metrics["incidence_estimate_d1"])))
        if self.include_d2:
            print("%-30s %s" % ("     D2 - sample:", str(self.metrics["incidence_sample_d2"])))
            print("%-30s %s" % ("     D2 - estimate:", str(self.metrics["incidence_estimate_d2"])))
        if self.include_c0:
            print("%-30s %s" % ("     C0:", str(self.metrics["incidence_c0"])))
        if self.include_c1:
            print("%-30s %s" % ("     C1:", str(self.metrics["incidence_c1"])))
        for l in self.l_n:
            print("%-30s %s" % ("     l_" + str(l) + ":", str(self.metrics["incidence_l_" + str(l)])))

    def to_dataFrame(self) -> DataFrame:
        """
        returns the diversity and completeness profile of the current observations as a data frame
        :returns: a data frame view of the Diversity and Completeness Profile
        """
        return pd.DataFrame.from_dict(self.metrics)
