from enum import Enum
from typing import Callable, List, Dict, Any

import pandas as pd
import pm4py
from pandas import DataFrame
from pm4py.objects.log.obj import EventLog, Trace
from tqdm import tqdm

from special4pm.estimation.metrics import get_singletons, get_doubletons, completeness, coverage, \
    sampling_effort_abundance, sampling_effort_incidence, hill_number_asymptotic, entropy_exp, simpson_diversity

class metric_names(Enum):
    NO_OBSERVATIONS_ABUNDANCE = "abundance_no_observations"
    NO_OBSERVATIONS_INCIDENCE = "incidence_no_observations"
    # TODO finalize


class MetricManager(dict):
    # TODO convert to dataclass
    """
    Manages metrics for abundance and incidence models.
    """

    def __init__(self, d0: bool, d1: bool, d2: bool, c0: bool, c1: bool, l_n: list) -> None:
        # reference sample stats
        super().__init__()
        self.reference_sample_abundance = {}
        self.reference_sample_incidence = {}

        self.incidence_current_total_species_count = 0
        self.abundance_current_total_species_count = 0
        self.incidence_sample_size = 0
        self.abundance_sample_size = 0
        self.current_spatial_aggregation = 0

        self["abundance_no_observations"] = [0]
        self["incidence_no_observations"] = [0]
        self["abundance_sum_species_counts"] = [0]
        self["incidence_sum_species_counts"] = [0]
        self["degree_of_aggregation"] = [0]
        self["abundance_singletons"] = [0]
        self["incidence_singletons"] = [0]
        self["abundance_doubletons"] = [0]
        self["incidence_doubletons"] = [0]

        if d0:
            self["abundance_sample_d0"] = [0]
            self["incidence_sample_d0"] = [0]
            self["abundance_estimate_d0"] = [0]
            self["incidence_estimate_d0"] = [0]

        if d1:
            self["abundance_sample_d1"] = [0]
            self["incidence_sample_d1"] = [0]
            self["abundance_estimate_d1"] = [0]
            self["incidence_estimate_d1"] = [0]

        if d2:
            self["abundance_sample_d2"] = [0]
            self["incidence_sample_d2"] = [0]
            self["abundance_estimate_d2"] = [0]
            self["incidence_estimate_d2"] = [0]

        if c0:
            self["abundance_c0"] = [0]
            self["incidence_c0"] = [0]

        if c1:
            self["abundance_c1"] = [0]
            self["incidence_c1"] = [0]

        for l in l_n:
            self["abundance_l_" + str(l)] = [0]
            self["incidence_l_" + str(l)] = [0]


class SpeciesEstimator:
    """
    A class for the estimation of diversity and completeness profiles of trace-based species definitions
    """

    def __init__(self, d0: bool = True, d1: bool = True, d2: bool = True, c0: bool = True,
                 c1: bool = True,
                 l_n: list = [.9, .95, .99], step_size: int | None = None):
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

        self.metrics = {}
        self.species_retrieval = {}

    def register(self, species_id: str, function: Callable) -> None:
        self.species_retrieval[species_id] = function
        self.metrics[species_id] = MetricManager(self.include_d0, self.include_d1, self.include_d2, self.include_c0,
                                                 self.include_c1, self.l_n)

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
            for species_id in self.species_retrieval.keys():
                for tr in tqdm(data, "Profiling Log for " + species_id):
                    self.add_observation(tr, species_id)
                    # if step size is set, update metrics after <step_size> many traces
                    if self.step_size is None:
                        continue
                    elif self.metrics[species_id].incidence_sample_size % self.step_size == 0:
                        self.update_metrics(species_id)
                self.update_metrics(species_id)

                #self.apply(tr)
            return
        if isinstance(data, Trace):
            pass
        #    #for species_id in self.species_retrieval.keys():
        #    self.add_observation(data, species_id)
        #    # if step size is set, update metrics after <step_size> many traces
        #    if self.step_size is None:
        #        return
        #    elif self.metrics[species_id].incidence_sample_size % self.step_size == 0:
        #        self.update_metrics(species_id)
        else:
            raise RuntimeError('Cannot apply data of type ' + str(type(data)))

    def add_observation(self, observation: Trace, species_id: str) -> None:
        """
        adds a single observation
        :param observation: the trace observation
        """
        # retrieve species from current observation
        species_abundance = self.species_retrieval[species_id](observation)
        species_incidence = set(species_abundance)
        self.metrics[species_id].trace_retrieved_species_abundance = species_abundance
        self.metrics[species_id].trace_retrieved_species_incidence = species_incidence

        # update species abundances/incidences
        for s in species_abundance:
            self.metrics[species_id].reference_sample_abundance[s] = self.metrics[
                                                                         species_id].reference_sample_abundance.get(s,
                                                                                                                    0) + 1

        for s in species_incidence:
            self.metrics[species_id].reference_sample_incidence[s] = self.metrics[
                                                                         species_id].reference_sample_incidence.get(s,
                                                                                                                    0) + 1

        # update current number of observation for each model
        self.metrics[species_id].abundance_sample_size = self.metrics[species_id].abundance_sample_size + len(
            species_abundance)
        self.metrics[species_id].incidence_sample_size = self.metrics[species_id].incidence_sample_size + 1

        # update current sum of all observed species for each model
        self.metrics[species_id].abundance_current_total_species_count = \
            self.metrics[species_id].abundance_current_total_species_count + len(species_abundance)
        self.metrics[species_id].incidence_current_total_species_count = \
            self.metrics[species_id].incidence_current_total_species_count + len(
                species_incidence)

        #update current degree of spatial aggregation
        self.metrics[species_id].current_spatial_aggregation = 1 - (
                self.metrics[species_id].incidence_current_total_species_count / self.metrics[
            species_id].abundance_current_total_species_count)

    def update_metrics(self, species_id: str) -> None:
        """
        updates the diversity and completeness profiles based on the current observations
        """
        # update number of observations so far
        self.metrics[species_id]["abundance_no_observations"].append(self.metrics[species_id].abundance_sample_size)
        self.metrics[species_id]["incidence_no_observations"].append(self.metrics[species_id].incidence_sample_size)

        #update number of species seen so far
        self.metrics[species_id]["abundance_sum_species_counts"].append(
            self.metrics[species_id].abundance_current_total_species_count)
        self.metrics[species_id]["incidence_sum_species_counts"].append(
            self.metrics[species_id].incidence_current_total_species_count)

        #update degree of spatial aggregation
        self.metrics[species_id]["degree_of_aggregation"].append(self.metrics[species_id].current_spatial_aggregation)

        #update singleton and doubleton counts
        self.metrics[species_id]["abundance_singletons"].append(
            get_singletons(self.metrics[species_id].reference_sample_abundance))
        self.metrics[species_id]["incidence_singletons"].append(
            get_singletons(self.metrics[species_id].reference_sample_incidence))

        self.metrics[species_id]["abundance_doubletons"].append(
            get_doubletons(self.metrics[species_id].reference_sample_abundance))
        self.metrics[species_id]["incidence_doubletons"].append(
            get_doubletons(self.metrics[species_id].reference_sample_incidence))

        #update diversity profile
        if self.include_d0:
            self.__update_d0(species_id)
        if self.include_d1:
            self.__update_d1(species_id)
        if self.include_d2:
            self.__update_d2(species_id)

        #update completeness profile
        if self.include_c0:
            self.__update_c0(species_id)
        if self.include_c1:
            self.__update_c1(species_id)

        #update estimated sampling effort for target completeness
        for l in self.l_n:
            self.__update_l(l, species_id)

    def __update_d0(self, species_id: str) -> None:
        """
        updates D0 (=species richness) based on the current observations
        """
        #update sample metrics
        self.metrics[species_id]["abundance_sample_d0"].append(len(self.metrics[species_id].reference_sample_abundance))
        self.metrics[species_id]["incidence_sample_d0"].append(len(self.metrics[species_id].reference_sample_incidence))

        #update estimated metrics
        self.metrics[species_id]["abundance_estimate_d0"].append(
            hill_number_asymptotic(0, self.metrics[species_id].reference_sample_abundance,
                                   self.metrics[species_id].abundance_sample_size))
        self.metrics[species_id]["incidence_estimate_d0"].append(
            hill_number_asymptotic(0, self.metrics[species_id].reference_sample_incidence,
                                   self.metrics[species_id].incidence_sample_size, abundance=False))

    def __update_d1(self, species_id: str) -> None:
        """
        updates D1 (=exponential of Shannon entropy) based on the current observations
        """
        #update sample metrics
        self.metrics[species_id]["abundance_sample_d1"].append(
            entropy_exp(self.metrics[species_id].reference_sample_abundance))
        self.metrics[species_id]["incidence_sample_d1"].append(
            entropy_exp(self.metrics[species_id].reference_sample_incidence))

        #update estimated metrics
        self.metrics[species_id]["abundance_estimate_d1"].append(
            hill_number_asymptotic(1, self.metrics[species_id].reference_sample_abundance,
                                   self.metrics[species_id].abundance_sample_size))
        self.metrics[species_id]["incidence_estimate_d1"].append(
            hill_number_asymptotic(1, self.metrics[species_id].reference_sample_incidence,
                                   self.metrics[species_id].incidence_sample_size, abundance=False))

    def __update_d2(self, species_id: str) -> None:
        """
        updates D2 (=Simpson Diversity Index) based on the current observations
        """
        #update sample metrics
        self.metrics[species_id]["abundance_sample_d2"].append(
            simpson_diversity(self.metrics[species_id].reference_sample_abundance))
        self.metrics[species_id]["incidence_sample_d2"].append(
            simpson_diversity(self.metrics[species_id].reference_sample_incidence))

        #update estimated metrics
        self.metrics[species_id]["abundance_estimate_d2"].append(
            hill_number_asymptotic(2, self.metrics[species_id].reference_sample_abundance,
                                   self.metrics[species_id].abundance_sample_size))
        self.metrics[species_id]["incidence_estimate_d2"].append(
            hill_number_asymptotic(2, self.metrics[species_id].reference_sample_incidence,
                                   self.metrics[species_id].incidence_sample_size, abundance=False))

    def __update_c0(self, species_id: str) -> None:
        """
        updates C0 (=completeness) based on the current observations
        """
        self.metrics[species_id]["abundance_c0"].append(
            completeness(self.metrics[species_id].reference_sample_abundance))
        self.metrics[species_id]["incidence_c0"].append(
            completeness(self.metrics[species_id].reference_sample_incidence))

    def __update_c1(self, species_id: str) -> None:
        """
        updates C1 (=coverage) based on the current observations
        """
        self.metrics[species_id]["abundance_c1"].append(
            coverage(self.metrics[species_id].reference_sample_abundance,
                     self.metrics[species_id].abundance_sample_size))
        self.metrics[species_id]["incidence_c1"].append(
            coverage(self.metrics[species_id].reference_sample_incidence,
                     self.metrics[species_id].incidence_sample_size))

    def __update_l(self, g: float, species_id: str) -> None:
        """
        updates l_g (=expected number additional observations for reaching completeness g) based on the current
        observations
        :param g: desired  completeness
        """
        self.metrics[species_id]["abundance_l_" + str(g)].append(
            sampling_effort_abundance(g, self.metrics[species_id].reference_sample_abundance,
                                      self.metrics[species_id].abundance_sample_size))
        self.metrics[species_id]["incidence_l_" + str(g)].append(
            sampling_effort_incidence(g, self.metrics[species_id].reference_sample_incidence,
                                      self.metrics[species_id].incidence_sample_size))

    def print_metrics(self) -> None:
        """
        prints the Diversity and Completeness Profile of the current observations
        """
        for species_id in self.species_retrieval:
            print("### " + species_id + " ###")
            print("### SAMPLE STATS ###")
            print("Abundance")
            print("%-30s %s" % ("     No Observations:", str(self.metrics[species_id]["abundance_no_observations"])))
            print("%-30s %s" % (
                "     Total Species Count:", str(self.metrics[species_id]["abundance_sum_species_counts"])))
            print("%-30s %s" % ("     Singletons:", str(self.metrics[species_id]["abundance_singletons"])))
            print("%-30s %s" % ("     Doubletons:", str(self.metrics[species_id]["abundance_doubletons"])))
            print("Incidence")
            print("%-30s %s" % ("     No Observations:", str(self.metrics[species_id]["incidence_no_observations"])))
            print("%-30s %s" % (
                "     Total Species Count:", str(self.metrics[species_id]["incidence_sum_species_counts"])))
            print("%-30s %s" % ("     Singletons:", str(self.metrics[species_id]["incidence_singletons"])))
            print("%-30s %s" % ("     Doubletons:", str(self.metrics[species_id]["incidence_doubletons"])))
            print("%-30s %s" % ("Degree of Aggregation:", str(self.metrics[species_id]["degree_of_aggregation"])))
            print()
            print("### DIVERSITY AND COMPLETENESS PROFILE ###")
            print("Abundance")
            if self.include_d0:
                print("%-30s %s" % ("     D0 - sample:", (self.metrics[species_id]["abundance_sample_d0"])))
                print("%-30s %s" % ("     D0 - estimate:", str(self.metrics[species_id]["abundance_estimate_d0"])))
            if self.include_d1:
                print("%-30s %s" % ("     D1 - sample:", str(self.metrics[species_id]["abundance_sample_d1"])))
                print("%-30s %s" % ("     D1 - estimate:", str(self.metrics[species_id]["abundance_estimate_d1"])))
            if self.include_d2:
                print("%-30s %s" % ("     D2 - sample:", str(self.metrics[species_id]["abundance_sample_d2"])))
                print("%-30s %s" % ("     D2 - estimate:", str(self.metrics[species_id]["abundance_estimate_d2"])))
            if self.include_c0:
                print("%-30s %s" % ("     C0:", str(self.metrics[species_id]["abundance_c0"])))
            if self.include_c1:
                print("%-30s %s" % ("     C1:", str(self.metrics[species_id]["abundance_c1"])))
            for l in self.l_n:
                print("%-30s %s" % ("     l_" + str(l) + ":", str(self.metrics[species_id]["abundance_l_" + str(l)])))
            print("Incidence")
            if self.include_d0:
                print("%-30s %s" % ("     D0 - sample:", (self.metrics[species_id]["incidence_sample_d0"])))
                print("%-30s %s" % ("     D0 - estimate:", str(self.metrics[species_id]["incidence_estimate_d0"])))
            if self.include_d1:
                print("%-30s %s" % ("     D1 - sample:", str(self.metrics[species_id]["incidence_sample_d1"])))
                print("%-30s %s" % ("     D1 - estimate:", str(self.metrics[species_id]["incidence_estimate_d1"])))
            if self.include_d2:
                print("%-30s %s" % ("     D2 - sample:", str(self.metrics[species_id]["incidence_sample_d2"])))
                print("%-30s %s" % ("     D2 - estimate:", str(self.metrics[species_id]["incidence_estimate_d2"])))
            if self.include_c0:
                print("%-30s %s" % ("     C0:", str(self.metrics[species_id]["incidence_c0"])))
            if self.include_c1:
                print("%-30s %s" % ("     C1:", str(self.metrics[species_id]["incidence_c1"])))
            for l in self.l_n:
                print("%-30s %s" % ("     l_" + str(l) + ":", str(self.metrics[species_id]["incidence_l_" + str(l)])))
            print()

    def to_dataFrame(self) -> DataFrame:
        """
        returns the diversity and completeness profile of the current observations as a data frame
        :returns: a data frame view of the Diversity and Completeness Profile
        """
        return pd.DataFrame([[i, j, ix, v]
                             for i in self.metrics.keys()
                             for j in self.metrics[i].keys()
                             for ix, v in enumerate(self.metrics[i][j])
                             ], columns=["species", "metric", "observation", "value"]
                            )
