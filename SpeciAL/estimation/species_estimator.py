from typing import Callable, List, Dict, Union, Any

import pandas as pd
import pm4py
from pm4py.objects.log.obj import Trace, EventLog

from SpeciAL.estimation.metrics import get_singletons, get_doubletons, hill_number_asymptotic, entropy_exp, \
    simpson_diversity, completeness, coverage, sampling_effort_abundance, sampling_effort_incidence


class MetricManager:
    """
    Manages metrics for abundance and incidence models.
    """
    def __init__(self, metrics: List[str]) -> None:
        self.data: Dict[str, Dict[str, List[Any]]] = {
            metric: {
                "abundance_sample": [],
                "abundance_estimate": [],
                "incidence_sample": [],
                "incidence_estimate": []
            } for metric in metrics
        }

    def update_metric(self, metric: str, data_type: str, value: Any, model_type: str) -> None:
        """
        Update the metric with a new value.

        :param metric: The name of the metric (e.g., 'd0', 'c0').
        :param data_type: The type of data (e.g., 'sample', 'estimate').
        :param value: The value to append to the metric list.
        :param model_type: The model type ('abundance' or 'incidence').
        """
        key = f"{model_type}_{data_type}"
        if metric in self.data and key in self.data[metric]:
            self.data[metric][key].append(value)

    def get_metrics(self) -> Dict[str, Dict[str, List[Any]]]:
        """
        Get all metrics data.

        :return: The dictionary containing all metrics.
        """
        return self.data


class SpeciesEstimator:
    """
    Estimate diversity and completeness profiles of trace-based species definitions.
    """

    def __init__(self, species_retrieval_function: Callable[[Any], List[str]],
                 d0: bool = True, d1: bool = False, d2: bool = False,
                 c0: bool = True, c1: bool = True,
                 l_n: List[float] = [.9, .95, .99], step_size: Union[int, None] = None) -> None:
        """
        Initialize the SpeciesEstimator.

        :param species_retrieval_function: Function mapping a trace to a list of corresponding species.
        :param d0: Include D0 (species richness).
        :param d1: Include D1 (exponential Shannon entropy).
        :param d2: Include D2 (Simpson diversity index).
        :param c0: Include C0 (completeness).
        :param c1: Include C1 (coverage).
        :param l_n: List of desired completeness values for estimation.
        :param step_size: Number of added traces after which the profiles are updated.
        """
        self.step_size = step_size
        self.species_retrieval = species_retrieval_function

        self.metrics_to_include = []
        if d0: self.metrics_to_include.append("d0")
        if d1: self.metrics_to_include.append("d1")
        if d2: self.metrics_to_include.append("d2")
        if c0: self.metrics_to_include.append("c0")
        if c1: self.metrics_to_include.append("c1")

        self.metric_manager = MetricManager(self.metrics_to_include)

        self.reference_sample_abundance: Dict[str, int] = {}
        self.reference_sample_incidence: Dict[str, int] = {}

        self.incidence_current_total_species_count = 0
        self.abundance_current_total_species_count = 0
        self.incidence_sample_size = 0
        self.abundance_sample_size = 0
        self.current_spatial_aggregation = 0

        self.include_d0 = d0
        self.include_d1 = d1
        self.include_d2 = d2
        self.include_c0 = c0

        self.include_c1 = c1
        self.l_n = l_n

    def profile_log(self, log: Union[EventLog, pd.DataFrame]) -> None:
        """
        Add all observations of an event log and update diversity and completeness profiles.

        :param log: The event log containing the trace observations.
        """
        if isinstance(log, pd.DataFrame):
            print("Provided log is in DataFrame format. Converting to EventLog type.")
            log = pm4py.convert_to_event_log(log)
        for tr in log:
            self.add_observation(tr)
            if self.step_size is not None and self.incidence_sample_size % self.step_size == 0:
                self.update_metrics()
        self.update_metrics()

    def add_observation(self, observation: Trace) -> None:
        """
        Add a single observation.

        :param observation: The trace observation.
        """
        trace_retrieved_species_abundance = self.species_retrieval(observation)
        trace_retrieved_species_incidence = set(trace_retrieved_species_abundance)

        for s in trace_retrieved_species_abundance:
            self.reference_sample_abundance[s] = self.reference_sample_abundance.get(s, 0) + 1

        for s in trace_retrieved_species_incidence:
            self.reference_sample_incidence[s] = self.reference_sample_incidence.get(s, 0) + 1

        self.abundance_sample_size += len(trace_retrieved_species_abundance)
        self.incidence_sample_size += 1

        self.abundance_current_total_species_count += len(trace_retrieved_species_abundance)
        self.incidence_current_total_species_count += len(trace_retrieved_species_incidence)

        self.current_spatial_aggregation = 1 - self.incidence_current_total_species_count / self.abundance_current_total_species_count

    def update_metrics(self) -> None:
        """
        Update the diversity and completeness profiles based on the current observations.
        """
        # Update number of observations
        self.metric_manager.update_metric("abundance_no_observations", "sample", self.abundance_sample_size,
                                          "abundance")
        self.metric_manager.update_metric("incidence_no_observations", "sample", self.incidence_sample_size,
                                          "incidence")

        # Update number of species seen so far
        self.metric_manager.update_metric("abundance_sum_species_counts", "sample",
                                          self.abundance_current_total_species_count, "abundance")
        self.metric_manager.update_metric("incidence_sum_species_counts", "sample",
                                          self.incidence_current_total_species_count, "incidence")

        # Update degree of spatial aggregation
        self.metric_manager.update_metric("degree_of_aggregation", "sample", self.current_spatial_aggregation,
                                          "abundance")

        # Update singleton and doubleton counts
        self.metric_manager.update_metric("abundance_singletons", "sample",
                                          get_singletons(self.reference_sample_abundance), "abundance")
        self.metric_manager.update_metric("incidence_singletons", "sample",
                                          get_singletons(self.reference_sample_incidence), "incidence")

        self.metric_manager.update_metric("abundance_doubletons", "sample",
                                          get_doubletons(self.reference_sample_abundance), "abundance")
        self.metric_manager.update_metric("incidence_doubletons", "sample",
                                          get_doubletons(self.reference_sample_incidence), "incidence")

        # Update diversity profiles
        if self.include_d0:
            self.metric_manager.update_metric("d0", "sample", len(self.reference_sample_abundance), "abundance")
            self.metric_manager.update_metric("d0", "sample", len(self.reference_sample_incidence), "incidence")
            self.metric_manager.update_metric("d0", "estimate",
                                              hill_number_asymptotic(0, self.reference_sample_abundance,
                                                                     self.abundance_sample_size), "abundance")
            self.metric_manager.update_metric("d0", "estimate",
                                              hill_number_asymptotic(0, self.reference_sample_incidence,
                                                                     self.incidence_sample_size, abundance=False),
                                              "incidence")

        if self.include_d1:
            self.metric_manager.update_metric("d1", "sample", entropy_exp(self.reference_sample_abundance), "abundance")
            self.metric_manager.update_metric("d1", "sample", entropy_exp(self.reference_sample_incidence), "incidence")
            self.metric_manager.update_metric("d1", "estimate",
                                              hill_number_asymptotic(1, self.reference_sample_abundance,
                                                                     self.abundance_sample_size), "abundance")
            self.metric_manager.update_metric("d1", "estimate",
                                              hill_number_asymptotic(1, self.reference_sample_incidence,
                                                                     self.incidence_sample_size, abundance=False),
                                              "incidence")

        if self.include_d2:
            self.metric_manager.update_metric("d2", "sample", simpson_diversity(self.reference_sample_abundance),
                                              "abundance")
            self.metric_manager.update_metric("d2", "sample", simpson_diversity(self.reference_sample_incidence),
                                              "incidence")
            self.metric_manager.update_metric("d2", "estimate",
                                              hill_number_asymptotic(2, self.reference_sample_abundance,
                                                                     self.abundance_sample_size), "abundance")
            self.metric_manager.update_metric("d2", "estimate",
                                              hill_number_asymptotic(2, self.reference_sample_incidence,
                                                                     self.incidence_sample_size, abundance=False),
                                              "incidence")

        # Update completeness profiles
        if self.include_c0:
            self.metric_manager.update_metric("c0", "sample", completeness(self.reference_sample_abundance),
                                              "abundance")
            self.metric_manager.update_metric("c0", "sample", completeness(self.reference_sample_incidence),
                                              "incidence")

        if self.include_c1:
            self.metric_manager.update_metric("c1", "sample",
                                              coverage(self.reference_sample_abundance, self.abundance_sample_size),
                                              "abundance")
            self.metric_manager.update_metric("c1", "sample",
                                              coverage(self.reference_sample_incidence, self.incidence_sample_size),
                                              "incidence")

        # Update estimated sampling effort for target completeness
        for l in self.l_n:
            self.metric_manager.update_metric(f"l_{l}", "estimate",
                                              sampling_effort_abundance(l, self.reference_sample_abundance,
                                                                        self.abundance_sample_size), "abundance")
            self.metric_manager.update_metric(f"l_{l}", "estimate",
                                              sampling_effort_incidence(l, self.reference_sample_incidence,
                                                                        self.incidence_sample_size), "incidence")

    def print_metrics(self) -> None:
        """
        Print the diversity and completeness profile of the current observations.
        """
        metrics = self.metric_manager.get_metrics()
        print("### SAMPLE STATS ###")
        for metric, data in metrics.items():
            print(f"{metric.upper()} Metrics:")
            for key, values in data.items():
                print(f"  {key}: {values}")
            print()

    def to_dataframe(self) -> pd.DataFrame:
        """
        Return the diversity and completeness profile of the current observations as a DataFrame.

        :returns: A DataFrame view of the Diversity and Completeness Profile.
        """
        metrics = self.metric_manager.get_metrics()
        return pd.DataFrame.from_dict({(i, j): metrics[i][j]
                                       for i in metrics.keys()
                                       for j in metrics[i].keys()})
