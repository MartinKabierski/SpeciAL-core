from functools import partial

import pm4py
from matplotlib import pyplot as plt

from special.estimation import species_estimator, species_retrieval
from special.simulation.simulation import simulate_model
from special.visualization.visualization import plot_rank_abundance, plot_diversity_series, plot_diversity_series_all, \
    plot_diversity_sample_vs_estimate, plot_diversity_profile, plot_completeness_profile, plot_expected_sampling_effort


def profile_log(log, name):
    """Estimates completeness profiles and species richness for different species retrieval functions on the provided
    log, printing obtained metrics over all observations in the log each 100 traces"""
    print("Profiling log " + name)

    estimator = species_estimator.SpeciesEstimator(d0=True, d1=True, d2=True, c0=True, c1=True, step_size=1000)
    estimator.register("1-gram", partial(species_retrieval.retrieve_species_n_gram, n=1))
    estimator.register("2-gram", partial(species_retrieval.retrieve_species_n_gram, n=2))
    estimator.register("3-gram", partial(species_retrieval.retrieve_species_n_gram, n=3))
    estimator.register("4-gram", partial(species_retrieval.retrieve_species_n_gram, n=4))
    estimator.register("5-gram", partial(species_retrieval.retrieve_species_n_gram, n=5))

    estimator.apply(log)
    estimator.print_metrics()
    estimator.to_dataFrame().to_csv("test.csv", index=False)

    #TODO only show these and remove file name from function, also unifiy call signatures
    plot_rank_abundance(estimator, "2-gram", "test_abundance.pdf", abundance=True)
    plot_rank_abundance(estimator, "2-gram", "test_incidence.pdf", abundance=False)
    plot_diversity_sample_vs_estimate(estimator, "2-gram", ["d0","d1","d2"] , "example.pdf",False)
    plot_diversity_series(estimator, "2-gram", "d0", "d0_example.pdf",False)
    plot_diversity_series(estimator, "2-gram", "d1", "d1_example.pdf",False)
    plot_diversity_series(estimator, "2-gram", "d2", "d2_example.pdf",False)
    plot_diversity_series_all(estimator, "2-gram", ["d0","d1","d2"], "all_example.pdf",False)
    plot_diversity_profile(estimator, "2-gram", "example.pdf",False)
    plot_completeness_profile(estimator, "2-gram", "completeness_example.pdf",False)
    plot_diversity_profile(estimator, "2-gram", "example.pdf",False)
    plot_expected_sampling_effort(estimator, "2-gram", "example.pdf",False)




#Estimating species richness of an event log
PATH_TO_XES = "logs/BPI_Challenge_2019.xes"
log = pm4py.read_xes(PATH_TO_XES)

profile_log(log, "example_log")

#Estimating species richness of a process model by simulating 1 log of size 500
#NO_LOGS = 1
#LOG_SIZE = 5000
#PATH_TO_XES = "logs/BPI_Challenge_2012.xes"
#log = pm4py.read_xes(PATH_TO_XES)
#mod, i, f = pm4py.discover_petri_net_inductive(log, noise_threshold=0.2)
#log_ind = simulate_model(mod, i, f, NO_LOGS, LOG_SIZE)
#profile_log(log_ind[0], "example_model")
