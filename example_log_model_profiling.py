from functools import partial

import pm4py

from special.estimation import species_estimator, species_retrieval
from special.simulation.simulation import simulate_model


def profile_log(log, name):
    """Estimates completeness profiles and species richness for different species retrieval functions on the provided
    log, printing obtained metrics over all observations in the log each 100 traces"""
    print("Profiling log " + name)

    estimator = species_estimator.SpeciesEstimator(d0=True, d1=True, d2=True, c0=True, c1=True, step_size=100)
    estimator.register("1-gram", partial(species_retrieval.retrieve_species_n_gram, n=1))
    estimator.register("2-gram", partial(species_retrieval.retrieve_species_n_gram, n=2))
    estimator.register("3-gram", partial(species_retrieval.retrieve_species_n_gram, n=3))
    estimator.register("4-gram", partial(species_retrieval.retrieve_species_n_gram, n=4))
    estimator.register("5-gram", partial(species_retrieval.retrieve_species_n_gram, n=5))

    estimator.apply(log)
    estimator.print_metrics()
    estimator.to_dataFrame().to_csv("test.csv", index=False)


#Estimating species richness of an event log
PATH_TO_XES = "logs/BPI_Challenge_2012.xes"
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
