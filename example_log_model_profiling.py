from functools import partial

import pm4py

from SpeciAL.estimation import species_estimator, species_retrieval
from SpeciAL.simulation.simulation import simulate_model


def profile_log(log, name):
    """Estimates completeness profiles and species richness for different species retrieval functions on the provided
    log, printing obtained metrics over all observations in the log each 100 traces"""
    print("Profiling log " + name)

    #used species retrieval functions
    estimators = \
        {
            #Different Species Retrieval Functions, metrics are updated each 100 traces
            "1-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=1),
                                                         step_size=100),
            "2-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=2),
                                                         step_size=100),
            "3-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=3),
                                                         step_size=100),
            "4-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=4),
                                                         step_size=100),
            "5-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=5),
                                                         step_size=100),
            "trace_variants": species_estimator.SpeciesEstimator(species_retrieval.retrieve_species_trace_variant,
                                                                 step_size=100),


            #Directly Follows Relations, metrics are calculated once after obsering the whole log
            "2-gram_complete_log": species_estimator.SpeciesEstimator(
                partial(species_retrieval.retrieve_species_n_gram, n=2)),


            #Directly Follows Relations, metrics are calculated once after obsering the whole log, explicit selection of profiles and metrics to compute
            "2-gram_all_metrics": species_estimator.SpeciesEstimator(
                partial(species_retrieval.retrieve_species_n_gram, n=2), d0=True, d1=True, d2=True, c0=True, c1=True,
                l_n=[0.9, 0.95, 0.99]),
        }

    #profile log using each of the retrieval functions
    for est_id, est in estimators.items():
        est.profile_log(log)

        #print obtained metrics
        print("###### " + name, est_id + " ######")
        est.print_metrics()
        print()


#Estimating species richness of an event log
PATH_TO_XES = "/Users/christianimenkamp/Documents/Data-Repository/Community/sepsis/Sepsis Cases - Event Log.xes"
log = pm4py.read_xes(PATH_TO_XES)  #, return_legacy_log_object=True)

profile_log(log, "example_log")

#Estimating species richness of a process model by simulating 1 log of size 500
NO_LOGS = 1
LOG_SIZE = 5000
PATH_TO_XES = "/Users/christianimenkamp/Documents/Data-Repository/Community/sepsis/Sepsis Cases - Event Log.xes"
log = pm4py.read_xes(PATH_TO_XES)
mod, i, f = pm4py.discover_petri_net_inductive(log, noise_threshold=0.2)
log_ind = simulate_model(mod, i, f, NO_LOGS, LOG_SIZE)
profile_log(log_ind[0], "example_model")
