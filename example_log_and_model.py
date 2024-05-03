from functools import partial

import pandas as pd
import pm4py

import species_estimator
import species_retrieval
from simulation import simulate_model


def profile_log(log, name):
"""Estimates completeness profiles and species richness for different species retrieval functions on the provided log, printing obtained metrics over all observations in the log"""
    print("Profiling log " + name)

    #used species retrieval functions
    estimators = \
        {
            "1-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=1), step_size=10,hill="q0")),
            "2-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=2), step_size=10,hill="q0"),
            "3-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=3), step_size=10,hill="q0"),
            "4-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=4), step_size=10,hill="q0"),
            "5-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=5), step_size=10,hill="q0"),
            "trace_variants": species_estimator.SpeciesEstimator(species_retrieval.retrieve_species_trace_variant, step_size=10,hill="q0"),
        }
        
    #profile log using each of the retrieval functions
    for est_id, est in estimators.items():
        print(name, est_id)
        est.profile_log(log)

    #for each retrieval function, print all obtained metrics
    for est_id, est in estimators.items():
        print("##### "+ est_id + " #####")
        for metric in species_estimator.METRICS_q1:
            print(metric + ": " + str(getattr(est, metric)))
        print()


#Estimating species richness of an event log
PATH_TO_XES = "./logs/BPI_Challenge_2012.xes"
log = pm4py.read_xes(PATH_TO_XES, return_legacy_log_object=True)
profile_log(log, "example_log_richness")


#Estimating species richness of a process model generating 1 log of size 500
NO_LOGS = 1
LOG_SIZE = 5000
PATH_TO_XES = "./logs/BPI_Challenge_2012.xes"
log = pm4py.read_xes(PATH_TO_XES, return_legacy_log_object=True)
mod, i, f = pm4py.discover_petri_net_inductive(log, noise_threshold=0.2)
log_ind = simulate_model(mod, i, f, NO_LOGS, LOG_SIZE)
profile_log(log_ind[0], "example_model_richness")

