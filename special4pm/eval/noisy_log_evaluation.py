import copy
import math

import species_estimator
import species_retrieval
from functools import partial

from plots import plot_rank_abundances_multiple_estimators
from species_estimator import SpeciesEstimator
from src.simulation import simulate_model
import pm4py
import random

#plt.style.use('seaborn-v0_8-ticks')


def profile_logs(log_dict):
    estimations = {"1-gram": [], "2-gram": [], "3-gram": [], "4-gram": [], "5-gram": [], "trace_variants": []}
    metrics = species_estimator.METRICS

    for log_id, log in log_dict.items():
        estimators = \
            {
                "1-gram": SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=1),
                                           quantify_all=True),
                "2-gram": SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=2),
                                           quantify_all=True),
                "3-gram": SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=3),
                                           quantify_all=True),
                "4-gram": SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=4),
                                           quantify_all=True),
                "5-gram": SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=5),
                                           quantify_all=True),
                "trace_variants": SpeciesEstimator(species_retrieval.retrieve_species_trace_variant, quantify_all=True)
            }

        for key, est in estimators.items():
            print("Profiling log "+log_id+" "+key)
            est.apply(log)
            estimations[key].append(est)

    for est_id, estimation in estimations.items():
        print(est_id)
        for log_id, est in zip(log_dict.keys(), estimation):
            print(log_id)
            for metric in metrics:
                print(metric +": "+str(getattr(est, metric)))

            print("Total Abundances: " + str(est.total_number_species_abundances))
            print("Total Incidences: " + str(est.total_number_species_incidences))
            print("Degree of Aggregation: " + str(est.degree_spatial_aggregation))
            print("Observations Incidence: " + str(est.number_observations_incidence))
            print("Species Counts: ")
            [print(x, end=" ") for x in est.reference_sample_incidence.values()]
            print()

        plot_rank_abundances_multiple_estimators(log_dict.keys(), estimation, est_id)






net, im, fm = pm4py.read_pnml("../../models/base_net.pnml")
#pm4py.view_petri_net(net, im, fm)

simulated_log = simulate_model(net, im, fm, 1, 1000)

events = []
for i,tr in enumerate(simulated_log[0]):
    for j,ev in enumerate(tr):
        events.append((i,j))
print("Total Number of events: "+ str(len(events)),"1%: "+str(len(events)/100),"0.1%: "+str(len(events)/1000),"0.01%: "+str(len(events)/10000))

#log consisting of 0.01% noise
log_a = copy.deepcopy(simulated_log[0])
evs = random.sample(events,math.ceil(len(events)/10000))
for ev in evs:
    log_a[ev[0]][ev[1]]["concept:name"] = str(random.getrandbits(128))

#log consisting of 0.1% noise
log_b = copy.deepcopy(simulated_log[0])
evs = random.sample(events,math.ceil(len(events)/1000))
for ev in evs:
    log_b[ev[0]][ev[1]]["concept:name"] = str(random.getrandbits(128))

#log consisting of 1% noise
log_c = copy.deepcopy(simulated_log[0])
evs = random.sample(events,math.ceil(len(events)/100))
for ev in evs:
    log_c[ev[0]][ev[1]]["concept:name"] = str(random.getrandbits(128))

profile_logs({"Original":simulated_log[0],"0.01%":log_a,"0.1%":log_b,"1%":log_c})