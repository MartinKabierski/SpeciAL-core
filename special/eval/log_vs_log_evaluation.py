import copy
from functools import partial

import pandas as pd
import pm4py
from pm4py.objects.log.obj import EventLog

import species_estimator
import species_retrieval
from plots import plot_rank_abundance


def profile_log(log, name):
    estimators = \
        {
            "1-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=1),
                                                         quantify_all=True),
            "2-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=2),
                                                         quantify_all=True),
            "3-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=3),
                                                         quantify_all=True),
            "4-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=4),
                                                         quantify_all=True),
            "5-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=5),
                                                         quantify_all=True),
            "trace_variants": species_estimator.SpeciesEstimator(species_retrieval.retrieve_species_trace_variant,
                                                                 quantify_all=True),
            "est_act_1" : species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_timed_activity, interval_size=1), quantify_all=True),
            "est_act_5" : species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_timed_activity, interval_size=5), quantify_all=True),
            "est_act_30" : species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_timed_activity, interval_size=30), quantify_all=True),
            "est_act_exp" : species_estimator.SpeciesEstimator(partial(
                species_retrieval.retrieve_timed_activity_exponential), quantify_all=True)
        }
    metrics_stats = {}
    print("Profiling log " + name)
    for est_id, est in estimators.items():
        print(est_id)
        est.apply(log)
        metrics_stats[est_id] = {}

    print("Saving metrics in csv file")
    for est_id, est in estimators.items():
        for metric in species_estimator.METRICS:
            #print(metric)
            metrics_stats[est_id][metric] = getattr(est, metric)[0]
            print(metric + ": " + str(getattr(est, metric)))
        print(est_id)
        plot_rank_abundance(est, str(name) + " - " + str(est_id))
        print()
        print("Total Abundances: " + str(est.total_number_species_abundances))
        print("Total Incidences: " + str(est.total_number_species_incidences))
        print("Degree of Aggregation: " + str(est.degree_spatial_aggregation))
        print("Observations Incidence: " + str(est.number_observations_incidence))
        print("Species Counts: ")
        [print(x, end=" ") for x in est.reference_sample_incidence.values()]
        print()
    df_stats = pd.DataFrame.from_dict(metrics_stats)
    df_stats.to_csv("results/" + str(name) + "_metrics.csv")



log = pm4py.read_xes("../../logs/Sepsis_Cases_-_Event_Log.xes", return_legacy_log_object=True)

log_pre_admission = copy.deepcopy(log)
log_post_admission = copy.deepcopy(log)

log_young = EventLog(attributes=log.attributes, extensions=log.extensions, omni_present=log.omni_present,
                                classifiers=log.classifiers, properties=log.properties)
log_old = EventLog(attributes=log.attributes, extensions=log.extensions, omni_present=log.omni_present,
                                classifiers=log.classifiers, properties=log.properties)


for t in range(0, len(log)):
    #print(log[t])
    if log[t][0]['Age']>=60:
        log_old.append(log[t])
    else:
        log_young.append(log[t])
    for idx, e in enumerate(log[t]):
        if "Admission" in e["concept:name"]:
            if len(log_pre_admission[t][:idx + 1])>0:
                log_pre_admission[t] = log_pre_admission[t][:idx + 1]

            if len(log_post_admission[t][idx + 1:]) > 0:
                log_post_admission[t] = log_post_admission[t][idx + 1:]
            #print([e["concept:name"] for e in log_a[t]])
            #print([e["concept:name"] for e in log_pre_admission[t]])
            #print([e["concept:name"] for e in log_post_admission[t]])
            #print()
            break


profile_log(log_pre_admission, "log_vs_log_eval_sepsis_cases_pre_admission")
profile_log(log_post_admission, "log_vs_log_eval_sepsis_cases_post_admission")
profile_log(log_young, "log_vs_log_eval_sepsis_cases_age_less_60")
profile_log(log_old, "log_vs_log_eval_sepsis_cases_age_geq_60")
