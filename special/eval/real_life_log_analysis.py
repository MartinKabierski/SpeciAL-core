import species_estimator
import species_retrieval
from functools import partial
import pm4py
import pandas as pd

from src.visualizations.plots import plot_rank_abundance


#plt.style.use('seaborn-v0_8-white')
#plt.style.use('seaborn-v0_8-ticks')

def profile_log(log_path, name):
    log = pm4py.read_xes(log_path, return_legacy_log_object=True)
    print("##### "+name+" #####")
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
    print("Profiling log")
    for est_id,est in estimators.items():
        print(name,est_id)
        est.apply(log)
        metrics_stats[est_id]={}

    print("Saving metrics in csv file")

    for est_id, est in estimators.items():
        print()
        print(name,est_id)
        for metric in species_estimator.METRICS:
            metrics_stats[est_id][metric] = getattr(est, metric)[-1]
            print(metric +": "+str(getattr(est, metric)))

        plot_rank_abundance(est, str(name) + " - " + str(est_id))
        print("Total Abundances: " + str(est.total_number_species_abundances))
        print("Total Incidences: " + str(est.total_number_species_incidences))
        print("Degree of Aggregation: " + str(est.degree_spatial_aggregation))
        print("Observations Incidence: " + str(est.number_observations_incidence))
        print("Species Counts: ")
        [print(x, end=" ") for x in est.reference_sample_incidence.values()]
        print()
    df_stats = pd.DataFrame.from_dict(metrics_stats)
    df_stats.to_csv("results/" + str(name) + "_metrics.csv")

#TODO rarefaction based on bootstrapping! Becomes prohibitively slow for q>=1
    # print("Building Rarefaction and Extrapolation Curves")
    # for est_id, est in estimators.items():
    #     q0, q1, q2 = rarefy_extrapolate_all(est, abundance_data=False, data_points=20)
    #     print(q0[0])
    #     print(q0[1])
    #
    #     plt.rcParams['figure.figsize'] = [6, 5]
    #     plt.rcParams['xtick.labelsize'] = 14
    #     plt.rcParams['ytick.labelsize'] = 14
    #
    #     plt.plot(q0[1][:22], q0[0][:22], label="q=0", color='#F8766D')
    #     plt.plot(q0[1][21:], q0[0][21:], linestyle="--", color='#F8766D')
    #
    #     plt.plot(q1[1][:22], q1[0][:22], label="q=1", color='#00BA38')
    #     plt.plot(q1[1][21:], q1[0][21:], linestyle="--", color='#00BA38')
    #
    #     plt.plot(q2[1][:22], q2[0][:22], label="q=2", color='#619CFF')
    #     plt.plot(q2[1][21:], q2[0][21:], linestyle="--", color='#619CFF')
    #
    #     plt.title(name, fontsize=22)
    #     plt.xlabel("Sample Size", fontsize=18)
    #     plt.ylabel("Diversity", fontsize=18)
    #     plt.xticks([0, est.number_observations_incidence, 2 * est.number_observations_incidence],[0, est.number_observations_incidence, 2 * est.number_observations_incidence])
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig("figures/real_eval_curves_" + name + "_" + est_id +  "_.pdf", format="pdf")
    #
    #     plt.show()
    #save all metrics to csv


profile_log("../../logs/BPI_Challenge_2012.xes", "real_eval_BPI 2012")

profile_log("../../logs/Road_Traffic_Fines_Management_Process.xes", "real_eval_Road Traffic Fines")

profile_log("../../logs/Sepsis_Cases_-_Event_Log.xes", "real_eval_Sepsis Cases")

profile_log("../../logs/BPI_Challenge_2018.xes", "real_eval_BPI 2018")

profile_log("../../logs/BPI_Challenge_2019.xes", "real_eval_BPI 2019")