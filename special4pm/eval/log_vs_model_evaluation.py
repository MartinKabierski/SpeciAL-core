from functools import partial

import pandas as pd
import pm4py

import species_estimator
import species_retrieval
from src.simulation import simulate_model


def profile_log(log, name):
    print("Profiling log " + name)

    estimators = \
        {
            "1-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=1),
                                                         quantify_all=True, hill="q0"),
            "2-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=2),
                                                         quantify_all=True, hill="q0"),
            "3-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=3),
                                                         quantify_all=True, hill="q0"),
            "4-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=4),
                                                         quantify_all=True, hill="q0"),
            "5-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=5),
                                                         quantify_all=True, hill="q0"),
            "trace_variants": species_estimator.SpeciesEstimator(species_retrieval.retrieve_species_trace_variant,
                                                                 quantify_all=True, hill="q0")
        }
    metrics_stats = {}
    for est_id, est in estimators.items():
        print(name, est_id)
        est.apply(log)
        metrics_stats[est_id] = {}

    print("Saving metrics in csv file")
    for est_id, est in estimators.items():
        print(name, est_id)
        for metric in species_estimator.METRICS_q1:
            metrics_stats[est_id][metric] = getattr(est, metric)[-1]
            print(metric + ": " + str(getattr(est, metric)))
        #plot_rank_abundance(est, str(name) + " - " + str(est_id))
        print()
        print("Total Abundances: " + str(est.total_number_species_abundances))
        print("Total Incidences: " + str(est.total_number_species_incidences))
        print("Degree of Aggregation: " + str(est.degree_spatial_aggregation))
        print("Observations Incidence: " + str(est.number_observations_incidence))
        #print("Species Counts: ")
        #[print(x, end=" ") for x in est.reference_sample_incidence.values()]
        print()
    df_stats = pd.DataFrame.from_dict(metrics_stats)
    df_stats.to_csv("results/" + str(name) + "_metrics.csv")



log = pm4py.read_xes("../../logs/BPI_Challenge_2012.xes", return_legacy_log_object=True)

print("###1###")
mod, i, f = pm4py.discover_petri_net_alpha(log)
print(i,f)
pm4py.save_vis_petri_net(mod, i, f, "BPI_2012_ALPHA_pm4py.pdf")
log_d_alpha = simulate_model(mod, i, f, 1, 5000)
print(log_d_alpha[0])
profile_log(log_d_alpha[0], "log_vs_model_evalB_PI_2012_alpha_pm4py")
print()

#print("###1.5###")
#mod, i, f = pm4py.discover_petri_net_alpha_plus(log)
#print(i,f)
#pm4py.save_vis_petri_net(mod, i, f, "ALPHA_PLUS_pm4py.pdf")
#log_d_alpha = simulate_model(mod, i, f, 1, 5000)
#profile_log(log_d_alpha[0], "log_vs_model_eval_alpha_plus_pm4py")

print("###2###")
mod2, i, f = pm4py.discover_petri_net_inductive(log, noise_threshold=0.2)
print(i,f)
pm4py.save_vis_petri_net(mod2, i, f, "BPI_2012_INDUCTIVE_infrequent_0.2_pm4py.pdf")
log_d_ind = simulate_model(mod2, i, f, 1, 5000)
profile_log(log_d_ind[0], "log_vs_model_eval_BPI_2012_inductive_infrequent_0.2_pm4py")

print("###2.05###")
mod2, i, f = pm4py.discover_petri_net_inductive(log, noise_threshold=0.4)
print(i,f)
pm4py.save_vis_petri_net(mod2, i, f, "BPI_2012_INDUCTIVE_infrequent_0.4_pm4py.pdf")
log_d_ind = simulate_model(mod2, i, f, 1, 5000)
profile_log(log_d_ind[0], "log_vs_model_eval_BPI_2012_inductive_infrequent_0.4_pm4py")

print("###2.1###")
mod2, i, f = pm4py.discover_petri_net_inductive(log, noise_threshold=0.6)
print(i,f)
pm4py.save_vis_petri_net(mod2, i, f, "BPI_2012_INDUCTIVE_infrequent_0.6_pm4py.pdf")
log_d_ind = simulate_model(mod2, i, f, 1, 5000)
profile_log(log_d_ind[0], "log_vs_model_eval_BPI_2012_inductive_infrequent_0.6_pm4py")


print("###2.2###")
mod2, i, f = pm4py.discover_petri_net_inductive(log, noise_threshold=0.8)
print(i,f)
pm4py.save_vis_petri_net(mod2, i, f, "BPI_2012_INDUCTIVE_infrequent_0.8_pm4py.pdf")
log_d_ind = simulate_model(mod2, i, f, 1, 5000)
profile_log(log_d_ind[0], "log_vs_model_eval_BPI_2012_inductive_infrequent_0.8_pm4py")

print("###3.1###")
mod3, i, f = pm4py.discover_petri_net_ilp(log, alpha=0.8)
print(i,f)
pm4py.save_vis_petri_net(mod3, i, f, "BPI_2012_ILP_0.8_pm4py.pdf")
log_d_ind2 = simulate_model(mod3, i, f, 1, 5000)
profile_log(log_d_ind2[0], "log_vs_model_eval_BPI_2012_ilp_0.8_pm4py")

print("###3.1###")
mod3, i, f = pm4py.discover_petri_net_ilp(log, alpha=0.7)
print(i,f)
pm4py.save_vis_petri_net(mod3, i, f, "BPI_2012_ILP_0.7_pm4py.pdf")
log_d_ind2 = simulate_model(mod3, i, f, 1, 5000)
profile_log(log_d_ind2[0], "log_vs_model_eval_BPI_2012_ilp_0.7_pm4py")

print("###3.2###")
mod3, i, f = pm4py.discover_petri_net_ilp(log, alpha=0.6)
print(i,f)
pm4py.save_vis_petri_net(mod3, i, f, "BPI_2012_ILP_0.6_pm4py.pdf")
log_d_ind2 = simulate_model(mod3, i, f, 1, 5000)
profile_log(log_d_ind2[0], "log_vs_model_eval_BPI_2012_ilp_0.6_pm4py")

print("###3.3###")
mod3, i, f = pm4py.discover_petri_net_ilp(log, alpha=0.4)
print(i,f)
pm4py.save_vis_petri_net(mod3, i, f, "BPI_2012_ILP_0.4_pm4py.pdf")
log_d_ind2 = simulate_model(mod3, i, f, 1, 5000)
profile_log(log_d_ind2[0], "log_vs_model_eval_BPI_2012_ilp_0.4_pm4py")



print("###3.4###")
mod3, i, f = pm4py.discover_petri_net_ilp(log, alpha=0.2)
print(i,f)
pm4py.save_vis_petri_net(mod3, i, f, "BPI_2012_ILP_0.2_pm4py.pdf")
log_d_ind2 = simulate_model(mod3, i, f, 1, 5000)
profile_log(log_d_ind2[0], "log_vs_model_eval_BPI_2012_ilp_0.2_pm4py")

print("###3.5###")
mod3, i, f = pm4py.discover_petri_net_ilp(log, alpha=1.0)
print(i,f)
pm4py.save_vis_petri_net(mod3, i, f, "BPI_2012_ILP_1.0_pm4py.pdf")
log_d_ind2 = simulate_model(mod3, i, f, 1, 5000)
profile_log(log_d_ind2[0], "log_vs_model_eval_BPI_2012_ilp_1.0_pm4py")


print("###4###")
mod3, i, f = pm4py.discover_petri_net_inductive(log)
print(i,f)
pm4py.save_vis_petri_net(mod3, i, f, "BPI_2012_INDUCTIVE_0.0_pm4py.pdf")
log_d_ind2 = simulate_model(mod3, i, f, 1, 5000)
profile_log(log_d_ind2[0], "log_vs_model_eval_BPI_2012_inductive_pm4py")

print("###5###")
net, im, fm = pm4py.read_pnml("../../models/bpi_2012_alpha.pnml")
print(im,fm)
pm4py.save_vis_petri_net(net, im, fm, "BPI_2012_ALPHA.pdf")
log_alpha = simulate_model(net, im, fm, 1, 5000)
profile_log(log_alpha[0], "log_vs_model_eval_BPI_2012_alpha")

print("###6###")
net, im, fm = pm4py.read_pnml("../../models/bpi_2012_im.pnml")
print(im,fm)
pm4py.save_vis_petri_net(net, im, fm,"BPI_2012_INDUCTIVE.pdf")
log_im = simulate_model(net, im, fm, 1, 5000)
profile_log(log_im[0], "log_vs_model_eval_BPI_2012_inductive")

print("###7###")
net, im, fm = pm4py.read_pnml("../../models/bpi_2012_im_infrequent.pnml")
print(im,fm)
pm4py.save_vis_petri_net(net, im, fm, "BPI_2012_INDUCTIVE_infrequent.pdf")
log_im_infrequent = simulate_model(net, im, fm, 1, 5000)
profile_log(log_im_infrequent[0], "log_vs_model_eval_BPI_2012_inductive_infrequent")

#print("###8###")
#net, im, fm = pm4py.read_pnml("models/bpi_2012_im_incomplete.pnml")
#pm4py.save_vis_petri_net(net, im, fm, "INDUCTIVE_incomplete.pdf")
#log_im_incomplete = simulate_model(net, im, fm, 1, 5000)
#profile_log(log_im_incomplete[0], "log_vs_model_eval_inductive_incomplete")

