from functools import partial

from pandas import DataFrame

import species_retrieval
import species_estimator

from plots import plot_all_stats
from src.simulation import simulate_model
import pm4py
import copy

STEP_SIZE = 10
#plt.style.use('seaborn-v0_8-ticks')


def profile_log(log, name, tick_locations=[20,40,60,80,100], tick_names=[200,400,600,800,1000]):
    estimators = \
        {
            "1-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=1), step_size=10),
            "2-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=2), step_size=10),
            "3-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=3), step_size=10),
            #"4-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=4), step_size=10),
            #"5-gram": species_estimator.SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=5), step_size=10),
            #"trace_variants": species_estimator.SpeciesEstimator(species_retrieval.retrieve_species_trace_variant, step_size=10),
        }

    df_stats = {}

    for key, est in estimators.items():
        print("Profiling Log "+str(name)+" "+str(key))
        est.apply(log)
        df_stats[key]={}

    for est_id, est in estimators.items():
        print(name, est_id)
        for metric in species_estimator.METRICS:
            print(metric +": "+str(getattr(est, metric)[-1]))
            df_stats[est_id][metric]=getattr(est, metric)[-1]

        print("Total Abundances: " + str(est.total_number_species_abundances))
        print("Total Incidences: " + str(est.total_number_species_incidences))
        print("Degree of Aggregation: " + str(est.degree_spatial_aggregation))
        print("Observations Incidence: " + str(est.number_observations_incidence))
        print("Species Counts: ")
        [print(x, end=" ") for x in est.reference_sample_incidence.values()]
        print()
        print()
        plot_all_stats([est], "drift_eval_" + name + "_" + est_id, tick_locations, tick_names)

    df = DataFrame.from_dict(df_stats)
    df.to_csv("results/drift_eval_"+name+"_metrics.csv")



# plot functions
# def plot_rank_abundance(estimation, name):
#     plt.rcParams['figure.figsize'] = [8, 3]
#     plt.rcParams['xtick.labelsize'] = 18
#     plt.rcParams['ytick.labelsize'] = 18
#
#     no_species = len(estimation.reference_sample_abundance.values())
#     plt.plot(sorted(list(estimation.reference_sample_abundance.values()), reverse=True))
#     plt.fill_between(np.linspace(0, no_species, no_species, endpoint=False),
#                      sorted(list(estimation.reference_sample_abundance.values()), reverse=True),
#                      [0 for _ in range(no_species)], alpha=0.5)
#     plt.title("Species Occurences (Abundance)", fontsize=26)
#     plt.xlabel("Species Rank", fontsize=22)
#     plt.ylabel("Occurences", fontsize=22)
#     plt.xticks([0, no_species - 1], [0, no_species - 1])
#     plt.yticks([0, max(estimation.reference_sample_abundance.values()) - 1],
#                [0, max(estimation.reference_sample_abundance.values()) - 1])
#
#     plt.tight_layout()
#     plt.savefig("figures/drift_eval_" + name + "_curve_abundance.pdf", format="pdf")
#     plt.show()
#
#     no_species = len(estimation.reference_sample_incidence.values())
#     plt.plot(sorted(list(estimation.reference_sample_incidence.values()), reverse=True))
#     plt.fill_between(np.linspace(0, no_species, no_species, endpoint=False),
#                      sorted(list(estimation.reference_sample_incidence.values()), reverse=True),
#                      [0 for _ in range(no_species)], alpha=0.5)
#     plt.title("Species Occurence (Incidence)", fontsize=26)
#     plt.xlabel("Species Rank", fontsize=22)
#     plt.ylabel("Occurences", fontsize=22)
#     plt.xticks([0, no_species - 1], [0, no_species - 1])
#     plt.yticks([0, max(estimation.reference_sample_incidence.values()) - 1],
#                [0, max(estimation.reference_sample_incidence.values()) - 1])
#
#     plt.tight_layout()
#     plt.savefig("figures/drift_eval_" + name + "_curve_incidence.pdf", format="pdf")
#     plt.show()


# def plot_all_stats(estimation, tick_locations, tick_names, name):
#     labels = estimation.observation_ids
#     no_data_points = len(labels)
#     # no_repetitions = len(estimation)
#     no_observations_abundance = estimation.number_observations_abundance
#     no_observations_incidence = estimation.number_observations_incidence
#
#     obs_stats = estimation.observed_species_count
#     chao2_abundance_stats = estimation.chao2_abundance
#     chao2_incidence_stats = estimation.chao2_incidence
#
#     D1_abundance_sample_stats = estimation.D1_sample_abundance
#     D1_incidence_sample_stats = estimation.D1_sample_incidence
#     D1_abundance_estimated_stats = estimation.D1_estimated_abundance
#     D1_incidence_estimated_stats = estimation.D1_estimated_incidence
#
#     D2_abundance_sample_stats = estimation.D2_sample_abundance
#     D2_incidence_sample_stats = estimation.D2_sample_incidence
#     D2_abundance_estimated_stats = estimation.D2_estimated_abundance
#     D2_incidence_estimated_stats = estimation.D2_estimated_incidence
#
#     coverage_abundance_stats = estimation.coverage_abundance
#     coverage_incidence_stats = estimation.coverage_incidence
#
#     completeness_abundance_stats = estimation.completeness_abundance
#     completeness_incidence_stats = estimation.completeness_incidence
#
#     plt.rcParams['figure.figsize'] = [28, 5]
#     plt.rcParams['xtick.labelsize'] = 18
#     plt.rcParams['ytick.labelsize'] = 18
#
#     # abundance plots
#     # NOTE: we only collect stats after trace level, thus we plot abundance-based stats using incidence observation counts
#     # as number of observations instead of abundance observation counts
#
#     f, (ax1, ax2, ax3), = plt.subplots(nrows=1, ncols=3, sharey='all')
#     y_max = (max(max(obs_stats), max(chao2_abundance_stats), max(chao2_incidence_stats))) + 1
#
#     ax1.plot(obs_stats, label="Observed")
#     ax1.plot(chao2_abundance_stats, label="Estimated")
#     ax1.set_title("q=0 (Abundance)", fontsize=26)
#     ax1.set_xlabel("Sample Size", fontsize=22)
#     ax1.set_ylabel("Hill number", fontsize=22)
#     ax1.set_xticks(tick_locations)
#     ax1.set_xticklabels(tick_names)
#     ax1.tick_params(labelleft=True)
#     ax1.set_ylim(bottom=0, top=y_max)
#     ax1.legend()
#
#     ax1.figure.savefig('TEST.pdf')
#
#     ax2.plot(D1_abundance_sample_stats, label="Observed")
#     ax2.plot(D1_abundance_estimated_stats, label="Estimated")
#     ax2.set_title("q=1 (Abundance)", fontsize=26)
#     ax2.set_xlabel("Sample Size", fontsize=22)
#     ax2.set_ylabel("Hill number", fontsize=22)
#     ax2.set_xticks(tick_locations)
#     ax2.set_xticklabels(tick_names)
#     ax2.tick_params(labelleft=True)
#     ax2.set_ylim(bottom=0, top=y_max)
#     ax2.legend()
#
#     ax3.plot(D2_abundance_sample_stats, label="Observed")
#     ax3.plot(D2_abundance_estimated_stats, label="Estimated")
#     ax3.set_title("q=2 (Abundance)", fontsize=26)
#     ax3.set_xlabel("Sample Size", fontsize=22)
#     ax3.set_ylabel("Hill number", fontsize=22)
#     ax3.set_xticks(tick_locations)
#     ax3.set_xticklabels(tick_names)
#     ax3.tick_params(labelleft=True)
#     ax3.set_ylim(bottom=0, top=y_max)
#     ax3.legend()
#
#     plt.tight_layout()
#     f.savefig("figures/drift_eval_" + name + "_abundance.pdf", format="pdf")
#     plt.show()
#
#     f, (ax4, ax5, ax6) = plt.subplots(nrows=1, ncols=3, sharey='all')
#
#     ax4.plot(obs_stats, label="Observed")
#     ax4.plot(chao2_incidence_stats, label="Estimated")
#     ax4.set_title("q=0 (Incidence)", fontsize=26)
#     ax4.set_xlabel("Sample Size", fontsize=22)
#     ax4.set_ylabel("Hill number", fontsize=22)
#     ax4.set_xticks(tick_locations)
#     ax4.set_xticklabels(tick_names)
#     ax4.tick_params(labelleft=True)
#     ax4.set_ylim(bottom=0, top=y_max)
#     ax4.legend()
#
#     ax5.plot(D1_incidence_sample_stats, label="Observed")
#     ax5.plot(D1_incidence_estimated_stats, label="Estimated")
#     ax5.set_title("q=1 (Incidence)", fontsize=26)
#     ax5.set_xlabel("Sample Size", fontsize=22)
#     ax5.set_ylabel("Hill number", fontsize=22)
#     ax5.set_xticks(tick_locations)
#     ax5.set_xticklabels(tick_names)
#     ax5.tick_params(labelleft=True)
#     ax5.set_ylim(bottom=0, top=y_max)
#     ax5.legend()
#
#     ax6.plot(D2_incidence_sample_stats, label="Observed")
#     ax6.plot(D2_incidence_estimated_stats, label="Estimated")
#     ax6.set_title("q=2 (Incidence)", fontsize=26)
#     ax6.set_xlabel("Sample Size", fontsize=22)
#     ax6.set_ylabel("Hill number", fontsize=22)
#     ax6.set_xticks(tick_locations)
#     ax6.set_xticklabels(tick_names)
#     ax6.tick_params(labelleft=True)
#     ax6.set_ylim(bottom=0, top=y_max)
#     ax6.legend()
#
#     plt.tight_layout()
#     f.savefig("figures/drift_eval_" + name + "_incidence.pdf", format="pdf")
#     plt.show()
#
#     plt.rcParams['figure.figsize'] = [9, 4]
#
#
#     plt.plot(completeness_abundance_stats, label="Completeness")
#     plt.plot(coverage_abundance_stats, label="Coverage")
#     plt.title("Completeness & Coverage", fontsize=26)
#     plt.xlabel("Sample Size", fontsize=22)
#     plt.ylabel("Completeness/Coverage", fontsize=22)
#     plt.xticks(tick_locations, tick_names)
#     plt.ylim(bottom=0)
#     plt.legend()
#
#     plt.tight_layout()
#     f.savefig("figures/drift_eval_" + name + "_coverage_completeness_abundance.pdf", format="pdf")
#     plt.show()
#
#     plt.plot(completeness_incidence_stats, label="Completeness")
#     plt.plot(coverage_incidence_stats, label="Coverage")
#     plt.title("Completeness & Coverage", fontsize=26)
#     plt.xlabel("Sample Size", fontsize=22)
#     plt.ylabel("Completeness/Coverage", fontsize=22)
#     plt.xticks(tick_locations, tick_names)
#     plt.ylim(bottom=0)
#
#     plt.legend()
#
#     plt.tight_layout()
#     f.savefig("figures/drift_eval_" + name + "_coverage_completeness_incidence.pdf", format="pdf")
#     plt.show()


BASE_LENGTH = 1000
MODIFIED_LENGTH = 1000

net, im, fm = pm4py.read_pnml("../../models/base_net.pnml")
#pm4py.view_petri_net(net, im, fm)

base_log = simulate_model(net, im, fm, 1, BASE_LENGTH)

net, im, fm = pm4py.read_pnml("../../models/net_added_e.pnml")
#pm4py.view_petri_net(net, im, fm)

added_e_log = simulate_model(net, im, fm, 1, MODIFIED_LENGTH)

net, im, fm = pm4py.read_pnml("../../models/net_removed_b.pnml")
#pm4py.view_petri_net(net, im, fm)

removed_b_log = simulate_model(net, im, fm, 1, MODIFIED_LENGTH)

net, im, fm = pm4py.read_pnml("../../models/net_equal_prob.pnml")
#pm4py.view_petri_net(net, im, fm)

equal_prob_log = simulate_model(net, im, fm, 1, MODIFIED_LENGTH)

profile_log(base_log[0], "base_model", [0, BASE_LENGTH/STEP_SIZE],[0, BASE_LENGTH])

profile_log(added_e_log[0], "added_activity", [0, MODIFIED_LENGTH/STEP_SIZE],[0, MODIFIED_LENGTH])

profile_log(removed_b_log[0], "removed_activity", [0, MODIFIED_LENGTH/STEP_SIZE],[0, MODIFIED_LENGTH])

profile_log(equal_prob_log[0], "equal_probability", [0, MODIFIED_LENGTH/STEP_SIZE],[0, MODIFIED_LENGTH])



a = copy.deepcopy(base_log[0])
b = copy.deepcopy(added_e_log[0])
for x in b:
    a.append(x)
print(type(b))
print(len(a))
profile_log(a, "base_to_added",[0,100,200], [0,1000,2000])#

a = copy.deepcopy(base_log[0])
b = copy.deepcopy(removed_b_log[0])
for x in b:
    a.append(x)
print(type(b))
print(len(a))
profile_log(a, "base_to_removed",[0,100,200], [0,1000,2000])#

a = copy.deepcopy(base_log[0])
b = copy.deepcopy(equal_prob_log[0])
for x in b:
    a.append(x)
print(type(b))
print(len(a))
profile_log(a, "base_to_equal",[0,100,200], [0,1000,2000])#