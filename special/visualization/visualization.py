from enum import Enum

import matplotlib.pyplot as plt
import numpy as np

from special.estimation.species_estimator import SpeciesEstimator


class plot_params(Enum):
    WIDTH = 9


def plot_rank_abundance(estimator: SpeciesEstimator, species_id: str, file_name: str, abundance: bool):
    #plt.style.use('seaborn-v0_8-ticks')

    plt.rcParams['figure.figsize'] = [9, 3.5]
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20

    reference_sample = estimator.metrics[species_id].reference_sample_abundance if abundance \
        else estimator.metrics[species_id].reference_sample_incidence
    reference_values_sorted = sorted(list(reference_sample.values()), reverse=True)
    no_species = len(reference_sample)

    plt.fill_between(np.linspace(0, no_species, no_species, endpoint=False),
                     reference_values_sorted,
                     [0 for _ in range(no_species)], alpha=0.4)
    plt.plot(reference_values_sorted)
    plt.xlabel("Species Rank", fontsize=24)
    plt.ylabel("Occurrences", fontsize=24)
    plt.xticks([0, no_species - 1], [1, no_species])
    plt.yticks([0, max(reference_values_sorted)],
               [0, max(reference_values_sorted)])

    plt.tight_layout()
    plt.savefig(file_name, format="pdf")
    plt.close()


def plot_diversity_sample_vs_estimate(estimator: SpeciesEstimator, species_id: str, metrics: [str], file_name: str,
                                      abundance: bool):
    '''
    Plots the obtained time series of sample-based diversity vs asymptotic diversity
    :param estimator:
    :param species_id:
    :param metrics:
    :param file_name:
    :param abundance:
    :return:
    '''
    plt.rcParams['figure.figsize'] = [3 * 9, 5]
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    f, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey="row", layout="constrained")

    for metric, ax, title in zip(metrics, (ax1, ax2, ax3), ("D0", "D1", "D2")):
        key_sample = "abundance_sample_" + metric if abundance else "incidence_sample_" + metric
        key_estimate = "abundance_estimate_" + metric if abundance else "incidence_estimate_" + metric

        series_sample = estimator.metrics[species_id][key_sample] if abundance else estimator.metrics[species_id][
            key_sample]
        series_estimate = estimator.metrics[species_id][key_estimate] if abundance else estimator.metrics[species_id][
            key_estimate]

        series_observations_ids = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
            estimator.metrics[species_id][
                "incidence_no_observations"]

        no_data_points = len(series_sample)
        series_ticks = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
            estimator.metrics[species_id]["incidence_no_observations"]

        ax.set_xticks([0, no_data_points], [0, series_observations_ids[-1]])
        ax.plot(series_sample, label="observed")
        ax.plot(series_estimate, label="estimated")
        ax.set_title(title, fontsize=24)
        ax.set_xlabel("Sample Size", fontsize=24)
        ax.set_ylabel("Diversity", fontsize=24)
        ax.legend(fontsize=20)
    #plt.tight_layout()
    plt.savefig("diversity_series_" + file_name, format="pdf")
    return f


def plot_diversity_series(estimator: SpeciesEstimator, species_id: str, metric: str, file_name: str,
                          abundance: bool):
    '''
    Plots the time series of the specified sample_based diversity metric, adding the asymptotic diversity as an indicator
    :param estimator:
    :param species_id:
    :param metric:
    :param file_name:
    :param abundance:
    :return:
    '''
    plt.rcParams['figure.figsize'] = [5, 4]
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    f = plt.figure(layout="constrained")

    #for metric, ax, title in zip(metrics, (ax1, ax2, ax3), ("D0","D1","d2")):
    key_sample = "abundance_sample_" + metric if abundance else "incidence_sample_" + metric
    key_estimate = "abundance_estimate_" + metric if abundance else "incidence_estimate_" + metric

    series_sample = estimator.metrics[species_id][key_sample] if abundance else estimator.metrics[species_id][
        key_sample]
    series_estimate = estimator.metrics[species_id][key_estimate] if abundance else estimator.metrics[species_id][
        key_estimate]

    series_observations_ids = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
        estimator.metrics[species_id][
            "incidence_no_observations"]

    no_data_points = len(series_sample)
    series_ticks = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
        estimator.metrics[species_id]["incidence_no_observations"]

    plt.xticks([0, no_data_points], [0, series_observations_ids[-1]])

    plt.plot(series_sample, label=metric)
    plt.title(metric, fontsize=24)
    plt.xlabel("Sample Size", fontsize=24)
    plt.ylabel("Diversity", fontsize=24)
    plt.plot(no_data_points, series_estimate[-1], 'o', c="grey")
    #plt.annotate(metric + "=" + str(series_estimate[-1]), (no_data_points, series_estimate[-1]))
    #plt.tight_layout()
    plt.savefig("diversity_series_" + file_name, format="pdf")
    return f


def plot_diversity_series_all(estimator: SpeciesEstimator, species_id: str, metrics: [str], file_name: str,
                              abundance: bool):
    '''
    Plots all sample-based diversity metrics with their asymptotiv diversity
    :param estimator:
    :param species_id:
    :param metrics:
    :param file_name:
    :param abundance:
    :return:
    '''
    plt.rcParams['figure.figsize'] = [5, 4]
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    f = plt.figure(layout="constrained")

    for metric, title in zip(metrics, ("D0", "D1", "D2")):
        key_sample = "abundance_sample_" + metric if abundance else "incidence_sample_" + metric
        key_estimate = "abundance_estimate_" + metric if abundance else "incidence_estimate_" + metric

        series_sample = estimator.metrics[species_id][key_sample] if abundance else estimator.metrics[species_id][
            key_sample]
        series_estimate = estimator.metrics[species_id][key_estimate] if abundance else estimator.metrics[species_id][
            key_estimate]

        series_observations_ids = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
            estimator.metrics[species_id][
                "incidence_no_observations"]

        no_data_points = len(series_sample)
        series_ticks = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
            estimator.metrics[species_id]["incidence_no_observations"]

        plt.xticks([0, no_data_points], [0, series_observations_ids[-1]])

        plt.plot(series_sample, label=metric)
        plt.plot(no_data_points, series_estimate[-1], 'o', c="grey")
        #plt.annotate(metric + "=" + str(series_estimate[-1]), (no_data_points, series_estimate[-1]))
    plt.xlabel("Sample Size", fontsize=24)
    plt.ylabel("Diversity", fontsize=24)

    plt.tight_layout()
    plt.savefig("diversity_series_" + file_name, format="pdf")
    return f


def plot_diversity_profile(estimator: SpeciesEstimator, species_id: str, file_name: str,
                           abundance: bool):
    '''
    Plots the asymptotic diversity profile
    :param estimator:
    :param species_id:
    :param file_name:
    :param abundance:
    :return:
    '''
    plt.rcParams['figure.figsize'] = [5, 4]
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    f = plt.figure(layout="constrained")

    key = "abundance_estimate_" if abundance else "incidence_estimate_"

    profile = [estimator.metrics[species_id][key + "d0"][-1],
               estimator.metrics[species_id][key + "d1"][-1],
               estimator.metrics[species_id][key + "d2"][-1]]

    plt.xticks([0, 1, 2], ["q=0", "q=1", "q=2"])
    plt.plot(profile)
    plt.xlabel("Order q", fontsize=24)
    plt.ylabel("Diversity", fontsize=24)

    plt.tight_layout()
    plt.savefig("diversity_profile_" + file_name, format="pdf")
    return f


def plot_completeness_profile(estimator: SpeciesEstimator, species_id: str, file_name: str,
                              abundance: bool):
    '''
    Plots the time series for both completeness and coverage
    :param estimator:
    :param species_id:
    :param file_name:
    :param abundance:
    :return:
    '''
    plt.rcParams['figure.figsize'] = [5, 4]
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    f = plt.figure(layout="constrained")

    key_completeness = "abundance_c0" if abundance else "incidence_c0"
    key_coverage = "abundance_c1" if abundance else "incidence_c1"

    series_completeness = estimator.metrics[species_id][key_completeness] if abundance else \
    estimator.metrics[species_id][
        key_completeness]
    series_coverage = estimator.metrics[species_id][key_coverage] if abundance else estimator.metrics[species_id][
        key_coverage]

    series_observations_ids = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
        estimator.metrics[species_id][
            "incidence_no_observations"]

    no_data_points = len(series_completeness)
    series_ticks = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
        estimator.metrics[species_id]["incidence_no_observations"]

    plt.xticks([0, no_data_points], [0, series_observations_ids[-1]])

    plt.plot(series_completeness, label="Completeness")
    plt.plot(series_coverage, label="Coverage")
    plt.legend()
    plt.title("Completeness Profile", fontsize=24)
    plt.xlabel("Sample Size", fontsize=24)
    plt.savefig("completeness_profile_" + file_name, format="pdf")
    return f


def plot_expected_sampling_effort(estimator: SpeciesEstimator, species_id: str, file_name: str,
                                  abundance: bool):
    '''
    Plots the time series for expected sampling efforts
    :param estimator:
    :param species_id:
    :param file_name:
    :param abundance:
    :return:
    '''
    plt.rcParams['figure.figsize'] = [5, 4]
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    f = plt.figure(layout="constrained")

    key = "abundance_l_" if abundance else "incidence_l_"
    for n in estimator.l_n:
        series_effort = estimator.metrics[species_id][key + str(n)] if abundance else estimator.metrics[species_id][
            key + str(n)]

        series_observations_ids = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
            estimator.metrics[species_id][
                "incidence_no_observations"]

        no_data_points = len(series_effort)
        series_ticks = estimator.metrics[species_id]["abundance_no_observations"] if abundance else \
            estimator.metrics[species_id]["incidence_no_observations"]

        plt.xticks([0, no_data_points], [0, series_observations_ids[-1]])

        plt.plot(series_effort, label="l=" + str(n))
    plt.legend()
    #plt.title("Completeness Profile", fontsize=24)
    plt.xlabel("Sample Size", fontsize=24)
    plt.savefig("effort_" + file_name, format="pdf")
    return f
