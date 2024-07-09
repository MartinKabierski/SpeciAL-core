# ---EXTENSION---#
# TODO bootstrap rarefaction and extrapolation!
# TODO measure generalization of discovered process models
import statistics

import pandas as pd

import src.estimation.species_estimator
import src.estimation.species_retrieval
from src.estimation import species_estimator, species_retrieval
from src.simulation.simulation import simulate_model
from src.visualizations.plots import plot_rank_abundance, plot_all_stats
from src.simulation import simulation
from src.estimation.species_estimator import SpeciesEstimator
from functools import partial
from tqdm import tqdm

import pm4py

ESTIMATORS = ["1-gram", "2-gram", "3-gram", "4-gram", "5-gram", "trace_variants"]

def summarize_metrics_reference_samples(estimations):
    metrics_stats = {}
    for metric in species_estimator.METRICS:
        metrics_reference_sample = sorted([getattr(est,metric)[-1] for est in estimations])
        metrics_stats[metric]=[
            min(metrics_reference_sample),
            metrics_reference_sample[2],
            statistics.mean(metrics_reference_sample),
            statistics.stdev(metrics_reference_sample),
            statistics.median(metrics_reference_sample),
            metrics_reference_sample[-2],
            max(metrics_reference_sample)
        ]
        print(metric + ": " + str(statistics.median(metrics_reference_sample)))
    return pd.DataFrame.from_dict(metrics_stats, orient='index', columns=["min","lower_ci","mean","stddev","median","upper_ci","max"])


def profile_logs(logs):
    #TODO make this a dictionary
    estimations_per_species = [[], [], [], [], [], []]
    for log in tqdm(logs, desc='Evaluating logs'):
        results = profile_log(log)
        for estimator, results in zip(results, estimations_per_species):
            results.append(estimator)
    return estimations_per_species


def profile_log(log):
    est_1gram = SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=1), step_size=20)
    est_2gram = SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=2), step_size=20)
    est_3gram = SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=3), step_size=20)
    est_4gram = SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=4), step_size=20)
    est_5gram = SpeciesEstimator(partial(species_retrieval.retrieve_species_n_gram, n=5), step_size=20)
    est_tv = SpeciesEstimator(species_retrieval.retrieve_species_trace_variant, step_size=20)


    estimators = [est_1gram, est_2gram, est_3gram, est_4gram, est_5gram, est_tv]

    # iterate over traces in log
    for est in estimators:
        est.apply(log)
    return estimators


def evaluate_model(path, name, repetitions, log_size):
    net, im, fm = pm4py.read_pnml(path)
    print(name)
    simulated_logs = simulate_model(net, im, fm, repetitions, log_size)
    results_per_species = profile_logs(simulated_logs)

    for result, estimator_name in zip(results_per_species, ESTIMATORS):
        print(name + " - " + estimator_name)
        #TODO summary statistics post-hoc to the csv as df.describe() and df.info()
        summary_df = summarize_metrics_reference_samples(result)
        summary_df.to_csv("results/baseline_eval_"+name+"_"+estimator_name+"_summary.csv")
        plot_rank_abundance(result[0], name + "_" + estimator_name)
        print()
        print("Total Abundances: " + str(result[0].total_number_species_abundances))
        print("Total Incidences: " + str(result[0].total_number_species_incidences))
        print("Degree of Aggregation: " + str(result[0].degree_spatial_aggregation))
        print("Observations Incidence: " + str(result[0].number_observations_incidence))
        print("Species Counts: ")
        [print(x, end=" ") for x in result[0].reference_sample_incidence.values()]
        print()

        plot_all_stats(result, name + "_" + estimator_name, [0, 50, 100], [0,1000,2000])


#evaluate_model("./models/model_1.pnml", "baseline_eval_model_1", 100, 2000)
#evaluate_model("./models/model_2.pnml", "baseline_eval_model_2", 100, 2000)
#evaluate_model("./models/model_4.pnml", "baseline_eval_model_4", 100, 1000)
#evaluate_model("./models/model_5.pnml", "baseline_eval_model_5", 100, 1000)
#evaluate_model("./models/model_6.pnml", "baseline_eval_model_6", 100, 1000)

#Used in Evaluation
evaluate_model("../../models/model_3.pnml", "baseline_eval_model_3", 100, 2000)
#evaluate_model("./models/model_7.pnml", "baseline_eval_model_7_LARGE", 100, 5000)
#evaluate_model("./models/morel_8.pnml", "baseline_eval_model_8_LARGE", 100, 5000)
