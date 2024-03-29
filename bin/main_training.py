import argparse
import os
from train_domain2go import train_domain2go
from em_algorithm import run_em_algorithm_theta_calc, run_em_algorithm_e_calc
from enrichment_analysis import run_enrichment_analysis
from eval_preds import *

parser = argparse.ArgumentParser(description='Domain2GO arguments')

parser.add_argument(
    '--em',
    type=str,
    default="skip",
    metavar='EM_algorithm_mode',
    help='')

parser.add_argument(
    '--enrichment',
    type=str,
    default="skip",
    metavar='enrichment_analysis_mode',
    help='')

parser.add_argument(
    '--cafa_eval',
    type=str,
    default="skip",
    metavar='cafa_evaluation',
    help='')


if __name__ == "__main__":
    args = parser.parse_args()

    initial_mapping, initial_mapping_w_scores, filtered_original_mapping = train_domain2go()

    # running EM algorithm
    
    em_subset_w_scores = initial_mapping_w_scores[initial_mapping_w_scores["s"] >= 0.1]
    em_subset = initial_mapping.merge(em_subset_w_scores, on=["GO", "Interpro"])[["GO", "Interpro", "Uniprot"]]

    if args.em == "only_theta":
        run_em_algorithm_theta_calc(em_subset, em_subset_w_scores, only_theta="yes")

    elif args.em == "full_mode":
        run_em_algorithm_e_calc(em_subset, em_subset_w_scores)

    else:
        pass

    # enrichment analysis
    # if user doesn't want to perform EM algorithm, 
    # "em_full_output.txt" file can be saved into output folder to perform enrichment analysis

    if args.enrichment == "skip":
        pass
    
    elif args.enrichment != "skip":
        if args.em == "skip":
            run_enrichment_analysis(initial_mapping_w_scores, full_mode = "yes")

        if args.em == "only_theta":
            run_enrichment_analysis(initial_mapping_w_scores, full_mode = "no")

        if args.em == "full_mode":
            run_enrichment_analysis(initial_mapping_w_scores, full_mode = "yes")

    if args.cafa_eval == "skip":
        pass

    else:
        print('-------------------Performing CAFA3 evaluation-------------------')
        # create_raw_predictions(initial_mapping_w_scores)
        preds_dict = load_single_model_predictions()
        benchmark = load_benchmark()
        overall_eval_all(benchmark, preds_dict)
        organism_specific_eval_all(benchmark, preds_dict)