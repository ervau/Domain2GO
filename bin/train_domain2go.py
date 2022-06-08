import argparse
import os
from data_prep import data_process, random_mapping_create, threshold_original_mapping
from ks_test import ks_test_calc
from enrichment_analysis import run_enrichment_analysis


# cwd = os.getcwd()
# input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
# output_file_path = "{}Domain2GO/output".format(cwd.split("Domain2GO")[0])



def train_domain2go():
# generating Domain2GO mappings
    initial_mapping, initial_mapping_w_scores = data_process()

    print("Generated initial mapping, starting randomization")
    # randomization
    random_mapping = random_mapping_create(initial_mapping)

    # performing KS-test
    ks_output_file = ks_test_calc(initial_mapping_w_scores, random_mapping)
    print("KS_test is done, results are saved: {}".format(ks_output_file))

    # filtering initial Domain2GO mappings based on KS-test results
    filtered_original_mapping, filtered_output_file = threshold_original_mapping(initial_mapping_w_scores)
    print("Filtered initial mappings based on KS-test results, finalized mappings are saved: {}".format(filtered_output_file))

    return initial_mapping, initial_mapping_w_scores, filtered_original_mapping




