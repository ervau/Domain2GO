import pandas as pd
import numpy as np
import os

cwd = os.getcwd()
input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
output_file_path = "{}Domain2GO/output".format(cwd.split("Domain2GO")[0])



top_ks = [100, 250, 750] +  list(range(1000,10000, 250)) + list(range(10000,76000, 1000))

def calculate_enrichment_analysis(score_dataframe, score_label, curated_data, manual_probability):
    
    ranked_mappings = score_dataframe.sort_values(by = score_label, ignore_index = True, ascending = False)
    
    enrichment_ratios = []

    for k in top_ks:
        enrichment_at_k = ranked_mappings.iloc[0:k].merge(curated_data, how = 'left', indicator = True)
        observed_at_k = len(enrichment_at_k[enrichment_at_k["_merge"] == "both"])
        observed_expected_ratio = observed_at_k/(manual_probability * k)
        enrichment_ratios.append(observed_expected_ratio)

    return enrichment_ratios

def run_enrichment_analysis(mapping_w_scores, full_mode = "no"):

    if full_mode == "no":
        em_output_file = "{}/em_algorithm/em_only_theta_output.txt".format(output_file_path)
        em_output = pd.read_csv(em_output_file, delimiter = "\t")
    
    else:
        em_output_file = "{}/em_algorithm/em_full_output.txt".format(output_file_path)
        em_output = pd.read_csv(em_output_file, delimiter = "\t")

    # column_names = ["S", "initial_theta", "theta", "E"]
    # result_df=pd.DataFrame(columns = column_names, index=top_ks)
    result_df=pd.DataFrame(index=top_ks)

    curated_data_file = "{}/interpro2go_2021.txt".format(input_file_path)
    curated_data = pd.read_csv(curated_data_file, delimiter = " ")
    curated_data.columns = ["GO", "Interpro"]

    subset_mapping = mapping_w_scores[mapping_w_scores["s"] >= 0.1]

    manual_merge=subset_mapping.merge(curated_data, how='left', indicator=True)[subset_mapping.merge(curated_data, how='left', indicator=True)["_merge"] == "both"]
    manual_probability=len(manual_merge)/len(subset_mapping)
    
    s_enrichment = calculate_enrichment_analysis(subset_mapping, ["s", "n"], curated_data, manual_probability)
    initial_theta_enrichment = calculate_enrichment_analysis(em_output, "initial_theta", curated_data, manual_probability)
    theta_enrichment = calculate_enrichment_analysis(em_output, "theta", curated_data, manual_probability)
    e_enrichment = []
    if full_mode == "yes":
        e_enrichment = calculate_enrichment_analysis(em_output, "E", curated_data, manual_probability)

    result_df["S"] = s_enrichment
    result_df["initial_theta"] = initial_theta_enrichment
    result_df["theta"] = theta_enrichment
    if full_mode == "yes":
        result_df["E"] = e_enrichment

    output_file = "{}/enrichment_analysis_results.txt".format(output_file_path)
    result_df.to_csv(output_file)
    print("Enrichment analysis is done, results are saved: {}".format(output_file))

    return s_enrichment, initial_theta_enrichment, theta_enrichment, e_enrichment

