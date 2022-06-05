import pandas as pd
import dask.dataframe as dd
import os
from score_calculation import random_score_cal, original_score_cal



cwd = os.getcwd()
input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
output_file_path = "{}Domain2GO/outputs".format(cwd.split("Domain2GO")[0])

# preprocessing GO data

def data_process():
    go_file = "{}/processed_go_data.txt".format(input_file_path)
    go_data = pd.read_csv(go_file, delimiter = " ", names=["GO", "Uniprot"])

    ipr_file = "{}/protein2ipr_onlydomains.txt".format(input_file_path)
    ipr_data = dd.read_csv(ipr_file, delimiter = " ", names=["Uniprot", "Interpro"])

    ipr_data = ipr_data[ipr_data['Uniprot'].isin(go_data['Uniprot'])]
    ipr_data=ipr_data.compute()
    new_go_data = go_data[go_data['Uniprot'].isin(ipr_data['Uniprot'])]
    merged = pd.merge(new_go_data, ipr_data, on='Uniprot')
    merged_w_scores = original_score_cal(merged)
    # merged = ipr_data.merge(go_data, on="Uniprot")
    # merged = merged.compute() 
    output_file = "{}/initial_mapping.txt".format(output_file_path)
    output_file_w_scores = "{}/initial_mapping_w_scores.txt".format(output_file_path)
    merged.to_csv(output_file, index=False)
    merged_w_scores.to_csv(output_file_w_scores, index=False)
    merged.reset_index(inplace=False)
    merged_w_scores.reset_index(inplace=False)

    return merged, merged_w_scores


def random_mapping_create(original_mapping):
    column_names = ["GO", "Interpro", "n", "n_go", "n_ip", "s"]
    random_mapping = pd.DataFrame(columns = column_names)

    ### SKOR HESAPLAMA FONKSİYONU YAZ, ONU BURADA KULLAN
    for x in range(1,11):
        random_iter = pd.DataFrame(columns = column_names)

        # randomization of go and interpro columns
        random_iter["GO"] = original_mapping['GO'].sample(frac=1, random_state= x).values
        random_iter['Interpro'] = original_mapping['Interpro'].sample(frac=1, random_state=10 + x).values

        # calculating scores
        random_iter_w_scores = random_score_cal(original_mapping, random_iter)
        random_mapping = random_mapping.append(random_iter_w_scores)

    output_file = "{}/randomized_mapping.txt".format(output_file_path)
    random_mapping.to_csv(output_file, index=False)
    random_mapping.reset_index(inplace=False)

    return random_mapping

def threshold_original_mapping(original_mapping_w_scores):
    filtered_original = original_mapping_w_scores.loc[(original_mapping_w_scores["n"] >= 2 ) & (original_mapping_w_scores["s"] > 0.2 )]
    output_file = "{}/filtered_original.txt".format(output_file_path)
    filtered_original.to_csv(output_file, index=False)

    return filtered_original, output_file
