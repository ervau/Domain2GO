import pandas as pd

column_names = ["GO", "Interpro", "n", "n_go", "n_ip", "s"]

def co_occurence_similarity(together, ip, go):
    return (2*together) / (go + ip)

def original_score_cal(mapping):
    mapping_w_scores = pd.DataFrame(columns = column_names)

    # number of proteins annotated with this go term 
    go_counts = mapping[["GO", "Uniprot"]].drop_duplicates()['GO'].value_counts()

    # number of proteins annotated with this interpro term 
    ip_counts = mapping[["Interpro", "Uniprot"]].drop_duplicates()['Interpro'].value_counts()

    # N (no_proteins_where_terms_annotated_together)
    mapping_w_n = mapping[["GO", "Interpro"]].groupby(mapping[["GO", "Interpro"]].columns.tolist(),as_index=False).size()
    
    n_go_column= []
    for i in mapping_w_n["GO"]:
        n_go_column.append(go_counts[i])
    
    n_ip_column = []
    for i in mapping_w_n["Interpro"]:
        n_ip_column.append(ip_counts[i])

    similarity_scores = []
    for i in range(len(mapping_w_n)):
        similarity_scores.append(co_occurence_similarity(mapping_w_n["size"][i],  n_ip_column[i], n_go_column[i]))

    mapping_w_scores["GO"] = mapping_w_n["GO"]
    mapping_w_scores["Interpro"] = mapping_w_n["Interpro"]
    mapping_w_scores["n"] = mapping_w_n["size"]
    mapping_w_scores["n_go"] = n_go_column
    mapping_w_scores["n_ip"] = n_ip_column
    mapping_w_scores["s"] = similarity_scores

    return mapping_w_scores

def random_score_cal(original_mapping, random_mapping):
    mapping_w_scores = pd.DataFrame(columns = column_names)

    go_counts = original_mapping[["GO", "Uniprot"]].drop_duplicates()['GO'].value_counts()
    ip_counts = original_mapping[["Interpro", "Uniprot"]].drop_duplicates()['Interpro'].value_counts()
    
    # no of proteins where terms are annotated together will be change due to randomization
    random_mapping_w_n = random_mapping[["GO", "Interpro"]].groupby(random_mapping[["GO", "Interpro"]].columns.tolist(),as_index=False).size()
    
    n_go_column = []
    for i in random_mapping_w_n["GO"]:
        n_go_column.append(go_counts[i])
    
    n_ip_column = []
    for i in random_mapping_w_n["Interpro"]:
        n_ip_column.append(ip_counts[i])
    
    similarity_scores = []
    for i in range(len(random_mapping_w_n)):
        similarity_scores.append(co_occurence_similarity(random_mapping_w_n["size"][i],  n_ip_column[i], n_go_column[i]))

    mapping_w_scores["GO"] = random_mapping_w_n["GO"]
    mapping_w_scores["Interpro"] = random_mapping_w_n["Interpro"]
    mapping_w_scores["n"] = random_mapping_w_n["size"]
    mapping_w_scores["n_go"] = n_go_column
    mapping_w_scores["n_ip"] = n_ip_column
    mapping_w_scores["s"] = similarity_scores

    return mapping_w_scores

