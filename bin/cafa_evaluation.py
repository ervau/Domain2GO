
import pandas as pd
import numpy as np
# from sklearn import metrics
import json
import os
import dask.dataframe as dd
from math import sqrt

cwd = os.getcwd()
input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
output_file_path = "{}Domain2GO/output".format(cwd.split("Domain2GO")[0])


cafa_files_path = "{}/cafa".format(input_file_path)
raw_predictions_path = "{}/cafa_raw_predictions".format(output_file_path)

ontologies = ["MFO", "CCO", "BPO"]
models = ["domain2go_s", "domain2go_e", "blast", "naive", "interpro2go_2016"]
cafa_to_domains = pd.read_csv("{}/CAFA3_targets_domains_mapped.txt".format(cafa_files_path))


def raw_predictions(model_df, model_name, score_col_name, assign_score = None, normalize = None, score_threshold = None):

    if score_threshold:
        model_scores = model_df[model_df[score_col_name] >= score_threshold] 
        predictions = pd.merge(cafa_to_domains, model_scores, on='Interpro')
    else:
        predictions = pd.merge(cafa_to_domains, model_df, on='Interpro')

    if assign_score:
        predictions["score"] = 1

    predictions= predictions.sort_values(by = score_col_name, ascending=False)
    predictions = predictions.drop_duplicates(subset = ["CAFA", "GO"], keep = "first")
    predictions = predictions[["CAFA", "GO", score_col_name]]

    if normalize:
        score_col = predictions[score_col_name]
        predictions[score_col_name] = (score_col-score_col.min())/(score_col.max()-score_col.min())
    
    predictions.to_csv("{}/{}.txt".format(raw_predictions_path, model_name ), sep= " ", index=False)
    print("Raw predictions are saved: {}/{}.txt".format(raw_predictions_path, model_name ))


def create_raw_predictions(initial_mapping_w_scores):

    print("Creating raw CAFA predictions. This will give raw predictions for all species in a single file for each model (Domain2GO-S, Domain2GO-E and InterPro2GO-2016). Species-specific parsed predictions are in input_data/cafa/predictions folder. Parsed and propagated prediction files are created by using CAFA Evaluation repo (more info in README)")

    raw_predictions(initial_mapping_w_scores, "domain2go_s", "s", score_threshold = 0.05)

    domain2go_e_model_df = pd.read_csv("{}/em_algorithm/em_full_output.txt".format(output_file_path), sep = "\t")
    raw_predictions(domain2go_e_model_df, "domain2go_e", "E", normalize = "yes")

    interpro2go_2016_model_df = pd.read_csv("{}/interpro2go_2016.txt".format(input_file_path), sep= " ", names=["Interpro", "GO"], header=None)
    raw_predictions(interpro2go_2016_model_df, "interpro2go_2016", "score", assign_score = "yes")


def read_parsed_predictions(model):
    prediction_file_path = "{}/predictions/{}_predictions/json".format(cafa_files_path, model)
    organisms_predictions = {}
    all_organisms_predictions = {}

    for ontology in ontologies:
        organisms_predictions[ontology] = []
        for file in os.listdir(prediction_file_path):
            if ontology in file:
                with open("{}/{}".format(prediction_file_path, file)) as json_file:
                    temp_dict = json.load(json_file)
                organisms_predictions[ontology].append(temp_dict)

        # all_organisms_predictions[ontology] = organisms_predictions[ontology][0].copy()
        all_organisms_predictions[ontology] = {}

        for organism in organisms_predictions[ontology]:
            all_organisms_predictions[ontology].update(organism)

    return all_organisms_predictions   


def import_cafa_benchmark():
    benchmark_file_path = "{}/benchmarks".format(cafa_files_path)
    benchmarks={}
    all_benchmark={}

    for ontology in ontologies:
        benchmarks[ontology]=[]
        for file in os.listdir(benchmark_file_path):
            if ontology in file and "karya" not in file:
                with open("{}/{}".format(benchmark_file_path, file)) as json_file:
                    temp_dict = json.load(json_file)
                benchmarks[ontology].append(temp_dict)
        
        all_benchmark[ontology]=benchmarks[ontology][0]["protein_annotations"].copy()
        for organism in benchmarks[ontology]:
            all_benchmark[ontology].update(organism["protein_annotations"])
        all_benchmark[ontology]["benchmark_ontology_term_count"]=benchmarks[ontology][0]["benchmark_ontology_term_count"]
    
    return all_benchmark

def precision(tp,fp):
    return (tp / (tp+fp))


def f_measure(pr,rc):
    return (2 * (pr*rc) / (pr+rc))


def mcc_calc(tp,tn,fp,fn):
    return (tp*tn - fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))


def accuracy(tp,tn,fp,fn):
    return (tp + tn) / (tp + tn + fp + fn)


def recall(tp,fn):
    return (tp / (tp+fn))


def score_calc(all_benchmark, ontology, thresholds, model):
    columns=['threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy']
    output_df=pd.DataFrame(columns = columns)


    for j in thresholds:
        tp=0
        fp=0
        fn=0
        tn=0
        for protein in all_benchmark[ontology].keys():
            if protein== "benchmark_ontology_term_count":
                continue
            else:
                predicted_terms = model[ontology].get(protein, {})
                predicted_annotations = {k for k, v in predicted_terms.items() if v >= j
                }

                benchmark_protein_annotation = set(all_benchmark[ontology].get(protein))
                if len(benchmark_protein_annotation) == 0:
                    continue
            
                #confusion matrix
                true_positive_terms = len(predicted_annotations & benchmark_protein_annotation)
                false_positive_terms = len(predicted_annotations - benchmark_protein_annotation)
                false_negative_terms = len(benchmark_protein_annotation - predicted_annotations)
                
                tp=tp+true_positive_terms
                fp=fp+false_positive_terms
                fn=fn+false_negative_terms
                tn=tn + (all_benchmark[ontology]["benchmark_ontology_term_count"] - (true_positive_terms + false_positive_terms + false_negative_terms))
        
        try:
            pr = precision(tp,fp)
        except ZeroDivisionError:
            pr = 0

        rc= recall(tp,fn)
        try:
            f_score= f_measure(pr,rc)
        except ZeroDivisionError:
            f_score=0
        
        try:
            mcc= mcc_calc(tp,tn,fp,fn)
        except ZeroDivisionError:
            mcc=0
        
        try:
            acc=accuracy(tp,tn,fp,fn)
        except ZeroDivisionError:
            acc=0

        output_df=output_df.append(pd.DataFrame([[j, ontology,tn, fp, fn, tp, pr, rc, f_score, mcc, acc]], columns= columns))
    
    return output_df


def partial_score_calc(all_benchmark, ontology, thresholds, model):
    columns=['threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy', 'coverage']
    output_df=pd.DataFrame(columns = columns)

    all_predicted_proteins= list(model["MFO"].keys()) + list(model["BPO"].keys()) + list(model["CCO"].keys())

    for j in thresholds:
        tp=0
        fp=0
        fn=0
        tn=0
        protein_count=0
        for protein in all_benchmark[ontology].keys():
            if protein== "benchmark_ontology_term_count" or protein not in all_predicted_proteins:
                continue
            else:
                protein_count=protein_count+1
                predicted_terms = model[ontology].get(protein, {})
                predicted_annotations = {k for k, v in predicted_terms.items() if v >= j
                }

                benchmark_protein_annotation = set(all_benchmark[ontology].get(protein))
                if len(benchmark_protein_annotation) == 0:
                    continue
            
                #confusion matrix
                true_positive_terms = len(predicted_annotations & benchmark_protein_annotation)
                false_positive_terms = len(predicted_annotations - benchmark_protein_annotation)
                false_negative_terms = len(benchmark_protein_annotation - predicted_annotations)
                
                tp=tp+true_positive_terms
                fp=fp+false_positive_terms
                fn=fn+false_negative_terms
                tn=tn + (all_benchmark[ontology]["benchmark_ontology_term_count"] - (true_positive_terms + false_positive_terms + false_negative_terms))
        
        try:
            pr = precision(tp,fp)
        except ZeroDivisionError:
            pr = 0

        rc= recall(tp,fn)
        try:
            f_score= f_measure(pr,rc)
        except ZeroDivisionError:
            f_score=0
        
        try:
            mcc= mcc_calc(tp,tn,fp,fn)
        except ZeroDivisionError:
            mcc=0
        
        try:
            acc=accuracy(tp,tn,fp,fn)
        except ZeroDivisionError:
            acc=0
        
        coverage=protein_count/(len(all_benchmark[ontology])-1)
        output_df=output_df.append(pd.DataFrame([[j, ontology,tn, fp, fn, tp, pr, rc, f_score, mcc, acc, coverage]], columns= columns))
        
    
    return output_df

def blast_domain2go_merge(all_benchmark, ontology, thresholds, blast_all_predictions, domain2go_all_predictions):
    columns=['threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy']
    output_df=pd.DataFrame(columns = columns)

    for j in thresholds:
        tp=0
        fp=0
        fn=0
        tn=0
        for protein in all_benchmark[ontology].keys():
            if protein== "benchmark_ontology_term_count":
                continue
            else:
                blast_predicted_terms = blast_all_predictions[ontology].get(protein, {})
                domain2go_predicted_terms = domain2go_all_predictions[ontology].get(protein, {})
                mutual_terms= blast_predicted_terms.keys() & domain2go_predicted_terms.keys()
                blast_only_terms= blast_predicted_terms.keys() - mutual_terms
                domain2go_only_terms = domain2go_predicted_terms.keys() - mutual_terms

                predicted_annotations= set()
                for b in blast_only_terms:
                    if blast_predicted_terms[b] >= j:
                        predicted_annotations.add(b)

                for d in domain2go_only_terms:
                    if domain2go_predicted_terms[d] >= j:
                        predicted_annotations.add(d)
                
                for m in mutual_terms:
                    if blast_predicted_terms[m] >= j or domain2go_predicted_terms[m]  >= j:
                        predicted_annotations.add(m)
            
                benchmark_protein_annotation = set(all_benchmark[ontology].get(protein))
                if len(benchmark_protein_annotation) == 0:
                    continue
        
                #confusion matrix
                true_positive_terms = len(predicted_annotations & benchmark_protein_annotation)
                false_positive_terms = len(predicted_annotations - benchmark_protein_annotation)
                false_negative_terms = len(benchmark_protein_annotation - predicted_annotations)

                tp=tp+true_positive_terms
                fp=fp+false_positive_terms
                fn=fn+false_negative_terms
                tn=tn + (all_benchmark[ontology]["benchmark_ontology_term_count"] - (true_positive_terms + false_positive_terms + false_negative_terms))
                
        try:
            pr = precision(tp,fp)
        except ZeroDivisionError:
            pr = 0

        rc= recall(tp,fn)

        try:
            f_score= f_measure(pr,rc)
        except ZeroDivisionError:
            f_score=0
        
        try:
            mcc= mcc_calc(tp,tn,fp,fn)
        except ZeroDivisionError:
            mcc=0
        
        try:
            acc=accuracy(tp,tn,fp,fn)
        except ZeroDivisionError:
            acc=0
        
        output_df=output_df.append(pd.DataFrame([[j, ontology,tn, fp, fn, tp, pr, rc, f_score, mcc, acc]], columns= columns))

    return output_df

def blast_domain2go_merge_partial(all_benchmark, ontology, thresholds, blast_all_predictions, domain2go_all_predictions):
    columns=['threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy', 'coverage']
    output_df=pd.DataFrame(columns = columns)

    all_predicted_proteins= list(blast_all_predictions["MFO"].keys()) + list(blast_all_predictions["BPO"].keys()) + list(blast_all_predictions["CCO"].keys()) +list(domain2go_all_predictions["MFO"].keys()) +list(domain2go_all_predictions["BPO"].keys()) + list(domain2go_all_predictions["CCO"].keys())

    for j in thresholds:
        tp=0
        fp=0
        fn=0
        tn=0
        protein_count=0
        for protein in all_benchmark[ontology].keys():
            if protein== "benchmark_ontology_term_count" or protein not in all_predicted_proteins:
                continue
            else:
                protein_count=protein_count+1
                blast_predicted_terms = blast_all_predictions[ontology].get(protein, {})
                domain2go_predicted_terms = domain2go_all_predictions[ontology].get(protein, {})
                mutual_terms= blast_predicted_terms.keys() & domain2go_predicted_terms.keys()
                blast_only_terms= blast_predicted_terms.keys() - mutual_terms
                domain2go_only_terms = domain2go_predicted_terms.keys() - mutual_terms

                predicted_annotations= set()
                for b in blast_only_terms:
                    if blast_predicted_terms[b] >= j:
                        predicted_annotations.add(b)

                for d in domain2go_only_terms:
                    if domain2go_predicted_terms[d] >= j:
                        predicted_annotations.add(d)
                
                for m in mutual_terms:
                    if blast_predicted_terms[m] >= j or domain2go_predicted_terms[m]  >= j:
                        predicted_annotations.add(m)
            
                benchmark_protein_annotation = set(all_benchmark[ontology].get(protein))
                if len(benchmark_protein_annotation) == 0:
                    continue
        
                #confusion matrix
                true_positive_terms = len(predicted_annotations & benchmark_protein_annotation)
                false_positive_terms = len(predicted_annotations - benchmark_protein_annotation)
                false_negative_terms = len(benchmark_protein_annotation - predicted_annotations)

                tp=tp+true_positive_terms
                fp=fp+false_positive_terms
                fn=fn+false_negative_terms
                tn=tn + (all_benchmark[ontology]["benchmark_ontology_term_count"] - (true_positive_terms + false_positive_terms + false_negative_terms))
                
        try:
            pr = precision(tp,fp)
        except ZeroDivisionError:
            pr = 0

        rc= recall(tp,fn)

        try:
            f_score= f_measure(pr,rc)
        except ZeroDivisionError:
            f_score=0
        
        try:
            mcc= mcc_calc(tp,tn,fp,fn)
        except ZeroDivisionError:
            mcc=0
        
        try:
            acc=accuracy(tp,tn,fp,fn)
        except ZeroDivisionError:
            acc=0

        coverage=protein_count/(len(all_benchmark[ontology])-1)
        output_df=output_df.append(pd.DataFrame([[j, ontology,tn, fp, fn, tp, pr, rc, f_score, mcc, acc,coverage]], columns= columns))

    return output_df

def fmax_calc(thresholds = np.arange(0.01,1,0.01)):
    columns=['model', 'ontology', 'evaluation_mode', 'fmax']
    output_df=pd.DataFrame(columns = columns)

    all_benchmark = import_cafa_benchmark()
    
    for model in models:
        model_all_predictions = read_parsed_predictions(model)
        for ontology in ontologies:
            fmax_score = max(score_calc(all_benchmark, ontology, thresholds, model_all_predictions)["f1"])
            output_df=output_df.append(pd.DataFrame([[model, ontology, 'full_mode', fmax_score]], columns= columns))
        
        for ontology in ontologies:
            fmax_score = max(partial_score_calc(all_benchmark, ontology, thresholds, model_all_predictions)["f1"])
            output_df=output_df.append(pd.DataFrame([[model, ontology, 'partial_mode', fmax_score]], columns= columns))

    blast_all_predictions = read_parsed_predictions("blast")
    domain2go_all_predictions = read_parsed_predictions("domain2go_s")
    for ontology in ontologies:
        fmax_score = max(blast_domain2go_merge(all_benchmark, ontology, thresholds, blast_all_predictions, domain2go_all_predictions)["f1"])
        output_df=output_df.append(pd.DataFrame([["domain2go-blast", ontology, 'full_mode', fmax_score]], columns= columns))

    for ontology in ontologies:
        fmax_score = max(blast_domain2go_merge_partial(all_benchmark, ontology, thresholds, blast_all_predictions, domain2go_all_predictions)["f1"])
        output_df=output_df.append(pd.DataFrame([["domain2go-blast", ontology, 'partial_mode', fmax_score]], columns= columns))
    
    output_file = "{}/cafa_fmax_results.txt".format(output_file_path)
    print("CAFA Fmax score results are saved: {}".format(output_file))
    output_df.to_csv(output_file)
