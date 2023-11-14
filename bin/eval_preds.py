
import pandas as pd
import numpy as np
import json
import os
from math import sqrt
from cafa_evaluation import *

cwd = os.getcwd()
input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
output_file_path = "{}Domain2GO/output".format(cwd.split("Domain2GO")[0])

cafa_files_path = "{}/cafa".format(input_file_path)
raw_predictions_path = "{}/cafa_raw_predictions".format(output_file_path)
cafa_to_domains = pd.read_csv("{}/CAFA3_targets_domains_mapped.txt".format(cafa_files_path))


def raw_predictions(model_df, model_name, score_col_name, assign_score = None, normalize = None, score_threshold = None):

    if score_threshold:
        model_scores = model_df[model_df[score_col_name] >= score_threshold] 
        predictions = pd.merge(cafa_to_domains, model_scores, on='Interpro')
    else:
        predictions = pd.merge(cafa_to_domains, model_df, on='Interpro')

    if assign_score:
        predictions["score"] = 1

    go_column_name = "GO" if model_name in ['domain2go_s', 'interpro2go_2016'] else "GO_ID"
    predictions= predictions.sort_values(by = score_col_name, ascending=False)
    predictions = predictions.drop_duplicates(subset = ["CAFA", go_column_name], keep = "first")
    predictions = predictions[["CAFA", go_column_name, score_col_name]]

    if normalize:
        score_col = predictions[score_col_name]
        predictions[score_col_name] = (score_col-score_col.min())/(score_col.max()-score_col.min())
    
    predictions.to_csv("{}/{}.txt".format(raw_predictions_path, model_name ), sep= " ", index=False)
    print("{} raw predictions are saved: {}/{}.txt".format(model_name[0].upper() + model_name[1:], raw_predictions_path, model_name ))

def create_raw_predictions(initial_mapping_w_scores):

    print("Creating raw CAFA predictions. This will give raw predictions for all species in a single file for each model (Domain2GO-S, Domain2GO-E and InterPro2GO-2016). Species-specific parsed predictions are in input_data/cafa/predictions folder. Parsed and propagated prediction files are created by using CAFA Evaluation repo (more info in README)")

    raw_predictions(initial_mapping_w_scores, "domain2go_s", "s", score_threshold = 0.05)

    domain2go_e_model_df = pd.read_csv("{}/em_algorithm/em_full_output.txt".format(output_file_path), sep = "\t")
    raw_predictions(domain2go_e_model_df, "domain2go_e", "E", normalize = "yes")

    interpro2go_2016_model_df = pd.read_csv("{}/interpro2go_2016.txt".format(input_file_path), sep= " ", names=["Interpro", "GO"], header=None)
    raw_predictions(interpro2go_2016_model_df, "interpro2go_2016", "score", assign_score = "yes")
    
def load_single_model_predictions(models=['domain2go_s', 'domain2go_e', 'interpro2go_2016', 'blast']):
    print('Loading predictions for all models...')
    pred_dict = {}
    for model in models:
        json_pred_file_path = os.path.join(cafa_files_path, 'predictions', f'{model}_predictions', 'json')
        pred_dict[model] = import_overall_preds(json_pred_file_path)
    
    return pred_dict

def load_benchmark():
    print('Loading CAFA3 benchmark...')
    parsed_bencmark_path = os.path.join(cafa_files_path, 'benchmarks')
    benchmark = import_benchmark(parsed_bencmark_path)
    
    return benchmark


def overall_eval_all(benchmark,  pred_dict, ontologies=["MFO", "CCO", "BPO"], thresholds=np.arange(0.01,1,0.01)):
    columns = ['ontology', 'evaluation_mode', 'model', 'fmax']
    output_df = pd.DataFrame(columns=columns)
    
    print("Performing overall evaluation for all models...")
    
    for ontology in ontologies:
        
        # single models
        for k, v in pred_dict.items():
            if k != 'blast':
                    model_ontology_fmax_full = max(overall_eval(ontology, thresholds, benchmark, v)["f1"])
                    model_ontology_fmax_partial = max(overall_eval(ontology, thresholds, benchmark, v, partial=True)["f1"])
                    output_df = output_df.append(pd.DataFrame([[ontology, "full", k, model_ontology_fmax_full]], columns=columns))
                    output_df = output_df.append(pd.DataFrame([[ontology, "partial", k, model_ontology_fmax_partial]], columns=columns))
    
        # merged models 
        s_blast_full_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['blast']])['f1'])
        s_blast_partial_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['blast']], partial=True)['f1'])
        output_df = output_df.append(pd.DataFrame([[ontology, "full", "s_blast", s_blast_full_fmax]], columns=columns))
        output_df = output_df.append(pd.DataFrame([[ontology, "partial", "s_blast", s_blast_partial_fmax]], columns=columns))
        
        s_e_blast_full_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['blast']])['f1'])
        s_e_blast_partial_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['blast']], partial=True)['f1'])
        output_df = output_df.append(pd.DataFrame([[ontology, "full", "s_e_blast", s_e_blast_full_fmax]], columns=columns))
        output_df = output_df.append(pd.DataFrame([[ontology, "partial", "s_e_blast", s_e_blast_partial_fmax]], columns=columns))
        
        s_e_ipr2go_blast_full_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['interpro2go_2016'], pred_dict['blast']])['f1'])
        s_e_ipr2go_blast_partial_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['interpro2go_2016'], pred_dict['blast']], partial=True)['f1'])
        output_df = output_df.append(pd.DataFrame([[ontology, "full", "s_e_ipr2go_blast", s_e_ipr2go_blast_full_fmax]], columns=columns))
        output_df = output_df.append(pd.DataFrame([[ontology, "partial", "s_e_ipr2go_blast", s_e_ipr2go_blast_partial_fmax]], columns=columns))
        
        s_e_ipr2go_full_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['interpro2go_2016']])['f1'])
        s_e_ipr2go_partial_fmax = max(eval_merged_model(ontology, thresholds, benchmark, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['interpro2go_2016']], partial=True)['f1'])
        output_df = output_df.append(pd.DataFrame([[ontology, "full", "s_e_ipr2go", s_e_ipr2go_full_fmax]], columns=columns))
        output_df = output_df.append(pd.DataFrame([[ontology, "partial", "s_e_ipr2go", s_e_ipr2go_partial_fmax]], columns=columns))
        
    output_df.to_csv("{}/overall_fmax_results.txt".format(output_file_path), index=False)
    print("Overall evaluation results are saved: {}/overall_fmax_results.txt".format(output_file_path))


def organism_specific_eval_all(benchmark, pred_dict, ontologies=["MFO", "CCO", "BPO"], thresholds=np.arange(0.01,1,0.01)):
    columns = ['model','organism', 'ontology', 'fmax']
    output_df = pd.DataFrame(columns=columns)
    print("Performing organism specific evaluation for all models...")

    for ontology in ontologies:
        organism_list = benchmark[ontology].keys()
        for organism in organism_list:
            if organism != 'benchmark_ontology_term_count':
                
                # single models
                for k,v in pred_dict.items():
                    try:
                        fmax = max(organism_spec_eval(ontology, thresholds, benchmark, v, organism)["f1"])
                    except KeyError:
                        fmax = '-'
                    output_df = output_df.append(pd.DataFrame([[k, organism, ontology, fmax]], columns=columns))
                
                # merged models
                
                s_e_blast_fmax = max(organism_specific_merged_eval(ontology, thresholds, benchmark, organism, [pred_dict['domain2go_s'], pred_dict['blast']])["f1"])
                output_df = output_df.append(pd.DataFrame([["s_blast", organism, ontology, s_e_blast_fmax]], columns=columns))
                
                s_e_blast_fmax = max(organism_specific_merged_eval(ontology, thresholds, benchmark, organism, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['blast']])["f1"])
                output_df = output_df.append(pd.DataFrame([["s_e_blast", organism, ontology, s_e_blast_fmax]], columns=columns))
                
                s_e_ipr2go_blast_fmax = max(organism_specific_merged_eval(ontology, thresholds, benchmark, organism, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['interpro2go_2016'], pred_dict['blast']])["f1"])
                output_df = output_df.append(pd.DataFrame([["s_e_ipr2go_blast", organism, ontology, s_e_ipr2go_blast_fmax]], columns=columns))
                
                s_e_ipr2go_fmax = max(organism_specific_merged_eval(ontology, thresholds, benchmark, organism, [pred_dict['domain2go_s'], pred_dict['domain2go_e'], pred_dict['interpro2go_2016']])["f1"])
                output_df = output_df.append(pd.DataFrame([["s_e_ipr2go", organism, ontology, s_e_ipr2go_fmax]], columns=columns))
                
    output_df.to_csv("{}/organism_specific_fmax_results.txt".format(output_file_path), index=False)
    print("Organism specific evaluation results are saved: {}/organism_specific_fmax_results.txt".format(output_file_path))
