import os
import json
import pandas as pd

def import_overall_preds(json_pred_file_path):

    preds = {}
    for ontology in ['MFO', 'CCO', 'BPO']:
        preds[ontology] = {}
        num_proteins = 0
        for file in os.listdir(json_pred_file_path):            
            if ontology in file:
                organism = file.split(f'_{ontology}')[0].split('_')[-1]
                organism_dict = json.load(open(os.path.join(json_pred_file_path, file)))
                preds[ontology][organism] = organism_dict
                num_proteins += len(organism_dict)
    return preds


def import_benchmark(parsed_bencmark_path):
    benchmark = {}
    for ontology in ['MFO', 'CCO', 'BPO']:
        benchmark[ontology] = {}
        num_proteins = 0
        for file in os.listdir(parsed_bencmark_path):
                if ontology in file and 'karya' not in file:
                    benchmark_ontology_term_count = json.load(open(os.path.join(parsed_bencmark_path, file)))['benchmark_ontology_term_count']
                    organism = file.split('_')[2]
                    organism_dict = json.load(open(os.path.join(parsed_bencmark_path, file)))['protein_annotations']
                    benchmark[ontology][organism] = organism_dict
                    num_proteins += len(organism_dict)
        benchmark[ontology]['benchmark_ontology_term_count'] = benchmark_ontology_term_count

    return benchmark


def precision(tp,fp):
    return (tp / (tp+fp))

def f_measure(pr,rc):
    return (2 * (pr*rc) / (pr+rc))

from math import sqrt

def mcc_calc(tp,tn,fp,fn):
    return (tp*tn - fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

def accuracy(tp,tn,fp,fn):
    return (tp + tn) / (tp + tn + fp + fn)

def recall(tp,fn):
    return (tp / (tp+fn))


def pred_merger(pred_list, ontology, overall=True):
    
    all_preds_list = []
    merged_preds = {}
    protein_set = []
    
    if overall:
        for i, pred in enumerate(pred_list):
            all_preds = {}
            for organism in pred_list[i][ontology].keys():
                all_preds = {**all_preds, **pred_list[i][ontology][organism]}
            protein_set.extend(list(all_preds.keys()))
            all_preds_list.append(all_preds)
    
    else:
        for i, pred in enumerate(pred_list):
            all_preds = {}
            all_preds = {**all_preds, **pred_list[i]}
            protein_set.extend(list(all_preds.keys()))
            all_preds_list.append(all_preds)
            
    protein_set = set(protein_set)
    
    for protein in protein_set:
        merged_preds[protein] = {}
        pred_list = []
        for pred in all_preds_list:
            pred_list.extend(list(pred.get(protein, {}).keys()))
        
        for term in pred_list:
            max_pred = 0
            for pred in all_preds_list:
                temp_pred = pred.get(protein, {}).get(term, 0)
                if temp_pred > max_pred:
                    max_pred = temp_pred
            if max_pred > 0:
                merged_preds[protein][term] = max_pred
            
    return merged_preds



def evaluator(benchmark, preds, threshold, benchmark_ontology_term_count):
    
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    predicted_proteins = 0
    
    for protein in benchmark.keys():
        if protein == 'benchmark_ontology_term_count':
            continue
        else:
            predicted_terms = preds.get(protein, {})
            if predicted_terms:
                predicted_proteins += 1
            predicted_annotations = {k for k, v in predicted_terms.items() if v >= threshold
                }
            
            benchmark_terms = set(benchmark.get(protein))
            
            if len(benchmark_terms) == 0:
                continue   
            
            tp += len(predicted_annotations & benchmark_terms)
            fp += len(predicted_annotations - benchmark_terms)
            fn += len(benchmark_terms - predicted_annotations)
            tn += benchmark_ontology_term_count - (tp + fp + fn)
            
    try:
        pr = precision(tp,fp)
    except ZeroDivisionError:
        pr = 0
        
    rc = recall(tp,fn)
    
    try:
        f1 = f_measure(pr,rc)
    except ZeroDivisionError:
        f1 = 0
        
    try:
        mcc = mcc_calc(tp,tn,fp,fn)
    except:
        mcc = 0
        
    try:
        acc = accuracy(tp,tn,fp,fn)
    except ZeroDivisionError:
        acc = 0
        
    coverage = predicted_proteins / len(benchmark.keys())
        
    return tn, fp, fn, tp, pr, rc, f1, mcc, acc, coverage

def overall_eval(ontology, thresholds, benchmark, preds, partial = False):
    columns = ['threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy']
    output_df = pd.DataFrame(columns=columns)

    # benchmark_ontology_term_count = benchmark[ontology]['benchmark_ontology_term_count']

    all_benchmark = {}
    all_benchmark[ontology] = {}

    for organism in benchmark[ontology].keys():
        if organism == 'benchmark_ontology_term_count':
            continue
        else:
            all_benchmark[ontology] = {**all_benchmark[ontology], **benchmark[ontology][organism]}
    
    all_preds = {}
    for organism in preds[ontology].keys():
        all_preds = {**all_preds, **preds[ontology][organism]}

    if partial:
        all_benchmark[ontology] = {k: v for k, v in all_benchmark[ontology].items() if k in all_preds.keys()}

    benchmark_ontology_terms = []
    for k, v in all_benchmark[ontology].items():
        benchmark_ontology_terms.extend(v)
    benchmark_ontology_term_count = len(set(benchmark_ontology_terms))
    
    for threshold in thresholds:
        
        tn, fp, fn, tp, pr, rc, f1, mcc, acc, coverage = evaluator(all_benchmark[ontology], all_preds, threshold, benchmark_ontology_term_count)

        output_df = output_df.append({'threshold': threshold, 'ontology': ontology, 'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp, 'precision': pr, 'recall': rc, 'f1': f1, 'mcc': mcc, 'accuracy': acc, 'coverage': coverage}, ignore_index=True)

    return output_df
            

def organism_spec_eval(ontology, thresholds, benchmark, preds, organism):
    columns = ['organism', 'threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy']
    output_df = pd.DataFrame(columns=columns)


    _benchmark = benchmark[ontology][organism].copy()
    benchmark_ontology_term_count = benchmark[ontology]['benchmark_ontology_term_count']
    preds = preds[ontology][organism]

    for threshold in thresholds:
        
        tn, fp, fn, tp, pr, rc, f1, mcc, acc, coverage = evaluator(_benchmark, preds, threshold, benchmark_ontology_term_count)

        output_df = output_df.append({'organism': organism, 'threshold': threshold, 'ontology': ontology, 'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp, 'precision': pr, 'recall': rc, 'f1': f1, 'mcc': mcc, 'accuracy': acc, 'coverage': coverage}, ignore_index=True)

    return output_df


def eval_merged_model(ontology, thresholds, benchmark, pred_list, partial =False):
    columns = ['threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy']
    output_df = pd.DataFrame(columns=columns)
    
    all_benchmark = {}
    all_benchmark[ontology] = {}

    for organism in benchmark[ontology].keys():
        if organism == 'benchmark_ontology_term_count':
            continue
        else:
            all_benchmark[ontology] = {**all_benchmark[ontology], **benchmark[ontology][organism]}

    merged_preds = pred_merger(pred_list, ontology)
    
    if partial:
        all_benchmark[ontology] = {k: v for k, v in all_benchmark[ontology].items() if k in merged_preds.keys()}

    benchmark_ontology_terms = []
    for k, v in all_benchmark[ontology].items():
        benchmark_ontology_terms.extend(v)
    benchmark_ontology_term_count = len(set(benchmark_ontology_terms))
    
    for threshold in thresholds:
        
        tn, fp, fn, tp, pr, rc, f1, mcc, acc, coverage = evaluator(all_benchmark[ontology], merged_preds, threshold, benchmark_ontology_term_count)

        output_df = output_df.append({'threshold': threshold, 'ontology': ontology, 'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp, 'precision': pr, 'recall': rc, 'f1': f1, 'mcc': mcc, 'accuracy': acc, 'coverage': coverage}, ignore_index=True)
    

    return output_df


def organism_specific_merged_eval(ontology, thresholds, benchmark, organism, pred_list):
    columns = ['organism', 'threshold', 'ontology', 'tn', 'fp', 'fn', 'tp', 'precision', 'recall', 'f1', 'mcc', 'accuracy']
    output_df = pd.DataFrame(columns=columns)

    _benchmark = benchmark[ontology][organism].copy()
    benchmark_ontology_term_count = benchmark[ontology]['benchmark_ontology_term_count']
    organism_spec_pred_list = [pred[ontology][organism] for pred in pred_list if organism in pred[ontology].keys()]
    
    merged_preds = pred_merger(organism_spec_pred_list, ontology, overall=False)
    
    for threshold in thresholds:
        
        tn, fp, fn, tp, pr, rc, f1, mcc, acc, coverage = evaluator(_benchmark, merged_preds, threshold, benchmark_ontology_term_count)

        output_df = output_df.append({'organism': organism, 'threshold': threshold, 'ontology': ontology, 'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp, 'precision': pr, 'recall': rc, 'f1': f1, 'mcc': mcc, 'accuracy': acc, 'coverage': coverage}, ignore_index=True)

    return output_df