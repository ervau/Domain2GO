import requests
from io import StringIO
from Bio import SeqIO
import os
import time
import pandas as pd

import argparse

parser = argparse.ArgumentParser(description='Protein function prediction based on Domain2GO mappings')
parser.add_argument('--mapping_path', required=True, help='Path to previously generated Domain2GO mappings or where to save new mappings if not previously generated')

args = parser.parse_args()


def find_domains():

    email = input("Please enter your email: ")
    protein_input = input("Please enter the protein sequence or fasta file location: ")
    print(f'Finding domains in sequence using InterProScan...')

    # read protein sequence
    # if '.fa' or '.fasta' in protein_input read fasta file

    if '.fa' in protein_input or '.fasta' in protein_input:
        # read fasta file
        fasta_sequences = SeqIO.parse(open(protein_input),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
    
    else:
        # read protein sequence
        sequence = protein_input

    # send request to interproscan api
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
        'Accept': 'text/plain',
    }

    data= {
       'email': email,
       'stype': 'p',
       'sequence': f'{sequence}'}


    job_id_response = requests.post('https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run', headers=headers, data=data)
    job_id = job_id_response.text

    print(f'InterProScan job id: {job_id}')

    # get results

    headers = {
    'Accept': 'application/json',
    }

    job_result_url = f'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json'
    job_result_response = requests.get(job_result_url, headers=headers)

    json_output = None
    if job_result_response.status_code == 200:
        print('InterProScan job done')
        json_output= job_result_response.json()['results']    
    elif job_result_response.status_code == 400:
        # wait for 30 seconds and try again
        time.sleep(60)

        job_result_response = requests.get(job_result_url, headers=headers)
        
        if job_result_response.status_code == 200:
            print('InterProScan job done')
            json_output= job_result_response.json()['results']    

        else:
            time.sleep(60)
            job_result_response = requests.get(job_result_url, headers=headers)
            if job_result_response.status_code == 200:
                print('InterProScan job done')
                json_output= job_result_response.json()['results']    
            else:
                print('Error in InterProScan job')
                print(job_result_response.text)

    entries = set()
    if json_output is not None:
        for elem in json_output['matches']:
            locations = {f"{i['start']}, {i['end']}" for i in elem['locations']}
            entry = elem['signature']['entry']
            if type(entry) == dict and entry['type'] == 'DOMAIN':
                entry_dict = {
                    'accession': entry['accession'],
                    'name': entry['name'],
                    'locations': locations
                }

                entries.add(entry_dict)
        
    if entries:
        print('Domains found')

        # create domains dataframe
        domains_df = pd.DataFrame(entries)

        return domains_df

    else:
        print('No domains found')
        return None
                        
    # generate protein function predictions based on domain2go mappings

def generate_function_predictions(domains_df, mapping_path):
    
    # read domain2go mappings
    domain2go_df = pd.read_csv(os.path.join(mapping_path, 'finalized_domain2go_mappings.txt'), sep='\t', header=None)

    # merge domain2go mappings with domains found in protein sequence
    merged_df = pd.merge(domains_df, domain2go_df, how='left', left_on='accession', right_on='Interpro')

    merged_df = merged_df[['accession', 'name', 'locations', 'GO', 's']]
    merged_df.columns = ['domain_accession', 'domain_name', 'domain_locations', 'go_id', 'probability']

    # save protein function predictions
    merged_df.to_csv(os.path.join(mapping_path, f'{uniprot_id}_function_predictions.txt'), sep='\t', index=False)


def main(args):

    domains_df = find_domains()

    if domains_df:
        # check if domain2go mappings exist
        domain2go_df_path = os.path.join(args.mapping_path, 'finalized_domain2go_mappings.txt')
        if not os.path.exists(domain2go_df_path):
            print('Domain2GO mappings not found, generating mappings')

            try:
                os.system(f'python main_training.py --em skip --enrichment skip --cafa_eval skip')
                generate_function_predictions(domains_df, args.mapping_path)

            except:
                print('Error in generating Domain2GO mappings')
                return None

        else:
            print(f'Domain2GO mappings found at {domain2go_df_path}, generating protein function predictions')
            generate_function_predictions(domains_df, mapping_path)


if __name__ == '__main__':
    main(args)