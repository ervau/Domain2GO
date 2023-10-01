import requests
from io import StringIO
from Bio import SeqIO
import os
import time
import pandas as pd

def find_domains():

    email = input("Please enter your email for InterProScan query: ")
    protein_input = input("Please enter the protein sequence or fasta file location: ")

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
        name = input("Please enter a name for the protein sequence: ")

    # send request to interproscan api
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
        'Accept': 'text/plain',
    }

    print(f'Finding domains in sequence using InterProScan. This may take a while...')

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
    
    json_output = None
    entries = dict()
    with requests.Session() as s:
        # try 10 times if not successful print error
        c=0
        while c<10:
            job_result_response = s.get(job_result_url, headers=headers)
            if job_result_response.status_code == 200:
                json_output= job_result_response.json()['results'][0]
                print('InterProScan job done')
                break
            else:
                time.sleep(60)
                c+=1

    if json_output is None:
        print('Error in InterProScan job')
        print(job_result_response.text)
        return None
    
    else:
        for elem in json_output['matches']:
            entry = elem['signature']['entry']

            location_list = [f"{i['start']}-{i['end']}" for i in elem['locations']]

            if type(entry) == dict and entry['type'] == 'DOMAIN':
                if entry['accession'] not in entries:
                    entries[entry['accession']] = {
                        'name': entry['name'],
                        # add locations as a list
                        'locations': location_list
                    }

                else:
                    try:
                        entries[entry['accession']]['locations'].extend(location_list)
                    except AttributeError:
                        entries[entry['accession']]['locations'] = entries[entry['accession']]['locations'].split(' ')
                        entries[entry['accession']]['locations'] = [i for i in entries[entry['accession']]['locations'] if i]
                        entries[entry['accession']]['locations'].extend(location_list)

                entries[entry['accession']]['locations'] = list(set(entries[entry['accession']]['locations']))
                entries[entry['accession']]['locations'] = sorted([i.split('-') for i in entries[entry['accession']]['locations']], key=lambda x: (int(x[0]), int(x[1])))
                entries[entry['accession']]['locations'] = ['-'.join(i) for i in entries[entry['accession']]['locations']]
                # entries[entry['accession']]['locations'] = '|'.join(entries[entry['accession']]['locations'])
        
    if entries:
        print('Domains found')

        # create domains dataframe
        domains_df = pd.DataFrame.from_dict(entries, orient='index').reset_index()
        domains_df['protein_name'] = name
        domains_df = domains_df[['protein_name', 'index', 'name', 'locations']]
        domains_df.columns = ['protein_name', 'domain_accession', 'domain_name', 'domain_locations']
        return domains_df

    else:
        print('No domains found')
        return None
                        
    # generate protein function predictions based on domain2go mappings

def generate_function_predictions(domains_df, mapping_path):
    
    # read domain2go mappings
    domain2go_df = pd.read_csv(os.path.join(mapping_path, 'finalized_domain2go_mappings.txt'))

    # merge domain2go mappings with domains found in protein sequence
    merged_df = pd.merge(domains_df, domain2go_df, left_on='accession', right_on='Interpro')

    # if merged_df is empty return
    if merged_df.empty:
        print('No function predictions found')
        return None
    
    else:
        merged_df['protein_name'] = domains_df['protein_name'].iloc[0]
        merged_df = merged_df[['protein_name', 'GO', 'locations', 's', 'accession', 'name',]]
        merged_df.columns = ['protein_name', 'GO_ID', 'domain_locations', 'probability', 'domain_accession', 'domain_name',]

        # save protein function predictions
        protein_name = domains_df['protein_name'].iloc[0]
        merged_df.to_csv(os.path.join(mapping_path, f'{protein_name}_function_predictions.txt'), index=False)
        print(f'Protein function predictions saved at {os.path.join(mapping_path, f"{protein_name}_function_predictions.txt")}')


def main():

    cwd = os.getcwd()
    mapping_path = "{}Domain2GO/output".format(cwd.split("Domain2GO")[0])

    domains_df = find_domains()

    if type(domains_df) == pd.DataFrame:
        # check if domain2go mappings exist
        domain2go_df_path = os.path.join(mapping_path, 'finalized_domain2go_mappings.txt')
        if not os.path.exists(domain2go_df_path):
            print('Domain2GO mappings not found, generating mappings')

            try:
                os.system(f'python main_training.py --em skip --enrichment skip --cafa_eval skip')
                generate_function_predictions(domains_df, mapping_path)

            except Exception as e:
                print('Error in generating Domain2GO mappings:')
                # print traceback
                print(e)
                return None

        else:
            print(f'Domain2GO mappings found at {domain2go_df_path}, generating protein function predictions')
            generate_function_predictions(domains_df, mapping_path)


if __name__ == '__main__':
    main()