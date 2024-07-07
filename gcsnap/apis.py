import requests
import json

class SwissProtAPI:
    @staticmethod
    def find_uniprot_in_swissmodel(uniprot_code: str) -> str:
        link = 'nan'
        if uniprot_code != 'nan':
            link = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
            json_link = json_link = '{}.json'.format(link)
            try:
                swissmodel_req = requests.get(json_link)
                if swissmodel_req.ok:
                    swiss_repository_data = swissmodel_req.text
                    swiss_repository_data = json.loads(swiss_repository_data)
                    if len(swiss_repository_data['result']) == 0 or len(swiss_repository_data['result']['uniprot_entries']) == 0:
                        link = 'nan*'
                    elif len(swiss_repository_data['result']['structures']) == 0:
                        link = 'nan'
                else:
                    link = 'nan*'
            except:
                link = 'nan*'
        return link
    

class AlphaFoldAPI:    
    @staticmethod
    def find_uniprot_in_alphafold_database(uniprot_code: str) -> str:
        link = 'nan'
        if uniprot_code != 'nan':
            link = 'https://alphafold.ebi.ac.uk/entry/{}'.format(uniprot_code)
            protein_link = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v3.pdb'.format(uniprot_code)
            try:
                afdb_req = requests.get(protein_link)
                if not afdb_req.ok:
                    link = 'nan*'
                else:
                    link = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
            except:
                link = 'nan*'
        return link
    
    
class EbiAPI:
    @staticmethod
    def get_uniprot_annotations(uniprot_code: str, previous_annotations: str = '') -> str:

        uniprot_annotations = 'nan'
        if uniprot_code != 'nan':
            uniprot_accession = uniprot_code.split('_')[0]
            uniprot_link = 'https://www.ebi.ac.uk/proteins/api/proteins/{}'.format(uniprot_accession)
            try:
                uniprot_req = requests.get(uniprot_link, headers={ "Accept" : "application/json"})
                if uniprot_req.ok:
                    uniprot_data = uniprot_req.text
                    uniprot_data = json.loads(uniprot_data)
                    uniprot_annotations = EbiAPI.parse_uniprot_data(uniprot_data, previous_annotations = previous_annotations)
            except:
                uniprot_annotations = 'nan'        
        return uniprot_annotations
    
    @staticmethod
    def parse_uniprot_data(uniprot_data, previous_annotations = ''):
        if previous_annotations == '':
            uniprot_annotations = {'TM_topology': '', 'GO_terms': [], 'Keywords': [], 'Function_description': ''}
        else:
            uniprot_annotations = previous_annotations

        if 'features' in uniprot_data:
            for feature in uniprot_data['features']:
                if feature['type'] == 'TRANSMEM' and uniprot_annotations['TM_topology'] == '':
                    tm_topology = feature['description']

                    if 'Helical' in tm_topology:
                        tm_topology = 'alpha-helical'
                    elif 'Beta' in tm_topology:
                        tm_topology = 'beta-stranded'
                    else:
                        tm_topology = 'transmembrane'

                    uniprot_annotations['TM_topology'] = tm_topology

        if 'dbReferences' in uniprot_data:
            for dbReference in uniprot_data['dbReferences']:
                if dbReference['type'] == 'GO':
                    go_term = dbReference['properties']['term']
                    if go_term not in uniprot_annotations['GO_terms']:
                        uniprot_annotations['GO_terms'].append(go_term)

        if 'keywords' in uniprot_data:
            for keyword in uniprot_data['keywords']:
                keyword_term = keyword['value']
                if keyword_term not in uniprot_annotations['Keywords'] and keyword_term != 'Reference proteome':
                    uniprot_annotations['Keywords'].append(keyword_term)

        if 'comments' in uniprot_data:
            for comment in uniprot_data['comments']:
                if comment['type'] == 'FUNCTION' and uniprot_annotations['Function_description'] == '':
                    uniprot_annotations['Function_description'] = comment['text'][0]['value']

        return uniprot_annotations    