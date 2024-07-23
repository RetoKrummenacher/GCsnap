import requests
import json

class SwissProtAPI:
    @staticmethod
    def find_uniprot_in_swissmodel(uniprot_code: str) -> str:
        """
        Find the link for information about a protein defined by its uniprot code in the SwissModel database.

        Args:
            uniprot_code (str): The uniprot code of the protein.

        Returns:
            str: The link to the SwissModel database entry of the protein.
        """        
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
        """
        Find the link to the 3D model of a protein defined by its uniprot code in the AlphaFold database.

        Args:
            uniprot_code (str): The uniprot code of the protein.

        Returns:
            str: The link to the 3D model of the protein in the AlphaFold database.
        """        
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
    def get_uniprot_annotations(uniprot_code: str, previous_annotations: dict = {}) -> dict:
        """
        Get annotations for a protein defined by its uniprot code from the EBI database.

        Args:
            uniprot_code (str): The uniprot code of the protein.
            previous_annotations (dict, optional): Previously fond annotations of the protein. Defaults to {}.

        Returns:
            dict: The annotations of the protein. 
        """
        uniprot_annotations = 'nan'
        if uniprot_code != 'nan':
            uniprot_accession = uniprot_code.split('_')[0]
            uniprot_link = 'https://www.ebi.ac.uk/proteins/api/proteins/{}'.format(uniprot_accession)
            try:
                uniprot_req = requests.get(uniprot_link, headers={ "Accept" : "application/json"})
                if uniprot_req.ok:
                    uniprot_data = uniprot_req.text
                    uniprot_data = json.loads(uniprot_data)
                    uniprot_annotations = EbiAPI.parse_uniprot_data(uniprot_data, 
                                                        previous_annotations = previous_annotations)
            except:
                uniprot_annotations = 'nan'        
        return uniprot_annotations
    
    @staticmethod
    def get_uniprot_annotations_batch(uniprot_codes: list, with_parsing: bool = False) -> dict:
        """
        Get annotations for a list of proteins defined by their uniprot codes from the EBI database.
        Searches in batches. Limit of up to 100 accessions per batch is not checked.

        Args:
            uniprot_codes (list): The list of uniprot codes of the proteins.
            with_parsing (bool, optional): Should the retrieved data be parsed. Defaults to False.

        Returns:
            dict: The annotations of the proteins, either parsed or not.
        """        
        uniprot_annotations = {}
        uniprot_accessions = [code.split('_')[0] for code in uniprot_codes if code != 'nan']
        # separator is %2C: https://www.ebi.ac.uk/proteins/api/doc/#!/proteins/search
        # the maximum number of accessions per is 100
        uniprot_accession_str = '%2C'.join(uniprot_accessions)
        
        if uniprot_accession_str:
            uniprot_link = 'https://www.ebi.ac.uk/proteins/api/proteins?accession={}'.format(uniprot_accession_str)
            try:
                # returns a list of results
                uniprot_req = requests.get(uniprot_link, headers={ "Accept" : "application/json"})
                if uniprot_req.ok:
                    uniprot_data = uniprot_req.json()
                    if with_parsing:
                        for data in uniprot_data:
                            accession = data.get('accession')
                            uniprot_annotations[accession] = EbiAPI.parse_uniprot_data(data)
                    else:
                        for data in uniprot_data:
                            # combine all data without parsing
                            accession = data.get('accession')
                            uniprot_annotations[accession] = data
                else:
                    uniprot_annotations = {accession: 'nan' for accession in uniprot_accessions}
            except:
                uniprot_annotations = {accession: 'nan' for accession in uniprot_accessions}
        else:
            uniprot_annotations = {accession: 'nan' for accession in uniprot_accessions}
        
        return uniprot_annotations    
    
    @staticmethod
    def parse_uniprot_data(uniprot_data: dict, previous_annotations: dict = {}) -> dict:
        """
        Parse the data of a protein from the EBI database and extract annotations.

        Args:
            uniprot_data (dict): The data of the protein.
            previous_annotations (dict, optional): Previously found annotations of the protein. Defaults to {}.

        Returns:
            dict: The annotations of the protein.
        """        
        if previous_annotations == {}:
            uniprot_annotations = {'TM_topology': '', 'GO_terms': [], 'Keywords': [],
                                   'Function_description': ''}
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
                if (keyword_term not in uniprot_annotations['Keywords'] 
                    and keyword_term != 'Reference proteome'):
                    uniprot_annotations['Keywords'].append(keyword_term)

        if 'comments' in uniprot_data:
            for comment in uniprot_data['comments']:
                if comment['type'] == 'FUNCTION' and uniprot_annotations['Function_description'] == '':
                    uniprot_annotations['Function_description'] = comment['text'][0]['value']

        return uniprot_annotations    