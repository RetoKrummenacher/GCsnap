import pandas as pd

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 
from gcsnap.utils import processpool_wrapper
from gcsnap.uniprot_dbs_dict import UniprotDict

from gcsnap.uniprot_api import submit_id_mapping
from gcsnap.uniprot_api import check_id_mapping_results_ready
from gcsnap.uniprot_api import get_id_mapping_results_link
from gcsnap.uniprot_api import get_id_mapping_results_search


class SequenceMappingOnline:
    def __init__(self, config: Configuration, target_list: list, to_type: str):
        self.config = config
        self.console = RichConsole()

        self.target_list = target_list
        self.to_type = to_type

        # extract needed values from config
        self.cores = self.config.arguments['n_cpu']['value']

    def run(self):
        self.create_uniprot_dbs_with_targets()
        parallel_args = self.create_parallel_input()
        mapping_df_list = processpool_wrapper(self.cores, parallel_args, self.do_mapping)
        # concatenate all mapping dfs to one
        self.mapping_df = pd.concat(mapping_df_list, ignore_index=True)
        # extract all results for self.to_type from that df
        self.create_result_list()
           
    def create_uniprot_dbs_with_targets(self) -> None:
        uniprot_dict = UniprotDict()
        uniprot_dict.assign_target_to_type(self.target_list)
        # fill uniprot_dbs with the targets
        self.uniprot_dict = uniprot_dict.uniprot_dbs

    def create_result_list(self):
        self.result_list = self.mapping_df[self.to_type].tolist()

    def create_parallel_input(self, api_limit: int = 10000) -> list[tuple]:
        parallel_input_args = []

        # 1. Create a parallel workload for each uniprot type
        # This assures that we can do all of the same type with one request
        for key, values in self.uniprot_dict.items():
            if len(values['targets']) != 0:
                from_type = key
                ids = values['targets']
                from_db = values['from_dbs']
                to_db = self.uniprot_dict[self.to_type]['to_dbs']
                
                parallel_input_args.append((ids, from_db, to_db, from_type, self.to_type))
                
        # 2. Split the ids into sublists of more than api_limit
        # This assures that we do not exceed the api limit
        refinded_parallel_input_args = []
        additions = []
        for parallel_arg in parallel_input_args:
            ids, from_db, to_db, from_type, _ = parallel_arg    

            if len(ids) > api_limit:
                ids_sublists += self.split_recursive(ids, api_limit)
                additions += [(sublist, from_db, to_db, from_type, self.to_type)
                                        for sublist in ids_sublists]
            else:
                refinded_parallel_input_args.append(parallel_arg)

        parallel_input_args = refinded_parallel_input_args + additions
        return parallel_input_args
    
    def split_recursive(self, ids: list, api_limit: int) -> list:
        # Base case: if the length of ids is less than or equal to api_limit, 
        # return a list containing ids
        if len(ids) <= api_limit:
            # the trick is to return [ids], a list of a list which we can iterate through
            return [ids]
        
        # Recursive case: split ids into two halves and recursively split each half
        mid = len(ids) // 2
        return self.split_recursive(ids[:mid], api_limit) + self.split_recursive(ids[mid:], api_limit)         

    def do_mapping(self, arg: tuple) -> pd.DataFrame:
        ids, from_db, to_db, from_type, to_type = arg
        
        mapping_api_results = self.run_mapping_api(ids, from_db, to_db)
        mapping_extracted = self.extract_mapping_from_api_results(mapping_api_results, from_type, to_type)
        mapping_df = self.create_mapping_df(mapping_extracted)
        
        return mapping_df

    def run_mapping_api(self, ids: list, from_db: str, to_db: str) -> list:  
        #ids, from_db, to_db, from_type, to_type = arg
        # call methods from gcsnap.uniprot_api.py
        job_id = submit_id_mapping(from_db=from_db, to_db=to_db, ids=ids)
        
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)['results']
            
        if isinstance(results, list):
            return results
        else:
            return None

    def extract_mapping_from_api_results(self, results: list, from_type:str , to_type:str) -> dict:
        # the UniProt Rest API returns different formats depending on what db is searched
        # for instance when mapping to EMBL-Genbank, the results is a list with 1 level dict:
            # [{'from': 'A0A0E3MFP2', 'to': 'AKA75089.1'},
            # {'from': 'A0A0E3MFP2', 'to': 'AKA77782.1'}]
        # mapping to UniProtKB returns a list with nested dict:
            # [{'from': 'P05067', 'to': {'entryType': 'UniProtKB reviewed (Swiss-Prot)', 'primaryAccession': 'P05067', 
            
        extracted = {'from_type': from_type,
                    'to_type': to_type,
                    'from_id': [],
                    'to_id': []}
                
        for result in results:
            extracted['from_id'].append(result['from'])
            if isinstance(results[0]['to'], dict):
                extracted['to_id'].append(result['to']['primaryAccession'])
            else:
                extracted['to_id'].append(result['to'])
            
        return extracted    
    
    def create_mapping_df(self, mapping_dict: dict) -> pd.DataFrame:        
        id_types = ['UniProtKB-AC','UniProtKB-ID','RefSeq','GeneID','Ensembl','UniParc','EMBL-CDS']
        mappings = {id_type: [] for id_type in id_types}
        
        from_type = mapping_dict['from_type']
        mappings[from_type] = mapping_dict['from_id']
        to_type = mapping_dict['to_type']
        mappings[to_type] = mapping_dict['to_id']
        
        # fill remaining cols with nan
        for id_type in list(set(id_types).difference(set([from_type,to_type]))):
            mappings[id_type] += ['nan'] * len(mapping_dict['from_id'])
                
        return pd.DataFrame.from_dict(mappings)
    
    def merge_mapping_dfs(self, mapping_df: pd.DataFrame, columns_to_merge: list = ['EMBL-CDS']) -> pd.DataFrame:
        # Merge DataFrames on 'UniProtKB-AC'
        merged_df = pd.merge(self.mapping_df, mapping_df[['UniProtKB-AC'] + columns_to_merge], 
                             on='UniProtKB-AC', how='left', suffixes=('', '_y'))
        # Update 'EMBL-CDS' in df1 with values from df2
        for column in columns_to_merge:
            merged_df[column] = merged_df[column].combine_first(merged_df[f'{column}_y'])
            # Drop the additional 'EMBL-CDS_y' column
            merged_df.drop(columns=[f'{column}_y'], inplace=True)
        self.mapping_df = merged_df.copy()



# import pickle
# from gcsnap.targets import Target 
# import os
# import sys

# # Add the directory containing gcsnap to the system path
# current_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
# sys.path.append(parent_dir)

# def test_picklability(instance):
#     for attr, value in instance.__dict__.items():
#         try:
#             pickle.dumps(value)
#             print(f"{attr} is picklable")
#         except Exception as e:
#             print(f"{attr} is not picklable: {e}")

# config = Configuration()
# config.parse_arguments()
# targets = Target(config)
# targets.run()
# for out_label in targets.targets_lists:
#     targets_list = targets.targets_lists[out_label]
#     seq_mapping_online = SequenceMappingOnline(config, targets_list, 'UniProtKB-AC')
#     test_picklability(seq_mapping_online)