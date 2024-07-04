import os
import json


class Syntenies:
    def __init__(self):
        self.all = {}

    def update(self, input_dict: dict) -> None:
        for k,v in input_dict.items():
            current = self.all.get(k, {})
            # merge dictionaries (|= in python 3.9+ in place merge)
            current |= v
            self.all[k] = current

    def get_all_ncbi_codes(self) -> list:
        return [code for sub_dict in self.all.values() 
                for code in sub_dict['flanking_genes'].get('ncbi_codes', [])]
    
    def get_key_value_list(self) -> list[tuple]:
        return [(target, target_dict) for target, target_dict in self.all.items()]       

    def write_to_file(self, file_name: str, file_path: str = None, msg: str = None) -> None:
        # using os.getcwd() as default path does not work, as its evaluated when the function is defined
        if file_path is None:
            file_path = os.getcwd()         
        with open(os.path.join(file_path, file_name), 'w') as file:
            json.dump(self.all, file, indent = 4)
        # print and log message
        if msg is not None:
            self.console.print_done(msg)               

    def read_from_file(self, file_name: str, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.getcwd()        
        with open(os.path.join(file_path, file_name), 'r') as file:
            self.all = json.load(file)    

    def write_to_fasta(self, file_name: str, file_path: str = None, exclude_pseudogenes: bool = False) -> str:
        # extract needed information from syntenies
        # There should be no performance difference between using dict.get() and dict[key]
        # at least when no default value is provided        
        lines_to_write = ['>{}|{}\n{}'.format(ncbi_code, 
                                self.all[target]['flanking_genes']['names'][i], 
                                self.all[target]['flanking_genes']['sequences'][i])
            for target in self.all.keys()
            for i, ncbi_code in enumerate(self.all[target]['flanking_genes']['ncbi_codes'])
            if self.all[target]['flanking_genes']['names'][i] != 'pseudogene' 
            or not exclude_pseudogenes]
        
        # write to fasta file
        if file_path is None:
            file_path = os.getcwd() 
        fasta_file = os.path.join(file_path, file_name)
        with open(fasta_file, 'w') as file:
            file.write('\n'.join(lines_to_write))

        return fasta_file

    # TODO: Seqeuence length is only needed for BLAST, but not for MMseqs
    # As we dropped BLAST support, this is redundant
    # def get_sequence_length(self) -> dict:
    #     return {ncbi_code : len(self.syntenies[target]['flanking_genes']['sequences'][i])
    #             for target in self.syntenies.keys()
    #             for i, ncbi_code in enumerate(self.syntenies[target]['flanking_genes']['ncbi_codes'])}                  

    @staticmethod
    def get_empty_flanking_genes() -> dict:
        return {'relative_starts' : [],
                'relative_ends' : [],
                'ncbi_codes': [],
                'starts': [],
                'ends': [],
                'directions': [],
                'names': [],
                'sequences': [],
                'species': None,
                'taxID': None,
                'families': [],
                'TM_annotations': []}