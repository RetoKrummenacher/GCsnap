import os
import json
import numpy as np
import pandas as pd
import statistics

from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration

class GenomicContext:
    """ 
    Methods and attributes to store and manipulate the genomic context information.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        out_label (str): The label of the output.
        syntenies (dict): The dictionary with the syntenies of the target genes.
        families (dict): The dictionary with the families assigned to the flanking genes.
        operon_types_summary (dict): The dictionary with the operon types summary.
        selected_operons (dict): The dictionary with the selected operons.
        most_populated_operon (str): The most populated operon type.
        taxonomy (dict): The dictionary with the taxonomy information.
        curr_targets (list): The list of current targets.
        n_max_operons (int): The maximum number of operons to select.
        current_targets (list): The list of current targets.
    """

    def __init__(self, config: Configuration, out_label: str = '') -> None:
        """
        Initialize the GenomicContext object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            out_label (str, optional): The label of the output. Defaults to ''.
        """        
        # empty attributes
        self.syntenies = {}
        self.families = {}
        self.operon_types_summary = {}
        self.selected_operons = {}
        self.most_populated_operon = ''    
        self.taxonomy = {}  
        self.curr_targets = []  
        # get config arguments
        self.n_max_operons = config.arguments['n_max_operons']['value']

        # set parameters
        self.out_label = out_label
        self.config = config

        self.console = RichConsole()

    def update_syntenies(self, input_dict: dict) -> None:
        """
        Update the syntenies attribute with new information.

        Args:
            input_dict (dict): The dictionary with the new syntenies information.
        """        
        for k,v in input_dict.items():
            current = self.syntenies.get(k, {})
            # merge dictionaries (|= in python 3.9+ in place merge)
            current |= v
            self.syntenies[k] = current

    def update_families(self, input_dict: dict) -> None:
        """
        Update the families attribute with new information.

        Args:
            input_dict (dict): The dictionary with the new families information.
        """        
        for k,v in input_dict.items():
            current = self.families.get(k, {})
            # merge dictionaries (|= in python 3.9+ in place merge)
            current |= v
            self.families[k] = current    

    def update_taxonomy(self, input_dict: dict) -> None: 
        """
        Update the taxonomy attribute with new information.

        Args:
            input_dict (dict): The dictionary with the new taxonomy information.
        """        
        for k,v in input_dict.items():
            current = self.taxonomy.get(k, {})
            # merge dictionaries (|= in python 3.9+ in place merge)
            current |= v
            self.taxonomy[k] = current                  

    def get_curr_targets(self) -> list:
        """
        Getter for the curr_targets attribute.

        Returns:
            list: The list of current targets.
        """        
        return self.curr_targets
    
    def get_families(self) -> dict:
        """
        Getter for the families attribute.

        Returns:
            dict: The dictionary with the families assigned to the flanking genes.
        """        
        return self.families
    
    def get_syntenies(self) -> dict:
        """
        Getter for the syntenies attribute.

        Returns:
            dict: The dictionary with the syntenies of the target genes.
        """        
        return self.syntenies
    
    def get_operon_types(self) -> dict:
        """
        Getter for the operon_types_summary attribute.

        Returns:
            dict: The dictionary with the operon types summary.
        """        
        return self.operon_types_summary
    
    def get_selected_operons(self) -> dict:
        """
        Getter for the selected_operons attribute.

        Returns:
            dict: The dictionary with the selected operons.
        """        
        return self.selected_operons
    
    def get_taxonomy(self) -> dict:
        """
        Getter for the taxonomy attribute.

        Returns:
            dict: The dictionary with the taxonomy information.
        """        
        return self.taxonomy
    
    def get_most_populated_operon(self) -> str:
        """
        Getter for the most_populated_operon attribute.

        Returns:
            str: The most populated operon type.
        """        
        return self.most_populated_operon

    def get_all_ncbi_codes(self) -> list:
        """
        Get all NCBI codes from the syntenies attribute.

        Returns:
            list: The list of all NCBI codes.
        """        
        return [code for sub_dict in self.syntenies.values() 
                for code in sub_dict['flanking_genes'].get('ncbi_codes', [])]
    
    def get_all_taxids(self) -> list:
        """
        Get all taxIDs from the syntenies attribute.

        Returns:
            list: The list of all taxIDs.
        """        
        return [sub_dict['flanking_genes']['taxID'] for sub_dict in self.syntenies.values()
                if sub_dict['flanking_genes']['taxID'] is not None]
    
    def get_syntenies_key_value_list(self) -> list:
        """
        Get the syntenies attribute as a list of key-value pairs.

        Returns:
            list: The list of key-value pairs.
        """        
        return [(k, v) for k, v in self.syntenies.items()]
    
    def create_and_write_families_summary(self) -> None:
        """
        Create and write the families summary to a text file.
        """        
        self.create_families_summary()
        self.write_families_summary_to_txt()

    def create_and_write_operon_types_summary(self) -> None:
        """
        Create and write the operon types summary to a text file.
        """        
        self.create_operon_types_summary()
        self.write_operon_types_summary_to_txt()

    def write_taxonomy_to_json(self, file_name: str, file_path: str = None) -> None:
        """
        Write the taxonomy information to a json file.

        Args:
            file_name (str): The name of the file to write.
            file_path (str, optional): The path of the file to write. Defaults to None using os.getcwd().
        """        
        if file_path is None:
            file_path = os.getcwd()         
        with open(file_name, 'w') as file:
            json.dump(self.taxonomy, file, indent = 4)   
        # print and log message
        self.console.print_done('Taxonomy information written to {}'.format(file_name))        
    
    def write_syntenies_to_json(self, file_name: str, file_path: str = None) -> None:
        """
        Write the syntenies information to a json file.

        Args:
            file_name (str): The name of the file to write.
            file_path (str, optional): The path of the file to write. Defaults to None using os.getcwd().
        """        
        # using os.getcwd() as default path does not work, as its evaluated when the function is defined
        if file_path is None:
            file_path = os.getcwd()         
        with open(os.path.join(file_path, file_name), 'w') as file:
            json.dump(self.syntenies, file, indent = 4)
        # print and log message
        self.console.print_done('Syntenies information written to {}'.format(file_name))      

    def write_families_to_json(self, file_name: str, file_path: str = None) -> None:
        """
        Write the families information to a json file.

        Args:
            file_name (str): The name of the file to write.
            file_path (str, optional): The path of the file to write. Defaults to None using os.getcwd().
        """        
        if file_path is None:
            file_path = os.getcwd()   
        with open(os.path.join(file_path, file_name), 'w') as file:
            json.dump(self.families, file, indent = 4)   
        # print and log message
        self.console.print_done('Families information written to {}'.format(file_name))     

    def write_selected_operons_to_json(self, file_name: str, file_path: str = None) -> None:
        """
        Write the selected operons to a json file.

        Args:
            file_name (str): The name of the file to write.
            file_path (str, optional): The path of the file to write. Defaults to None using os.getcwd().
        """        
        if file_path is None:
            file_path = os.getcwd()   
        with open(os.path.join(file_path, file_name), 'w') as file:
            json.dump(self.selected_operons, file, indent = 4)   
        # print and log message
        self.console.print_done('Selected operon written to {}'.format(file_name))     

    def write_to_fasta(self, file_name: str, file_path: str = None, exclude_pseudogenes: bool = False) -> str:
        """
        Write the flanking genes to a fasta file.

        Args:
            file_name (str): The name of the file to write.
            file_path (str, optional): The path of the file to write. Defaults to None using os.getcwd().
            exclude_pseudogenes (bool, optional): Whether to exclude pseudogenes. Defaults to False.

        Returns:
            str: _description_
        """        
        # extract needed information from syntenies
        # There should be no performance difference between using dict.get() and dict[key]
        # at least when no default value is provided        
        lines_to_write = ['>{}|{}\n{}\n'.format(ncbi_code, 
                                self.syntenies[target]['flanking_genes']['names'][i], 
                                self.syntenies[target]['flanking_genes']['sequences'][i])
            for target in self.syntenies.keys()
            for i, ncbi_code in enumerate(self.syntenies[target]['flanking_genes']['ncbi_codes'])
            if self.syntenies[target]['flanking_genes']['names'][i] != 'pseudogene' 
            or not exclude_pseudogenes]
        
        # write to fasta file
        if file_path is None:
            file_path = os.getcwd() 
        fasta_file = os.path.join(file_path, file_name)
        with open(fasta_file, 'w') as file:
            file.writelines(lines_to_write)

        return fasta_file

    def read_syntenies_from_json(self, file_name: str, file_path: str = None) -> None:
        """
        Read the syntenies information from a json file.

        Args:
            file_name (str): The name of the file to read.
            file_path (str, optional): The path of the file to read. Defaults to None using os.getcwd().
        """        
        if file_path is None:
            file_path = os.getcwd()        
        with open(os.path.join(file_path, file_name), 'r') as file:
            self.syntenies = json.load(file)  

    def read_families_from_json(self, file_name: str, file_path: str = None) -> None:
        """
        Read the families information from a json file.

        Args:
            file_name (str): The name of the file to read.
            file_path (str, optional): The path of the file to read. Defaults to None using os.getcwd().
        """        
        if file_path is None:
            file_path = os.getcwd()        
        with open(os.path.join(file_path, file_name), 'r') as file:
            self.families = json.load(file)              

    def get_fasta_order(self, exclude_pseudogenes: bool = False) -> list:
        """
        Get the order of the ncbi codes in the fasta file.

        Args:
            exclude_pseudogenes (bool, optional): Whether to exclude pseudogenes. Defaults to False.

        Returns:
            list: The list of ncbi codes in the same order as in the fasta file.
        """        
        # the same order as in the fasta file
        return [ncbi_code for target in self.syntenies.keys()
            for i, ncbi_code in enumerate(self.syntenies[target]['flanking_genes']['ncbi_codes'])
            if self.syntenies[target]['flanking_genes']['names'][i] != 'pseudogene' 
            or not exclude_pseudogenes]  

    def create_families_summary(self) -> None: 
        """
        Create the families summary output.
        """          
        with self.console.status('Create families summary'):
            for target in self.syntenies.keys():
                for i, family in enumerate(self.syntenies[target]['flanking_genes']['families']):
                    curr_ncbi_code = self.syntenies[target]['flanking_genes']['ncbi_codes'][i]

                    if family not in self.families:
                        self.families[family] = {'name': [], 'members': [], 'all_names': []}

                    if curr_ncbi_code not in self.families[family]['members']:

                        if family != 0:
                            self.families[family]['all_names'].append(self.syntenies[target]['flanking_genes']['names'][i])
                        else:
                            self.families[family]['all_names'] = ['Non-conserved']

                        self.families[family]['members'].append(curr_ncbi_code)

            for family in self.families.keys():
                if (len(set(self.families[family]['all_names'])) > 1 
                    and 'hypothetical protein' in set(self.families[family]['all_names'])):
                    self.families[family]['name'] = [name for name in self.families[family]['all_names'] 
                                                     if name != 'hypothetical protein']
                else:
                    self.families[family]['name'] = self.families[family]['all_names']
                
                try:
                    self.families[family]['name'] = statistics.mode(self.families[family]['name'])
                except:
                    self.families[family]['name'] = self.families[family]['name'][0]

            if -1 in self.families:
                n_pseudogenes = len(self.families[-1]['members'])
            else:
                n_pseudogenes = 0

            if 0 in self.families:
                n_nonconserved = len(self.families[0]['members'])
            else:
                n_nonconserved = 0

            # make it sorted
            self.families = dict(sorted(self.families.items()))

        msg = 'Found {} conserved protein families, {} pseudogenes and {} non-conserved protein coding regions'.format(
            len([i for i in self.families if i > 0]), n_pseudogenes, n_nonconserved)
        self.console.print_info(msg)

    def write_families_summary_to_txt(self) -> None:
        """
        Write the families summary to a text file.
        """        
        out_file = '{}_protein_families_summary.txt'.format(self.out_label)
        with open(out_file, 'w') as file:
            for family in sorted(list(self.families.keys())):
                if self.families[family]['name'] != 'Non-conserved':
                    file.write('\n ### Family: {} -> {}\n\n'.format(family, self.families[family]['name']))
                    for i, member in enumerate(self.families[family]['members']):
                        file.write('	 {}\t{}\n'.format(member, self.families[family]['all_names'][i]))

    def create_operon_types_summary(self) -> None:
        """
        Create the operon types summary output.
        """        
        with self.console.status('Create operons summary'):
            advanced = False
            if 'operon_PaCMAP' in self.syntenies[list(self.syntenies.keys())[0]]:
                advanced = True

            for target in self.syntenies:
                curr_operon_type = self.syntenies[target]['operon_type']
                if curr_operon_type not in self.operon_types_summary:
                    self.operon_types_summary[curr_operon_type] = {'target_members': [], 'operon_protein_families_structure': []}
                    if advanced:
                        self.operon_types_summary[curr_operon_type]['operon_PaCMAP'] = []
                        self.operon_types_summary[curr_operon_type]['operon_filtered_PaCMAP'] = []


                self.operon_types_summary[curr_operon_type]['target_members'].append(target)
                self.operon_types_summary[curr_operon_type]['operon_protein_families_structure'].append(
                    self.syntenies[target]['flanking_genes']['families'])
                if advanced:
                    self.operon_types_summary[curr_operon_type]['operon_PaCMAP'].append(
                        self.syntenies[target]['operon_PaCMAP'])
                    self.operon_types_summary[curr_operon_type]['operon_filtered_PaCMAP'].append(
                        self.syntenies[target]['operon_filtered_PaCMAP'])

            if advanced:
                for curr_operon_type in self.operon_types_summary:
                    centroid_coords = np.mean(self.operon_types_summary[curr_operon_type]
                                            ['operon_filtered_PaCMAP'], axis=0)
                    self.operon_types_summary[curr_operon_type]['operon_centroid_PaCMAP'] = list(centroid_coords) 

        msg = 'Found {} operon types (out of a total of {} input targets)'.format(
            len(self.operon_types_summary), len(self.syntenies))    
        self.console.print_info(msg)             

    def write_operon_types_summary_to_txt(self) -> None:
        """
        Write the operon types summary to a text file.
        """        
        out_file = '{}_operon_types_summary.txt'.format(self.out_label)
        with open(out_file, 'w') as file:
            for operon_type in self.operon_types_summary:
                file.write('\n ### Operon type: {}\n\n'.format(operon_type))
                for i, target in enumerate(self.operon_types_summary[operon_type]['target_members']):
                    file.write('	 {}\t{}\n'.format(target, 
                                self.operon_types_summary[operon_type]['operon_protein_families_structure'][i])) 

    def find_most_populated_operon_types(self) -> None:
        """
        Find the most populated operon types.
        """        
        with self.console.status('Find most populated operon types'):        
            operons_count_matrix = []
            for operon in self.operon_types_summary:
                operons_count_matrix.append([operon, len(
                    self.operon_types_summary[operon]['target_members'])])

            operons_count_matrix = pd.DataFrame(operons_count_matrix)
            operons_count_matrix = operons_count_matrix.sort_values(by = [1, 0], ascending = [False, True])	
            operons_count_matrix = np.array(operons_count_matrix)
            
            if len(operons_count_matrix) > self.n_max_operons:
                operons_count_matrix = operons_count_matrix[:self.n_max_operons+1]

            for i, line in enumerate(operons_count_matrix):
                label = 'GC Type {:05d}'.format(line[0])
                if i == 0:
                    self.most_populated_operon = label
                
                self.selected_operons[label] = self.operon_types_summary[line[0]]

        msg = 'Selected {} operon/genomic_context types, with most populated corresponding to {}'.format(
            len(self.selected_operons), self.most_populated_operon)   
        self.console.print_info(msg)    

        # write selected operons to json
        self.write_selected_operons_to_json('selected_operons.json')

    def write_summary_table(self, file_name: str, file_path: str = None) -> None:
        """
        Write the summary table to a text file (.tab)

        Args:
            file_name (str): The name of the file to write.
            file_path (str, optional): The path of the file to write. Defaults to None using os.getcwd().
        """        
        if file_path is None:
            file_path = os.getcwd()    

        lines_to_write = []
        # header
        if 'TM_annotations' in self.syntenies[list(self.syntenies.keys())[0]]['flanking_genes']:
            header_line = '\t'.join(['Operon type', 
                                     'Target', 
                                     'AssemblyId', 
                                     'Gene direction', 
                                     'Gene start', 
                                     'Gene end',
                                     'Relative gene start',
                                     'Relative gene end',
                                     'Protein family code',
                                     'EntrzID',
                                     'Protein name',
                                     'Transmembrane/Signal peptide prediction',
                                     'Superkingdom',
                                     'Phylum',
                                     'Class',
                                     'Order',
                                     'Genus',
                                     'Species'])
        else:
            # same without TM annotations
            header_line = '\t'.join(['Operon type',
                                     'Target',
                                     'AssemblyId',
                                     'Gene direction',
                                     'Gene start',
                                     'Gene end',
                                     'Relative gene start',
                                     'Relative gene end',
                                     'Protein family code',
                                     'EntrzID',
                                     'Protein name',
                                     'Superkingdom',
                                     'Phylum',
                                     'Class',
                                     'Order',
                                     'Genus',
                                     'Species'])
        # additional new line after header
        lines_to_write.append(header_line + '\n')
            
        # all paths from taxonomy
        tax_search_dict = self.create_taxonomy_search_dict()   

        # all targets and the corresponding operon type
        targets_operon_list = [(target, operon.split()[-2]) for operon in self.selected_operons 
                        for target in self.selected_operons[operon]['target_members']]
        
        for target, operon_type in targets_operon_list:
            # line space between the targets
            lines_to_write.append('\n')
            for i, prot_name in enumerate(self.syntenies[target]['flanking_genes']['names']):
                list_to_join = [operon_type, 
                                target, 
                                self.syntenies[target]['assembly_id'][1], 
                                self.syntenies[target]['flanking_genes']['directions'][i], 
                                str(self.syntenies[target]['flanking_genes']['starts'][i]), 
                                str(self.syntenies[target]['flanking_genes']['ends'][i]),
                                str(self.syntenies[target]['flanking_genes']['relative_starts'][i]), 
                                str(self.syntenies[target]['flanking_genes']['relative_ends'][i]), 
                                str(self.syntenies[target]['flanking_genes']['families'][i]), 
                                self.syntenies[target]['flanking_genes']['ncbi_codes'][i], 
                                self.syntenies[target]['flanking_genes']['names'][i]]
                line_to_write = '\t'.join(list_to_join)
                if 'TM_annotations' in self.syntenies[target]['flanking_genes']:
                    line_to_write += '\t' + self.syntenies[target]['flanking_genes']['TM_annotations'][i]
                # add taxonomy information by searching the dictionary
                line_to_write += '\t' + '\t'.join(tax_search_dict.get(target)) + '\n'
                lines_to_write.append(line_to_write)

        summary_file = os.path.join(file_path, file_name)
        with open(summary_file, 'w') as file:
            file.writelines(lines_to_write)

    def create_taxonomy_search_dict(self) -> dict:
        """
        Create a dictionary to search for taxonomy information.
        The information is a flat list version of the hierarchical taxonomy dictionary.

        Returns:
            dict: The dictionary to search for taxonomy information.
        """        
        flat_taxonomy = self.flatten_taxonomy(self.taxonomy)
        # create dictionary
        return {member : tax_list[:-1] for tax_list in flat_taxonomy 
                for member in tax_list[-1]['target_members']}

    def flatten_taxonomy(self, taxonomy: dict, parent_keys: list = []) -> list[list,dict]:
        """ 
        Recursvie method to flatten a nested dictionary. 

        Args:
            taxonomy (dict): The taxonomy dictionary to flatten.
            parent_keys (list, optional): List of parents. Defaults to [].

        Returns:
            list[list,dict]: List of lists representing a path from the root of the taxonomy 
            to a leaf node, including all the keys along the way. The last element of each list
            is a dictionary with the key 'target_members'.
        """        
        flat_list = []
        for key, value in taxonomy.items():
            if isinstance(value, dict) and 'target_members' not in value:
                flat_list.extend(self.flatten_taxonomy(value, parent_keys + [key]))
            else:
                flat_list.append(parent_keys + [key, value])
        return flat_list        

    def copy(self) -> 'GenomicContext':
        new_gc = GenomicContext(self.config, self.out_label)
        new_gc.syntenies = self.syntenies.copy()
        new_gc.families = self.families.copy()
        new_gc.operon_types_summary = self.operon_types_summary.copy()
        new_gc.selected_operons = self.selected_operons.copy()
        new_gc.most_populated_operon = self.most_populated_operon
        new_gc.taxonomy = self.taxonomy.copy()
        new_gc.curr_targets = self.curr_targets.copy()
        return new_gc    
                    
    @staticmethod
    def get_empty_flanking_genes() -> dict:  
        """
        Get an empty dictionary for the flanking genes.

        Returns:
            dict: The empty dictionary for the flanking genes.
        """        
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
                'families': []}