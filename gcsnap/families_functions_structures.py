from gcsnap.sequence_mapping_online import SequenceMappingOnline
from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.apis import SwissProtAPI, AlphaFoldAPI, EbiAPI

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import split_dict_chunks

class FamiliesFunctionsStructures:
    def __init__(self, config: Configuration, gc: GenomicContext):
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.get_pdb = config.arguments['get_pdb']['value']
        self.get_annotations = config.arguments['get_functional_annotations']['value']

        # set parameters
        self.families = gc.get_families()

        self.console = RichConsole()

    def get_annotations_and_structures(self) -> dict:
        return self.annotations_and_structures
    
    def run(self) -> None:
        if self.get_pdb or self.get_annotations:
            # find all family members
            # always exclude non-conserved families and pseudogenes as uniprot mapping
            all_members = [member_code for family in self.families.keys() 
                           for member_code in self.families[family]['members']
                           if 0 < family < 10000]
            
            # map all members to UniProtKB-AC
            mapping = SequenceMappingOnline(self.config, all_members, 
                                            'UniProtKB-AC', 'family members')
            mapping.run()
            mapping_dict = mapping.get_target_to_result_dict()
            mapping.log_failed()

            # split the dictionary into chunks
            dict_list = split_dict_chunks(self.families, self.cores)
            # create parallel arguments
            parallel_args = [(dict_, mapping_dict, self.get_pdb, self.get_annotations) 
                             for dict_ in dict_list]

            with self.console.status('Get functional annotations and structures'):
                result_list = processpool_wrapper(self.cores, parallel_args, self.run_each)
                # combine results
                self.annotations_and_structures = {k: v for dict_ in result_list for k, v in dict_.items()}
        else:
            self.annotations_and_structures = {}
            self.console.print_skipped_step('No annotations or structures retrieval requested')

    def run_each(self, args: tuple) -> dict:
        families, mapping_dict, get_pdb, get_annotations = args

        # all families that are not pseudogenes
        # a list of tuples (family, member)
        family_keys = [(key, sorted(families[key]['members'])) 
                        for key in sorted(list(families.keys()))
                        if 'function' not in families[key]]
        
        family_function = ''
        family_structure = ''
        family_uniprot_code = ''
        family_model_state = 'Model does not exist'

        for family, members in family_keys:

            # pseudo genes and non-conserved families
            if family == 0 or family == 10000:
                family_function = 'n.a.'
                family_structure = 'n.a.'
                family_uniprot_code = 'n.a.'
                family_model_state = 'n.a.'

            # conserved families
            else:
                family_function = ''
                family_structure = ''
                family_uniprot_code = ''
                family_model_state = 'Model does not exist'

                # it is done until a result is found
                for member in members:
                    # if not found, it will be nan
                    uniprot_code = mapping_dict.get(member,'nan')

                    if (family_uniprot_code == '' or family_structure == '' or 
                        (family_function == '' or family_function['Function_description'] == '')):
                    
                        if get_pdb and family_structure == '':
                            # check either in SwissModel or AlphaFold
                            curr_pdb = SwissProtAPI.find_uniprot_in_swissmodel(uniprot_code)
                            if 'nan' in curr_pdb:
                                curr_pdb = AlphaFoldAPI.find_uniprot_in_alphafold_database(uniprot_code)

                            if 'nan' not in curr_pdb:
                                family_structure = curr_pdb
                                family_uniprot_code = uniprot_code
                                family_model_state = 'Model exists'
                            elif uniprot_code != 'nan' and curr_pdb == 'nan*':
                                family_uniprot_code = uniprot_code

                        if (get_annotations and (family_function == '' 
                                                 or family_function['Function_description'] == '')):        
                            # get functional annotations
                            curr_uniprot_annotations = EbiAPI.get_uniprot_annotations(uniprot_code, 
                                                                    previous_annotations = family_function)
                            if curr_uniprot_annotations != 'nan':
                                family_function = curr_uniprot_annotations

                            if uniprot_code != 'nan' and family_structure == '' and family_uniprot_code == '':
                                family_uniprot_code = uniprot_code

                if family_uniprot_code == '':
                    family_model_state = 'Not possible to map'

            families[family]['function'] = family_function
            families[family]['structure'] = family_structure
            families[family]['uniprot_code'] = family_uniprot_code
            families[family]['model_state'] = family_model_state

        return families