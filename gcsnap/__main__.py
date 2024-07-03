import os
# TODO: For faster debugging
import json

from gcsnap.rich_console import RichConsole 
from gcsnap.configuration import Configuration 
from gcsnap.targets import Target 
from gcsnap.sequence_mapping_online import SequenceMappingOnline
from gcsnap.assemblies import Assemblies
from gcsnap.sequences import Sequences
from gcsnap.protein_families import ProteinFamilies

def main():

    starting_directory = os.getcwd()

    console = RichConsole()
    console.print_title()

    # 1. Parse configuration and arguments
    with console.status('Parsing CLI arguments and config.yaml'):
        config = Configuration()
        config.parse_arguments()

    # 2. parse targets
    with console.status('Parsing targets'):
        targets = Target(config)
        targets.run()

    # 3. Iterate over each target list
    for out_label in targets.targets_lists:

        # TODO: General crash resisance (much more difficult as no loop over ach target)
        # Right now, check existing files and check what is missing.

        # A. Create working directory
        # TODO: Change if out_label_suffix is still used
        working_dir = '{}/{}'.format(starting_directory, out_label)
        # if not os.path.isdir(working_dir):
        #     os.mkdir(working_dir)
        os.chdir(working_dir)

        # targets_list = targets.targets_lists[out_label]

        # # B. Map sequences to UniProtKB-AC and NCBI EMBL-CDS
        # # TODO: Change to whatever input from Joana regardin the mapping.
        # # a). Map all targets to UniProtKB-AC
        # mappingA = SequenceMappingOnline(config, targets_list, 'UniProtKB-AC')
        # mappingA.run()

        # # b) Map all to RefSeq
        # mappingB = SequenceMappingOnline(config, mappingA.get_codes(), 'RefSeq')
        # mappingB.run()
        # # merge them to A (only if A is not nan)
        # mappingA.merge_mapping_dfs(mappingB.mapping_df)

        # # c). Map all targets to NCBI EMBL-CDS
        # mappingC = SequenceMappingOnline(config, mappingA.get_codes(), 'EMBL-CDS')
        # mappingC.run()
        # # merge the two mapping results dataframes
        # mappingA.merge_mapping_dfs(mappingC.mapping_df)

        # # create targets and ncbi_columns and log not found targets
        # mappingA.finalize()
        # targets_and_ncbi_codes = mappingA.get_targets_and_ncbi_codes()  

        # # TODO: Temporary saving for faster debugging
        # parent_dir = os.path.dirname(os.getcwd())
        # file = os.path.join(parent_dir, 'debug_helper', 'targets_and_ncbi_codes.txt')
        # # # with open(file, 'w') as f:
        # # #     f.write(str(targets_and_ncbi_codes))
        # # read from file
        # with open(file, 'r') as f:
        #     targets_and_ncbi_codes = eval(f.read())

        # # C. Find assembly accession, download and parse assemblies
        # assemblies = Assemblies(config, targets_and_ncbi_codes)
        # assemblies.run()
        # flanking_genes = assemblies.get_flanking_genes()

        # # TODO: Temporary saving for faster debugging
        # parent_dir = os.path.dirname(os.getcwd())
        # file = os.path.join(parent_dir, 'debug_helper', 'flanking_genes.json')
        # # with open(file, 'w') as json_file:
        # #     json.dump(flanking_genes, json_file)
        # # read from file
        # with open(file, 'r') as json_file:
        #     flanking_genes = json.load(json_file)

        # # D. Add sequence information to flanking genes
        # sequences = Sequences(config, flanking_genes)
        # sequences.run()
        # genomic_context = sequences.get_genomic_context()

        # TODO: read faster debugging
        # read from file
        with open('genomic_context_information.json', 'r') as json_file:
            genomic_context = json.load(json_file)

        # E. Add protein families
        families = ProteinFamilies(config, genomic_context, out_label)



    # 4. 


    console.print_final()


    
if __name__ == '__main__':
    main()