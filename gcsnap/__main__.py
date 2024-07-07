import os
# TODO: For faster debugging
import json

from gcsnap.rich_console import RichConsole 
from gcsnap.configuration import Configuration 
from gcsnap.targets import Target 
from gcsnap.sequence_mapping_online import SequenceMappingOnline
from gcsnap.assemblies import Assemblies
from gcsnap.genomic_context import GenomicContext
from gcsnap.sequences import Sequences
from gcsnap.families import Families
from gcsnap.families_functions_structures import FamiliesFunctionsStructures
from gcsnap.operons import Operons
from gcsnap.taxonomy import Taxonomy
from gcsnap.tm_segments import TMsegments

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
        # all execution conditions depending on arguments from CLI or config.yaml
        # are handled in the classes themselves

        # TODO: For debugging
        gc = GenomicContext(config, out_label)

        # A. Create working directory
        working_dir = os.path.join(starting_directory, out_label)
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

        # # # TODO: Temporary saving for faster debugging
        # with open('targets_and_ncbi_codes.txt', 'w') as f:
        #     f.write(str(targets_and_ncbi_codes))
        # # # read from file
        # with open('targets_and_ncbi_codes.txt', 'r') as f:
        #     targets_and_ncbi_codes = eval(f.read())


        # # C. Find assembly accession, download and parse assemblies
        # assemblies = Assemblies(config, targets_and_ncbi_codes)
        # syntenies = assemblies.run()
        # # Datastructure to store all information
        # gc = GenomicContext(config, out_label)
        # gc.update_syntenies(assemblies.get_flanking_genes())

        # # TODO: Temporary saving for faster debugging
        # gc.write_syntenies_to_json('flanking_genes.json')
        # gc.read_syntenies_from_json('flanking_genes.json')


        # # D. Add sequence information to flanking genes
        # sequences = Sequences(config, gc)
        # sequences.run()
        # gc.update_syntenies(sequences.get_sequences())
        # gc.write_syntenies_to_json('genomic_context_information.json')

        # TODO: read faster debugging
        # read from file
        gc.read_syntenies_from_json('genomic_context_information.json')


        # Ea) Add protein families
        families = Families(config, gc, out_label)
        families.run()
        gc.update_syntenies(families.get_families())
        gc.create_and_write_families_summary()

        # Eb). Add functions and structures to families
        # execution conditions handeled in the class
        ffs = FamiliesFunctionsStructures(config, gc)
        ffs.run()
        gc.update_families(ffs.get_annotations_and_structures())
        gc.write_families_to_json('protein_families_summary.json')

        # TODO: For debugging
        gc.read_families_from_json('protein_families_summary.json')

    # 4. 


    console.print_final()


    
if __name__ == '__main__':
    main()