import os
import shutil
# TODO: For faster debugging
import pickle

from gcsnap.rich_console import RichConsole 
from gcsnap.configuration import Configuration 
from gcsnap.timing import Timing
from gcsnap.targets import Target 
from gcsnap.sequence_mapping import SequenceMapping
from gcsnap.assemblies import Assemblies
from gcsnap.genomic_context import GenomicContext
from gcsnap.sequences import Sequences
from gcsnap.families import Families
from gcsnap.families_functions_structures import FamiliesFunctionsStructures
from gcsnap.operons import Operons
from gcsnap.taxonomy import Taxonomy
from gcsnap.tm_segments import TMsegments
from gcsnap.figures import Figures

def main():

    starting_directory = os.getcwd()

    console = RichConsole()
    console.print_title()

    # 1. Parse configuration and arguments
    config = Configuration()
    config.parse_arguments()

    # 2. start timing
    timing = Timing()
    t_all = timing.timer('All steps 0-10')

    t_parse = timing.timer('Step 0: Parse Targets')
    # 2. parse targets
    targets = Target(config)
    targets.run()
    t_parse.stop()

    # 3. Iterate over each target list
    for out_label in targets.targets_lists:
        # all execution conditions depending on arguments from CLI or config.yaml
        # are handled in the classes themselves

        # A. Prework 
        working_dir = os.path.join(starting_directory, out_label)
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
        os.chdir(working_dir)
        targets_list = targets.targets_lists[out_label]
        if len(targets_list) < 2:
            console.print_warning('GCsnap was asked to analyze only one target')
            console.print_skipped_step('Skipping target {}'.format(out_label))
            continue
        else:
            console.print_working_on('Analyzing target list: {}'.format(out_label))
        # write configuration to log file
        config.write_configuration_yaml_log('{}_input_arguments.log'.format(out_label))
        # Datastructure to store all information
        gc = GenomicContext(config, out_label)
        # add targets to genomic context
        gc.targets = targets_list

        t_collect = timing.timer('Step 1: Collecting the genomic contexts')
        # B. Map sequences to UniProtKB-AC and NCBI EMBL-CDS
        # a). Map all targets to UniProtKB-AC
        mappingA = SequenceMapping(config, targets_list, 'UniProtKB-AC')
        mappingA.run()
        # b) Map all to RefSeq
        mappingB = SequenceMapping(config, mappingA.get_codes(), 'RefSeq')
        mappingB.run()
        # merge them to A (only if A is not nan)
        mappingA.merge_mapping_dfs(mappingB.mapping_df)
        # c). Map all targets to NCBI EMBL-CDS
        mappingC = SequenceMapping(config, mappingA.get_codes(), 'EMBL-CDS')
        mappingC.run()
        # merge the two mapping results dataframes
        mappingA.merge_mapping_dfs(mappingC.mapping_df)
        # create targets and ncbi_columns and log not found targets
        mappingA.finalize()
        targets_and_ncbi_codes = mappingA.get_targets_and_ncbi_codes()  

        # C. Find assembly accession, download and parse assemblies
        assemblies = Assemblies(config, targets_and_ncbi_codes)
        assemblies.run()
        gc.update_syntenies(assemblies.get_flanking_genes())

        # D. Add sequence information to flanking genes
        sequences = Sequences(config, gc)
        sequences.run()
        gc.update_syntenies(sequences.get_sequences())
        gc.write_syntenies_to_json('genomic_context_information.json')

        t_collect.stop()

        if not config.arguments['collect_only']['value']:

            t_family = timing.timer('Step 2: Finding protein families')
            # Ea) Add protein families
            families = Families(config, gc, out_label)
            families.run()
            gc.update_syntenies(families.get_families())
            gc.create_and_write_families_summary()

            t_family.stop()

            t_annotate_families = timing.timer('Step 3: Annotating functions and structures')
            # Eb). Add functions and structures to families
            # execution conditions handeled in the class
            ffs = FamiliesFunctionsStructures(config, gc)
            ffs.run()
            gc.update_families(ffs.get_annotations_and_structures())
            gc.write_families_to_json('protein_families_summary.json')

            t_annotate_families.stop()
        
            t_operons = timing.timer('Step 4-5: Finding operon/genomic_context')
            # F. Find and add operons
            operons = Operons(config, gc, out_label)
            operons.run()
            gc.update_syntenies(operons.get_operons())
            gc.create_and_write_operon_types_summary()
            gc.find_most_populated_operon_types()   

            t_operons.stop()

            t_taxonomy = timing.timer('Step 6: Mapping taxonomy')
            # G. Get taxonomy information
            taxonomy = Taxonomy(config, gc)
            taxonomy.run()
            gc.update_taxonomy(taxonomy.get_taxonomy())
            gc.write_taxonomy_to_json('taxonomy.json')     

            t_taxonomy.stop()
  
            t_tm = timing.timer('Step 7: Finding ALL proteins with transmembrane segments')
            # H. Annotate TM 
            tm = TMsegments(config, gc, out_label)
            tm.run()
            gc.update_syntenies(tm.get_annotations())

            t_tm.stop()

            # # TODO: For debugging
            # with open('gc.pkl', 'wb') as file:
            #     pickle.dump(gc, file)  

            # # TODO: For debugging        
            # with open('gc.pkl', 'rb') as file:
            #     gc = pickle.load(file) 
             
            t_figures = timing.timer('Step 8-9: Producing figures')
            # I. Produce genomic context figures
            figures = Figures(config, gc, out_label, starting_directory)      
            figures.run()

            t_figures.stop()

            t_output = timing.timer('Step 10: Write output')

            # G. Write output to summary file
            # TODO: Still a bug line 372 in gc:
            #     line_to_write += '\t' + '\t'.join(tax_search_dict.get(target)) + '\n'
            #               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            # TypeError: can only join an iterable
            # the return fomr dict.get() is None. check that dict
            #gc.write_summary_table('{}_summary_table.tab'.format(out_label))

            gc.write_families_to_json('protein_families_summary.json')
        
        else:
            console.print_skipped_step('GCsnap was asked to collect genomic context only. Will not proceed further.')
            t_output = timing.timer('Step 10: Write output')

        # J. Wrap up
        gc.write_syntenies_to_json('all_syntenies.json')
        t_output.stop()

        if config.arguments['overwrite_config']['value']:
            config.write_configuration_yaml()
        # copy log file to working direcotry
        shutil.copy(os.path.join(starting_directory,'gcsnap.log'), os.getcwd())

    t_all.stop()
    timing.to_csv('timing.csv')
    console.print_final()
    
if __name__ == '__main__':
    main()