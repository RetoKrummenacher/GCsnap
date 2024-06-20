import os

from gcsnap.rich_console import RichConsole 
from gcsnap.configuration import Configuration 
from gcsnap.targets import Target 
from gcsnap.sequence_mapping_online import SequenceMappingOnline

def main():

    starting_directory = os.getcwd()

    console = RichConsole()
    console.print_title()

    # 1. Parse configuration and arguments
    config = Configuration()
    config.parse_arguments()

    # 2. parse targets
    with console.status('Parsing targets'):
        targets = Target(config)
        targets.run()

    # 3. Iterate over each target list
    for out_label in targets.targets_lists:

        # A. Create working directory
        # TODO: Change if out_label_suffix is still used
        working_dir = '{}/{}'.format(starting_directory, out_label)
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
        os.chdir(working_dir)

        targets_list = targets.targets_lists[out_label]

        with console.status('Mapping sequences'):
            # B. Map all targets to UniProtKB-AC
            mappingB = SequenceMappingOnline(config, targets_list, 'UniProtKB-AC')
            mappingB.run()

            # C. Map all targets to NCBI EMBL-CDS
            mappingC = SequenceMappingOnline(config, mappingB.result_list, 'EMBL-CDS')
            mappingC.run()
            # merge the two mapping results
            mappingB.merge_mapping_dfs(mappingC.mapping_df)

        # D. Determin and extract assemblies
        


        # with console.progress('Count',10000000) as (progress, task_id):
        #     i = 0
        #     while i < 10000000:
        #         i += 1
        #         progress.update(task_id, advance=1)


    # 4. 


    console.print_final()


    
if __name__ == '__main__':
    main()