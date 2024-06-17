import os
import argparse

# pip install rich
from rich.console import Console
# pip install pyyaml
import yaml
# pip install ruamel.yaml
from ruamel.yaml.main import round_trip_dump as yaml_dump

from gcsnap.rich_console import RichConsole 


class CustomArgumentParser(argparse.ArgumentParser):

    def __init__(self, usage=None, epilog=None):
        super().__init__(usage=usage, epilog=epilog, add_help=False)

    def add_argument_from_config(self, config):
        for key, details in config.items():
            value_type = eval(details['type'])
            help_text = details['help']
            default_value = details['value']
            
            if value_type == bool:
                self.add_argument(f'--{key}', action='store_true' if default_value else 'store_false', help=help_text)
            else:
                self.add_argument(f'--{key}', type=value_type, default=default_value, help=help_text)


class Configuration:    
    def __init__(self):        
             
        self.path = os.path.dirname(__file__)
        self.indentation = 4 # space to use for indentation

        self.console = RichConsole()

        self.read_configuration_yaml()
        self.create_argument_parser()
        
    def read_configuration_yaml(self):       
        
        with open(os.path.join(self.path,'config.yaml'), 'r') as file:
            self.arguments = yaml.load(file, Loader=yaml.FullLoader)
            
    def write_configuration_yaml(self):
        if self.arguments['overwrite_config']['value']:
            out = self.arguments
            out.yaml_set_start_comment('Configuration file to set arguments for GCsnap.'
                                    'To change argument, change: value: entry.'
                                    'E.g. value: 1 to value: 2'
                                    '---------------------------------------\n')
            
            yaml_dump(os.path.join(self.path,out), indent=self.indentation, 
                    block_seq_indent=self.indentation)
        
    def parse_arguments(self):
        args = self.parser.parse_args()
        if hasattr(args, 'help') and args.help:
            self.console.print_help(self.parser)
            exit(0)

        # Ensure 'targets' argument is present
        if not hasattr(args, 'targets'):
            self.console.print_error('The following argument is required: --targets')
            exit(1)  # exit with error code 1

        # Update self.arguments dictionary with parsed values
        for arg in vars(args):
            if arg in self.arguments:
                self.arguments[arg]['value'] = getattr(args, arg)

        self.write_configuration_yaml()
        
    def create_argument_parser(self):
        usage = 'GCsnap --targets <targets> [Optional arguments]'
        epilog = 'Example: GCsnap --targets PHOL_ECOLI A0A0U4VKN7_9PSED A0'
        
        self.parser = CustomArgumentParser(usage=usage, epilog=epilog)

        # help argument, as we overwrote the default --help with rich printing
        self.parser.add_argument('-h', '--help', action='store_true', 
                                 help='Show this help message and exit')
        # consequences, we can't use required=True but check it in parse_arguments()

        # required arguments
        self.parser.add_argument('-t', '--targets', nargs='+', type=str, 
                                 help='List of input targets. Can be a list of fasta files,' + 
                                 'a list of text files encompassing a list of protein sequence identifiers,' + 
                                 'a list of protein sequence identifiers, or a mix of them')
        
        # Add arguments from YAML
        self.parser.add_argument_from_config(self.arguments)




 
# # GET INPUTS
# parser = argparse.ArgumentParser(prog = 'GCsnap v1.0.17', usage = 'GCsnap -targets <targets> -user_email <user_email> [options]', 
# 								 description = 'GCsnap is a python-based, local tool that generates interactive snapshots\nof conserved protein-coding genomic contexts.',
# 								 epilog = 'Example: GCsnap -targets PHOL_ECOLI A0A0U4VKN7_9PSED A0A0S1Y445_9BORD -user_email <user_email')
 
# requiredNamed = parser.add_argument_group('Required arguments')
# optionalNamed = parser.add_argument_group('Options with defaults')
 
# # required inputs
# requiredNamed.add_argument('-targets', dest='targets', nargs='+', required=True, help='List of input targets. Can be a list of fasta files, a list of text files encompassing a list of protein sequence identifiers, a list of protein sequence identifiers, or a mix of them')
# # optional inputs
# optionalNamed.add_argument('-user_email', dest='user_email', type=str, default = None, help='Email address of the user. May be required to access NCBI databases and is not used for anything else (default: None)')
# optionalNamed.add_argument('-ncbi_api_key', dest='ncbi_api_key', default = None,type=str, help='The key for NCBI API, which allows for up to 10 queries per second to NCBI databases. Shall be obtained after obtaining an NCBI account (default: None)')
# optionalNamed.add_argument('-get_taxonomy', dest='get_taxonomy', default = 'True',type=str, help='Boolean statement to get and map taxonomy information (default: True)')
# optionalNamed.add_argument('-cpu', dest='n_cpu', default = 1,type=int, help='Number of cpus to use (default: 1)')
# optionalNamed.add_argument('-n_flanking', dest='n_flanking', default = 4,type=int, help='Number of flanking sequences (to each side) to take (default: 4)')
# optionalNamed.add_argument('-n_flanking5', dest='n_flanking5', default = 4,type=int, help="Number of flanking sequences to take on the 5' (default: 4)")
# optionalNamed.add_argument('-n_flanking3', dest='n_flanking3', default = 4,type=int, help="Number of flanking sequences to take on the 3' (default: 4)")
# optionalNamed.add_argument('-exclude_partial', dest='exclude_partial', default = False,type=bool, help='Boolean statement to exclude partial operon/genomic_context blocks (default: False)\nIf turned off, partial cases will still be ignored to get the most common genomic features')
# optionalNamed.add_argument('-out_label', dest='out_label', default = 'default',type=str, help='The label to append to the out folder (default: "default"). Important when the input list corresponds to raw sequence identifiers.')
# optionalNamed.add_argument('-out_label_suffix', dest='out_label_suffix', default = '',type=str, help='A suffix to add to the out_label (default: "").')
# optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp',type=str, help='The temporary folder (default: /tmp). May be changed so that intermediary files (e.g., assembly files) are saved somewhere else.')
# optionalNamed.add_argument('-collect_only', dest='collect_only', default = False,type=bool, help='Boolean statement to make GCsnap collect genomic contexts only, without comparing them (default: False).')
# # operon clustering
# optionalNamed.add_argument('-n_max_operons', dest='n_max', default = 30,type=int, help='Maximum number of top most populated operon/genomic_context block types (default: 30)')
# optionalNamed.add_argument('-operon_cluster_advanced', dest='operon_cluster_advanced', default = False,type=bool, help='Boolean statement to use the operon clustering advanced mode (using PacMAP) (default: False)')
# optionalNamed.add_argument('-max_family_freq', dest='max_family_freq', default = 20,type=int, help='Maximum frequency of a family in the set of genomic cotexts found to be considered for advanced operon clustering (default: 20)')
# optionalNamed.add_argument('-min_family_freq', dest='min_family_freq', default = 2,type=int, help='Minimum frequency of a family in the set of genomic cotexts found to be considered for advanced operon clustering (default: 20)')
# #protein family identification
# optionalNamed.add_argument('-n_iterations', dest='num_iterations', default = 1,type=int, help='Number of iterations for all-against-all searches (default: 1). Required to define protein families.')
# optionalNamed.add_argument('-evalue', dest='max_evalue', default = 1e-3,type=float, help='Max e-value at which two sequences are considered to be homologous (default: 1e-3). Required to define protein families.')
# optionalNamed.add_argument('-coverage', dest='min_coverage', default = 70,type=float, help='Minimum coverage of target and subject a match needs to be so that two sequences are considered to be homologous (default: 70). Required to define protein families.')
# optionalNamed.add_argument('-base', dest='default_base', default = 10,type=int, help='Artificial distance value for two sequences that do not match with an E-value better than -evalue (default: 10).')
# optionalNamed.add_argument('-all-against-all_method', dest='clustering_method', default = 'psiblast',type=str, choices=['mmseqs', 'psiblast'], help='Method for clustering (default: psiblast)')
# optionalNamed.add_argument('-psiblast_location', dest='blast', default = 'psiblast',type=str, help='Location of psiBLAST (if not in path) (default: psiblast)')
# optionalNamed.add_argument('-mmseqs_location', dest='mmseqs', default = 'mmseqs',type=str, help='Location of MMseqs (if not in path) (default: mmseqs)')
# # figure making optional inputs
# optionalNamed.add_argument('-genomic_context_cmap', dest='genomic_context_cmap', default = 'Spectral',type=str, help='Color map (as of matplotlib) to assign colors to and plot the syntenic blocks (default: Spectral)')
# optionalNamed.add_argument('-out_format', dest='out_format', default = 'png',type=str, help='Output format of the core figures (default: png)')
# optionalNamed.add_argument('-print_color_summary', dest='print_color_summary', default = False, type=bool, help='Boolean statement to print the RGBA codes of the colors defined (default: False)')
# # annotation optional inputs
# optionalNamed.add_argument('-get_pdb', dest='get_pdbs', default = 'True', type=str, help='Boolean statement to get PDB information for representatives of the families found (default: True)\nTurn off to make it faster.')
# optionalNamed.add_argument('-get_functional_annotations', dest='get_functional_annotations', default = 'True' ,type=str, help='Boolean statement to find functional annotations for representatives of the families found (default: True)\nTurn off to make it faster.')
# optionalNamed.add_argument('-annotate_TM', dest='annotate_TM', default = False, type=bool, help='Boolean statement to find sequence features in the flanking genes (default: False)')
# optionalNamed.add_argument('-annotation_TM_mode', dest='annotation_TM_mode', default = 'uniprot', type=str, choices=['phobius', 'tmhmm', 'uniprot'], help='Method to use to find transmembrane segments (default: uniprot)')
# optionalNamed.add_argument('-annotation_TM_file', dest='annotation_TM_file', default = None, type=str, help='File with pre-computed transmembrane features. Only use when the targets correspond to a single project (no multiple fasta or text files) (default: None)')
# # interactive optional inputs
# optionalNamed.add_argument('-interactive', dest='interactive', default = 'True',type=str, help='Boolean statement to make the interactive html output (default: True). WARNING: It requires the Bokeh python package. It will check if it is installed')
# optionalNamed.add_argument('-gc_legend_mode', dest='gc_legend_mode', default = 'species',type=str, choices=['species', 'ncbi_code'], help='Mode of the genomic context legend (default: species)')
# optionalNamed.add_argument('-min_coocc', dest='min_coocc', default = 0.30,type=float,help='Minimum maximum co-occurrence of two genes to be connected in the graphs (default: 0.30)')
# optionalNamed.add_argument('-min_freq_accross_contexts', dest='min_family_freq_accross_contexts', default = 30,type=float,help='Minimum frequency of a family in a conserved genomic context type to be considered as a member (default: 30)')
# optionalNamed.add_argument('-sort_mode', dest='sort_mode', default = 'taxonomy',type=str, choices=['taxonomy', 'as_input', 'tree', 'operon'], help='Mode to sort the genomic contexts (default: taxonomy)')
# optionalNamed.add_argument('-in_tree', dest='in_tree', default = None, type=str, help='Input phylogenetic tree. Only use when the targets correspond to a single project (no multiple fasta or text files) (default: None)')
# optionalNamed.add_argument('-in_tree_format', dest='in_tree_format', default = "newick", type=str, help='Format of the input phylogenetic tree (default: newick)')
# # clans map optional inputs
# optionalNamed.add_argument('-clans_patterns', dest='clans_patterns', default = None,type=str, nargs='+', help='Patterns to identify the clusters to analyse. They will be used to select the individual clusters in the clans map to analyse (default: None).')
# optionalNamed.add_argument('-clans_file', dest='clans_file', default = None,type=str, help='Clans file. Used only for advanced interactive output representation (default: None or input clans file if the input is a clans file and -operon_cluster_advanced is set to True).')
 
