import os
import subprocess
import shutil

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.sequence_mapping import SequenceMapping
from gcsnap.apis import UniProtAPI

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import split_list_chunks

import logging
logger = logging.getLogger(__name__) # inherits configuration from main logger

class TMsegments:
    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str):
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.annotate_TM = config.arguments['annotate_TM']['value']
        self.annotate_mode = config.arguments['annotate_TM_mode']['value']
        self.annotate_file = config.arguments['annotate_TM_file']['value']

        # set parameters
        self.gc = gc
        self.synthenies = gc.get_synthenies()
        self.out_label = out_label  
        self.out_dir = os.path.join(os.getcwd(), f'{out_label}_TM_annotations') 

        # outfile if not specified in config.yaml
        if self.annotate_file is None:
            # name of fasta file without extension
            self.annotation_out_file = os.path.join(self.out_dir, '{}_{}.out'.format(
                    os.path.basename(self.fasta_file)[:-6], self.annotate_mode)) 
        else:
            # the name of the input file
            self.annotation_out_file = os.path.join(self.out_dir, '{}'.format(
                    os.path.basename(self.annotation_out_file)))                

        self.console = RichConsole()

    def get_annotations(self) -> dict:
        return self.annotations

    def run(self) -> None:
        # abort if no annotation requested
        if not self.annotate_TM:
            self.console.print_skipped_step('TM annotation set to False. \
                                            Transmembrane segments and signal peptides will not be searched')
            self.annotations = {}
            return
        
        # create fasta file and get order if sequences in fasta file
        self.fasta_file = self.gc.write_to_fasta('flanking_sequences.fasta', 
                                                self.out_dir, exclude_pseudogenes = True)  
        self.ncbi_code_order = self.gc.get_fasta_order(False) 

        if self.annotate_file is None:
            self.single_annotation()
        else:
            # copy the specified file to the output directory
            try:
                shutil.copy(self.annotate_file, self.annotation_out_file)
            # catch when it is not there
            except FileNotFoundError:
                self.console.print_error('Your specifed file {} could not be found.'. \
                                         format(self.annotate_file)) 
                msg =  'Check file at speciefed path {} run GCsnap by specifying --annotation-TM-file.'. \
                        format(self.annotate_file)
                self.console.print_hint(msg)
                exit(1)                   
        
        # at this point, the file does exit
        with self.console.status('Parsing annotation file and updating flanking genes'):
            self.parse_annotation_file()
            self.update_flanking_genes()

    def single_annotation(self) -> None:
        # do single tm annotation depending on mode
        with self.console.status('Annotating TM segments with {}'.format(self.annotate_mode)):
            if self.annotate_mode == 'phobius' or self.annotate_mode == 'tmhmm':
                self.tool_annotation()
            else:
                # map all members to UniProtKB-AC
                mapping = SequenceMapping(self.config, self.ncbi_code_order, 
                                            to_type = 'UniProtKB-AC', msg = 'TM segments')
                mapping.run()
                mapping_dict = mapping.get_target_to_result_dict()
                mapping.log_failed()

                # create parallel arguments
                parallel_args = [(sub_list, mapping_dict) 
                                for sub_list in split_list_chunks(self.ncbi_code_order, self.cores)]

                result_lists = processpool_wrapper(self.cores, parallel_args, self.uniprot_annotation)
                # combine results
                result_list = [result for sub_list in result_lists for result in sub_list]   

                # create writing list
                write_list = ['{}\t{}'.format(ncbi_code, tm_annotation) 
                              for ncbi_code, tm_annotation in result_list]   
                self.write_to_file('\n'.join(write_list)) 

    def tool_annotation(self) -> None:
        if self.annotate_mode == 'phobius':
            tool = 'phobius.pl'
        elif self.annotate_mode == 'tmhmm':
            tool = 'tmhmm'

        # TODO: Why checking if file exists first?
        if not os.path.isfile(self.annotation_out_file):
            try:
                stdout, stderr = self.run_command(tool)
                if len(stderr) > 0 or len(stdout) == 0:
                    raise FileNotFoundError
                # if no error raised
                self.write_to_file(stdout)
            except FileNotFoundError:
                self.console.print_error('No {} installation was found.'.format(self.annotate_mode)) 
                msg =  'Run {} online and run GCsnap by specifying --annotation-TM-file.'. \
                        format(self.annotate_mode)
                self.console.print_hint(msg)
                self.extended_instructions()
                exit(1)          

    def run_command(self, tool: str) -> tuple:
        # returns stdout,stderr
        command = [tool, 
                    self.fasta_file,
                    '-short']
        # TODO: is this -short for tmhmm as well?
        
        result = subprocess.run(command, capture_output=True, text=True, shell=True)        
        return result.stdout, result.stderr   
    
    def uniprot_annotation(self, args = tuple) -> list[tuple]:
        ncbi_codes, mapping_dict = args

        result_list = []

        for ncbi_code in ncbi_codes:
            uniprot_code = mapping_dict.get(ncbi_code, 'nan')
            # if not found, it will be nan
            if uniprot_code != 'nan':
                result_list.append((ncbi_code, 'nan'))
                continue

            # get uniprot annotations, returns a dictionary
            curr_uniprot_annotations = UniProtAPI.get_uniprot_annotations(uniprot_code)

            if curr_uniprot_annotations != 'nan':
                tm_annotation = curr_uniprot_annotations['TM_topology']
                if tm_annotation == '':
                    tm_annotation = 'nan'
            else:
                tm_annotation = 'nan'

            result_list.append((ncbi_code, 'nan'))
        return result_list
    
    def parse_annotation_file(self) -> None:
        # read full file to list
        with open(self.annotation_out_file, 'r') as file:
            lines = file.readlines()

        # there would be no need to check for any other mode, as the parser
        # with choices took care of that
        if self.annotate_mode == 'phobius':
            self.protein_annotations = self.parse_phobius_output(lines)
        elif self.annotate_mode == 'tmhmm':
            self.protein_annotations = self.parse_tmhmm_output(lines)
        elif self.annotate_mode == 'uniprot':
            self.protein_annotations = self.parse_tm_uniprot_output(lines)
        else:
            self.console.print_warning('Annotation mode not recognized.')
            self.protein_annotations = {}

    def parse_phobius_output(self, lines: list) -> dict:
        annotations = {}
        for line in lines:
            if 'PREDICTION' not in line:
                ncbi_code = line.split('|')[0]
                tm_pred = int(line.split()[1])
                sp_pred = line.split()[2]

                if tm_pred > 0:
                    annotations[ncbi_code] = 'TM' # it has a transmembrane helix -> it is transmembrane
                elif sp_pred == 'Y':
                    annotations[ncbi_code] = 'SP' # it has a signal peptide
                else:
                    annotations[ncbi_code] = '' # does not have any special feature
        return annotations

    def parse_tmhmm_output(self, lines: list) -> dict:
        annotations = {}
        for line in lines:
            if 'PREDICTION' not in line:
                ncbi_code = line.split('|')[0]
                tm_pred = line.split('Topology=')[-1].strip()

                if tm_pred != 'o' and tm_pred != 'i':
                    annotations[ncbi_code] = 'TM' # it has a transmembrane helix -> it is transmembrane
                else:
                    annotations[ncbi_code] = '' # does not have any special feature (signal peptides are not identified)
        return annotations     

    def parse_tm_uniprot_output(self, lines: list) -> dict:
        annotations = {}
        for line in lines:
            ncbi_code = line.split()[0]
            tm_pred = line.split()[-1].strip()

            if tm_pred != 'nan':
                annotations[ncbi_code] = 'TM' 
            else:
                annotations[ncbi_code] = '' 
        return annotations       
    
    def update_flanking_genes(self) -> None:
        self.annotations = {}
        for curr_target in self.syntenies:
            flanking_genes = self.syntenies[curr_target]['flanking_genes']
            flanking_genes['TM_annotations'] = []

            for ncbi_code in flanking_genes['ncbi_codes']:
                if ncbi_code in self.protein_annotation:
                    flanking_genes['TM_annotations'].append(self.protein_annotations[ncbi_code])
                else:
                    flanking_genes['TM_annotations'].append('')

            self.annotations[curr_target]['flanking_genes'] = flanking_genes        

    def write_to_file(self, data: str) -> None:
        with open(self.annotation_out_file, 'w') as file:
            file.write(data)

    def extended_instructions(self) -> None:
        msg_output_file = 'Output file: copy-paste output text below the line to a .txt file'
        msg_save_folder = 'Provide path to file via --annotation-TM-file either in CLI or config.yaml'         
        if self.annotate_mode == 'phobius':
            msg_server = 'Phobius webserver: http://phobius.sbc.su.se/'
            msg_run_mode = 'Run mode: "short"'
        elif self.annotate_mode == 'tmhmm':
            msg_server = 'TMHMM webserver: https://services.healthtech.dtu.dk/service.php?TMHMM-2.0'
            msg_run_mode = 'Run mode: "One line per protein"'

        for msg in [msg_server, msg_run_mode, msg_output_file, msg_save_folder]:
            self.console.print_hint(msg)
