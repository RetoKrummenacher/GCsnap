import os
import gzip # to work with .gz
import urllib.request
import time

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.assembly_links import AssemblyLinks
from gcsnap.entrez_query import EntrezQuery
from gcsnap.genomic_context import GenomicContext

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import WarningToLog

import logging
logger = logging.getLogger(__name__) # inherits configuration from main logger

class Assemblies:    
    def __init__(self, config: Configuration, mappings: list[tuple[str,str]]):                     
        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value'] 
        self.n_flanking5 = config.arguments['n_flanking5']['value']  
        self.n_flanking3 = config.arguments['n_flanking3']['value']
        self.exclude_partial = config.arguments['exclude_partial']['value'] 
        self.config = config

        self.console = RichConsole()

        # set path to store assembly summaries
        parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.assembly_dir = os.path.join(parent_path,'data','assemblies')
        self.create_folder()

        # input list with [(target, ncbi_code)]
        self.targets_and_ncbi_codes = mappings

    def get_flanking_genes(self) -> dict:
        return self.flanking_genes
    
    def run(self) -> None:
        # download and parse the assembly summary files
        self.load_summaries()
        # find the assembly accessions
        ncbi_codes = [target[1] for target in self.targets_and_ncbi_codes]
        self.find_accessions(ncbi_codes)

        with self.console.status('Download assemblies and extract flanking genes'):
            dict_list = processpool_wrapper(self.cores, self.targets_and_ncbi_codes, self.run_each)
            # combine results
            self.flanking_genes = {k: v for d in dict_list for k, v in d.items() 
                                   if v.get('flanking_genes') is not None}
            not_found = {k: v for d in dict_list for k, v in d.items() 
                         if v.get('flanking_genes') is None}
        if not_found:
            self.log_not_found(not_found)

    def run_each(self, args: tuple[str,str]) -> dict[str, dict]:
        # TODO: What to do if no accession/url is found?
        # adapt the two getter functions to return.
        # right now, it raises an exception which returns None for flanking genes
        # those are exluded from the flanking_genes dict and logged at the end.
        target, ncbi_code = args
        try:
            accession = self.get_assembly_accession(ncbi_code)
            assembly_url = self.get_assembly_url(accession)
            assembly_file, lines = self.download_and_read_gz_file(assembly_url)
            flanking_genes = self.parse_assembly(ncbi_code, lines)
            return {target: {'flanking_genes': flanking_genes,
                             'assembly_id':  [ncbi_code, accession, assembly_url]}}   
        except WarningToLog as e:
            # return None for flanking genes and message, logged later
            return {target: {'flanking_genes': None,
                             'msg': str(e)}}
        
    def create_folder(self) -> None:          
        if not os.path.exists(self.assembly_dir):
            os.makedirs(self.assembly_dir)

    def load_summaries(self) -> None:
        # get file with assembly links
        assembly_links = AssemblyLinks(self.config)
        assembly_links.run()
        self.links = assembly_links.get()   

    def find_accessions(self, ncbi_codes: list) -> None:
        # get file with assembly links, no logging as its done after run()
        entrez = EntrezQuery(self.config, ncbi_codes, db='protein', rettype='ipg', 
                             retmode='xml', logging=False)
        self.accessions = entrez.run()

    def get_assembly_accession(self, ncbi_code: str) -> str:
        accession = self.accessions.get(ncbi_code, 'unk')
        if accession == 'unk':
            raise WarningToLog('No assembly accession found for {}'.format(ncbi_code))
        return accession

    def get_assembly_url(self, assembly_accession: str) -> str:
        url = self.links.get(assembly_accession, 'unk')   
        if url == 'unk':
            raise WarningToLog('No url found for accession {}'.format(assembly_accession))
        return url      
        
    def download_and_read_gz_file(self, url: str, retries: int = 3) -> tuple:
        try:
            full_path = self.download_gz_file(url)
            content = self.read_gz_file(full_path)
        except (urllib.error.URLError, EOFError) as e:        
            if retries > 0:     
                time.sleep(1)
                full_path, content = self.download_and_read_gz_file(url, retries - 1)
            else:
                raise WarningToLog('Download failed for {} with error {}'.format(url, e))

        return full_path, content

    def download_gz_file(self, url: str) -> str:
        assembly_label = url.split('/')[-1]  # e.g., GCF_000260135.1_ASM26013v1
        assembly_file_gz = '{}_genomic.gff.gz'.format(assembly_label)
        full_path = os.path.join(self.assembly_dir, assembly_file_gz)

        # TODO: Skip this for experiments and developpment
        # Check if file already exists
        # if os.path.exists(full_path):
        #     return full_path

        # Download the .gz file and save it without uncompressing
        with urllib.request.urlopen(url + '/' + assembly_file_gz) as response:
            with open(full_path, 'wb') as out_file:
                out_file.write(response.read())

        return full_path

    def read_gz_file(self, file_path: str) -> list:
        with gzip.open(file_path, 'rt', encoding='utf-8') as file:
            content = file.read()
        return content.splitlines()                
                
    def delete_assemblies(self) -> None:
        for filename in os.listdir(self.assembly_dir):
            file_path = os.path.join(self.assembly_dir, filename)
            os.remove(file_path)

    def parse_assembly(self, ncbi_code: str, lines: list) -> dict:
        genomic_context_block = self.extract_genomic_context_block(ncbi_code, lines)
        return self.parse_genomic_context_block(ncbi_code, genomic_context_block)

    def extract_genomic_context_block(self, ncbi_code: str, lines: list) -> list:
         # line numbers where different scaffolds (sequence regions) start
        scaffold_positions = [0] + [index for index, val in enumerate(lines) if 
                            val.startswith('##sequence-region')] + [len(lines)]   
             
        # line number of target
        target_position = [index for index, val in enumerate(lines) if ncbi_code in val]
        
        if not target_position:
            raise WarningToLog('{} not found in'.format(ncbi_code))

        # select scaffold region with target
        region_start_end = [(scaffold_positions[i], scaffold_positions[i + 1]) for i in 
                range(len(scaffold_positions) - 1) if scaffold_positions[i] 
                <= target_position[0] <= scaffold_positions[i + 1]]
        
        # extract the scaffold and lines containing coding sequence (CDS) information
        scaffold = [line for line in lines[region_start_end[0][0]:region_start_end[0][1]] if 
                    (len(line.split('\t')) >= 3 and line.split('\t')[2] == 'CDS')]
        
        # target position in scaffold
        index_of_target = [index for index, val in enumerate(scaffold) 
                           if ncbi_code in val][0]
        
        # need to know direction to define what flanking genes to extract
        direction_of_target = scaffold[index_of_target].split('\t')[6]
        
        # define neighbor indizes depending on direction
        if direction_of_target == '+':
            # upper slice index can be out of bounds without error, lower must be non negative
            start = 0 if (index_of_target - self.n_flanking5) < 0 else index_of_target - self.n_flanking5
            end = index_of_target + self.n_flanking3 + 1
            # extract the genomic context 
            genomic_context_block = scaffold[start:end]
        else:
            # upper slice index can be out of bounds without error, lower must be non negative
            start = 0 if (index_of_target - self.n_flanking3) < 0 else index_of_target - self.n_flanking3
            end = index_of_target + self.n_flanking5 + 1
            # extract the genomic context 
            genomic_context_block = scaffold[start:end]
            # reverse if direction of target is '-'
            genomic_context_block = genomic_context_block[::-1]
        
        # exclude partials if desired
        if self.exclude_partial and len(genomic_context_block) < (self.n_flanking5 + self.n_flanking3 + 1):
            raise WarningToLog('Partial genomic block for {} excluded!'.format(ncbi_code))
        
        return genomic_context_block

    def parse_genomic_context_block(self, target_ncbi_code: str, genomic_context_block: list) -> dict:
        # result dictionary
        flanking_genes = GenomicContext.get_empty_flanking_genes()

        # parse the genomic context
        for line in genomic_context_block:   
            # one line of the genomic context block:
                # JQ926483.1	Genbank	CDS	1	668	.	+	0	
                # ID=cds-AFI40896.1;Parent=gene-VP1;Dbxref=NCBI_GP:AFI40896.1;
                # Name=AFI40896.1;end_range=668,.;gbkey=CDS;gene=VP1;partial=true;
                # product=RNA-dependent RNA polymerase;protein_id=AFI40896.1;start_range=.,1

            line_data = line.split('\t')            
            start = int(line_data[3])
            end = int(line_data[4])
            direction = line_data[6]
            
            if 'cds-' in line:
                ncbi_code = line_data[8].split('ID=cds-')[1].split(';')[0]
            elif 'Name=' in line:
                ncbi_code = line_data[8].split('Name=')[1].split(';')[0]
            else:
                ncbi_code = 'unk'

            if 'pseudo=' not in line and 'product=' in line and 'fragment' not in line:
                prot_name = line_data[8].split('product=')[1].split(';')[0]
            else:
                prot_name = 'pseudogene'

            # it means that this is some kind of fragmented gene (has introns?...) 
            # and so we have to collect the largest interval encompassing it
            if ncbi_code in flanking_genes['ncbi_codes'] and flanking_genes['ncbi_codes'][-1] == ncbi_code: 
                if start < flanking_genes['starts'][-1]:
                    flanking_genes['starts'][-1] = start
                if end > flanking_genes['ends'][-1]:
                    flanking_genes['ends'][-1] = end
            else:
                if '|' in ncbi_code:
                    ncbi_code = ncbi_code.replace('|','_')

                flanking_genes['ncbi_codes'].append(ncbi_code)
                flanking_genes['names'].append(prot_name)
                flanking_genes['starts'].append(start)
                flanking_genes['ends'].append(end)
                flanking_genes['directions'].append(direction)

        # index of target in flanking genes
        index_of_target = flanking_genes['ncbi_codes'].index(target_ncbi_code)

        # direction of target
        direction_of_target = flanking_genes['directions'][index_of_target]
        
        # add relative starts and ends depending on direction
        if direction_of_target == '+':
            for key in ['starts','ends']:
                lst = [e - flanking_genes['starts'][index_of_target] + 1 for e in flanking_genes[key]]
                flanking_genes['relative_{}'.format(key)] = lst  
        else:
            # order is reversed, old ends determin the starts and vice-versa
            # old end of target is the new base point
            lst = [e - flanking_genes['ends'][index_of_target] + 1 for e in flanking_genes['ends']]
            flanking_genes['relative_starts'] = lst  
            lst = [e - flanking_genes['ends'][index_of_target] + 1 for e in flanking_genes['starts']]
            flanking_genes['relative_ends'] = lst  
            
            # turn + and - directions
            flanking_genes['directions'] = ['+' if d == '-' else '-' for d in flanking_genes['directions']]
            
        return flanking_genes       

    def log_not_found(self, not_found: dict) -> None:
        message = 'No flanking genes found for {} target sequences.'.format(len(not_found))
        self.console.print_warning(message)
        for k,v in not_found.items():
            logger.warning('For target {}: {}'.format(k,v.get('msg'))) 
