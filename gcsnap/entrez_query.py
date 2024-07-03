import requests
import time
import datetime
import random
import numpy as np
import xml.etree.ElementTree as ET

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import WarningToLog

import logging
logger = logging.getLogger(__name__) # inherits configuration from main logger

class EntrezQuery:
    def __init__(self, config: Configuration, search_codes: list, db: str = 'protein', 
                 rettype: str = 'ipg', retmode: str = 'xml', logging: bool = True):
        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value']
        self.api_key = config.arguments['ncbi_api_key']['value']

        self.console = RichConsole()

        # set arguments
        self.db = db
        self.rettype = rettype
        self.retmode = retmode
        self.search_codes = search_codes
        self.logging = logging

        # set correct function depending on rettype
        if rettype == 'ipg':
            self.msg = 'Find identical protein groups'
            self.function = self.run_accession
            self.chunk_size = 50
        elif rettype == 'fasta':
            self.msg = 'Download sequences'
            self.function = self.run_sequence
            self.chunk_size = 100
        else:
            self.msg = 'Download taxonomy'
            self.function = self.run_taxonomy

    def run(self) -> dict:
        with self.console.status(self.msg):
            chunk_list = self.create_chunks(self.chunk_size)
            # returns a result tuple with (found_dict, None) or ({}, (WarningToLog, chunk))
            results_list = processpool_wrapper(self.cores, chunk_list, self.function)
            # combine all dictionaries
            found = {k: v for result_tuple in results_list for k, v 
                     in result_tuple[0].items() if result_tuple[1] is None}
            # get those with failed entrez requests
            errors = [(e, chunk) for result_tuple in results_list if result_tuple[1] 
                      is not None for e, chunk in result_tuple[1]]
 
        # log errors if desired
        failed_request_ids = []
        if errors and self.logging:
            failed_request_ids = self.log_request_error(errors)

        # log not found if desired
        not_found = [id for id in self.search_codes if id not in found.keys() 
                     and id not in failed_request_ids]
        if not_found and self.logging:
            self.log_not_found(not_found)

        return found
        
    def run_accession(self, chunk: list) -> tuple[dict,tuple]:
        try:
            root = self.query_chunk(chunk)
            return (self.extract_accessions(chunk, root) , None)       
        except WarningToLog as e:
            return ({}, (e, chunk))      
    
    def run_sequence(self, chunk: list) -> tuple[dict,tuple]:
        try:
            root = self.query_chunk(chunk)
            return (self.extract_sequences(root), None)
        except WarningToLog as e:
            return ({}, (e, chunk))

    def create_chunks(self, chunk_size: int = 50) -> list:
        return [self.search_codes[i:i + chunk_size] for i in range(0, len(self.search_codes), chunk_size)]

    def get_ncbi_base_url(self) -> str:
        return 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db='+self.db

    def query_chunk(self, chunk: list) -> ET.Element:
        ids = '&id=' + ','.join(chunk)
        rt = '&rettype=' + self.rettype
        rm = '&retmode=' + self.retmode

        # add api key if provided
        if self.api_key is not None:
            api_key = '&api_key=' + self.api_key
        else:
            api_key = ''

        # set timeout depending on chunk size
        #timeout = len(chunk_list) * 0.2
        timeout = None

        url = self.get_ncbi_base_url() + ids + rt + rm + api_key

        # Entrez has an API limit of 3 request per second without an API key
        # and 10 request per second with an API key. If no API key is specified
        # we need to handle '429 Too Many Requests error'
        initial_wait_time = 1  # Initial wait time in seconds
        max_wait_time = 10  # Maximum wait time in seconds
        factor = 2  # Exponential backoff factor

        while True:
            response = requests.get(url, timeout=timeout)
            if response.status_code == 200:
                xml_string = response.text
                return ET.fromstring(xml_string) 
            elif response.status_code == 429:
                # extract the time with wait suggestion from the response
                retry_after = response.headers.get('Retry-After')
                if retry_after:
                    try:
                        retry_after = int(retry_after)
                    except ValueError:
                        # If Retry-After is not a number, it might be a date
                        retry_after_date = requests.utils.parse_date(retry_after)
                        if retry_after_date:
                            wait_time = (retry_after_date - datetime.datetime.now()).total_seconds()
                    wait_time = int(retry_after)
                else:
                    # Decorrelated Jitter
                    wait_time = min(random.uniform(initial_wait_time, wait_time
                                                * factor), max_wait_time)
                time.sleep(wait_time)
            else:
                raise WarningToLog('Entrez request failed with status code {}'.format(response.status_code))
  
    def extract_accessions(self, chunk_list: list, root: ET.Element) -> dict:
        reports = root.findall('IPGReport') # one report for each id
        # ProteinList: full identical protein group for each id
        ipg_list = [(reports[i].find('ProteinList'), chunk_list[i]) for i in range(len(reports))
                    if reports[i].find('ProteinList') is not None]

        return_dict = {}
        for ipg, ncbi_code in ipg_list: # iterate though all ipgs
            protein_list = ipg.findall('Protein') # all proteins from the ipg group
            
            # all individual CDS of the full sequnce of the protein
            cds_lists = [(protein.find('CDSList').findall('CDS'), protein.get('accver')) 
                         for protein in protein_list if protein.find('CDSList') is not None]
            # flatten and exclude those without assembly accession
            cds_list = [(item,sublist[1]) for sublist in cds_lists for item in sublist[0]
                        if item.get('assembly') is not None]

            assembly_acc = None
            for cds, accver in cds_list:    
                if ncbi_code == accver: 
                    # the assembly accession for the desired ncbi_code
                    assembly_acc_tmp = cds.get('assembly')
                    if assembly_acc_tmp.startswith('GCF'):
                        # if the assembly accession is a GCF (RefSeq), use it directly
                        return_dict[ncbi_code] = assembly_acc_tmp
                        break
                    # if the assembly accession is a GCA (GenBank), continue search
                    assembly_acc = assembly_acc_tmp

            # if still None, don't add and catch later by logger
            if assembly_acc is not None:
                return_dict[ncbi_code] = assembly_acc
        return return_dict
    
    def extract_sequences(self, root: ET.Element) -> dict:
        # create dictionary for all found sequences
        sequences = {report.find('TSeq_accver').text : {'seq': report.find('TSeq_sequence').text,
                                'taxid': report.find('TSeq_taxid').text,
                                'species': report.find('TSeq_orgname').text} 
                    for report in root.findall('TSeq')}        
        return sequences
    
    def log_not_found(self, not_found: list) -> None:
        message = '{} entrez queries returned no results for db "{}" with rettype "{}"'.format(len(not_found), self.db, self.rettype)
        self.console.print_warning(message)
        for id in not_found:
            logger.warning('No entrez entry found for {} with db "{}" and rettype "{}"'.format(id, self.db, self.rettype)) 

    def log_request_error(self, errors: list) -> list:
        message = '{} entrez requests failed for db "{}" with rettype "{}"'.format(len(errors), self.db, self.rettype)
        self.console.print_warning(message)
        # all filed ncbi_codes
        all_failed = []
        for error, chunk in errors:
            logger.warning(error)
            for each in chunk:
                logger.warning('Request failed for {}'.format(each))
                all_failed.append(each)
        return all_failed