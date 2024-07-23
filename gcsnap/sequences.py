
import json

from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration
from gcsnap.entrez_query import EntrezQuery
from gcsnap.genomic_context import GenomicContext

from gcsnap.utils import processpool_wrapper

class Sequences:
    """ 
    Methods and attributes to get the sequences for the flanking genes of the target genes.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cores (int): The number of CPU cores to use.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        sequences (dict): The dictionary with the sequences of the flanking genes.
        console (RichConsole): The RichConsole object to print messages.
    """

    def __init__(self, config: Configuration, gc: GenomicContext):
        """
        Initialize the Sequences object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
        """        
        self.config = config
        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value'] 

        # set arguments
        self.gc = gc

        self.console = RichConsole()

    def get_sequences(self) -> dict:
        """
        Getter for the sequences attribute.

        Returns:
            dict: The dictionary with flanking genes and their sequences.
        """        
        return self.genomic_context        

    def run(self) -> None:
        """
        Run the assignment of sequences to the flanking genes.
            - Find sequences for all flanking genes.
            - Add sequences, tax id and species name to flanking genes.
        Uses parallel processing with the processpool_wrapper from utils.py.
        """        
        # Find sequnces for all ncbi codes
        self.find_sequences(self.gc.get_all_ncbi_codes())

        # Prepare a list of tuples (target, dict_for_target)
        # here each process gets one target: {} combination
        # henve we have many different processes
        parallel_args = self.gc.get_syntenies_key_value_list()

        with self.console.status('Add sequences, tax id and species name to flanking genes'):
            dict_list = processpool_wrapper(self.cores, parallel_args, self.run_each)
            # combine results
            self.genomic_context = {k: v for d in dict_list for k, v in d.items()}

    def run_each(self, args: tuple[str,dict]) -> dict:
        """
        Run the assignment of sequences to the flanking genes for one target used
        in parallel processing.

        Args:
            args (tuple[str,dict]): The arguments for the sequence assignment.
                First element is the target gene.
                Second element is the dictionary with the flanking genes of the target gene.

        Returns:
            dict: The dictionary with the flanking genes and their sequences.
        """        
        target, content_dict = args
        # update flanking genes with sequence
        sequences = [self.get_sequence(ncbi_code) for ncbi_code in content_dict['flanking_genes']['ncbi_codes']]
        content_dict['flanking_genes']['sequences'] = sequences

        # add species and taxid for target_ncbi code (first one in the list)
        target_ncbi = content_dict['assembly_id'][0]
        # species in contained twice in the dict
        content_dict['flanking_genes']['species'] = self.get_species(target_ncbi)
        content_dict['species'] = content_dict['flanking_genes']['species']
        content_dict['flanking_genes']['taxID'] = self.get_taxid(target_ncbi)

        return {target: content_dict}

    def find_sequences(self, ncbi_codes: list) -> None:
        """
        Find sequences for all flanking genes using the EntrezQuery class.

        Args:
            ncbi_codes (list): The list of NCBI codes to find sequences for.
        """        
        # get the information for all ncbi codes
        # Noteworthy: While for accessions there were at most as many series as targets
        # here it is done for all flanking genes (1 + n_flanking_3 + n_flanking_5)
        entrez = EntrezQuery(self.config, ncbi_codes, db='protein', rettype='fasta', 
                             retmode='xml', logging=True)
        self.sequences = entrez.run()

    def get_sequence(self, ncbi_code: str) -> str:
        """
        Get the sequence for a ncbi code. 
        If the sequence is not found, a fake sequence is returned.

        Args:
            ncbi_code (str): The NCBI code.

        Returns:
            str: The sequence for the NCBI code.
        """        
        entry = self.sequences.get(ncbi_code, {})        
        return entry.get('seq', 'FAKESEQUENCEFAKESEQUENCEFAKESEQUENCEFAKESEQUENCE')

    def get_taxid(self, ncbi_code: str) -> str:
        """
        Get the tax id for a ncbi code.

        Args:
            ncbi_code (str): The NCBI code.

        Returns:
            str: The tax id for the NCBI code.
        """        
        entry = self.sequences.get(ncbi_code, {})        
        return entry.get('taxid','')
    
    def get_species(self, ncbi_code: str) -> str:
        """
        Get the species name for a ncbi code.

        Args:
            ncbi_code (str): The NCBI code.

        Returns:
            str: The species name for the NCBI code.
        """        
        entry = self.sequences.get(ncbi_code, {})        
        return entry.get('species','')    
    