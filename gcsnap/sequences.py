import os
import json

from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.db_handler_sequences import SequenceDBHandler
from gcsnap.parallel_tools import ParallelTools
from gcsnap.utils import split_dict_chunks


class Sequences:
    """ 
    Methods and attributes to get the sequences for the flanking genes of the target genes.
    They are stored in a structure like:
        data-path as defined in config.yaml or via CLI
        ├── genbank
        │   └── data
        │       └── GCA_000001405.15_genomic.gff.gz
        ├── refseq
        │   └── data
        │       └── GCF_000001405.38_genomic.gff.gz
        ├── db
        │   └── assemblies.db
        │   └── mappings.db
        │   └── sequences.db
        │   └── rankedlineage.dmp

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        n_nodes (int): The number of nodes to use for parallel processing.
        n_cpu (int): The number of CPUs to use for parallel processing.
        n_worker_chunks (int): The number of chunks to split the workload into.
        database_path (str): The path to the database.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        sequences_dict (dict): The dictionary with the sequences of the flanking genes.
        console (RichConsole): The RichConsole object to print messages.
        sequences (dict): The results dictionary with the flanking genes and their sequences.
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
        self.n_nodes = config.arguments['n_nodes']['value']
        self.n_cpu = config.arguments['n_cpu_per_node']['value']
        self.n_worker_chunks = config.arguments['n_worker_chunks']['value']
        self.database_path = os.path.join(config.arguments['data_path']['value'],'db') 
        self.gc = gc

        # in how many sizes to split the workload
        # 1 reserved for the main process
        self.work_split = ((self.n_nodes * self.n_cpu) - 1) * 2

        self.console = RichConsole()

    def get_sequences(self) -> dict:
        """
        Getter for the sequences attribute.

        Returns:
            dict: The dictionary with flanking genes and their sequences.
        """        
        return self.sequences        

    def run(self) -> None:
        """
        Run the assignment of sequences to the flanking genes.
            - Find sequences for all flanking genes.
            - Add sequences, tax id and species name to flanking genes.
        Uses parallel processing with the parallel_wrapper from ParallelTools.
        """        
        # Get syntenies
        syntneies = self.gc.get_syntenies()

        # extract target and neeeded flanking genes ncbi codes
        part_syntenies = {k: v['flanking_genes']['ncbi_codes'] for k, v in syntneies.items()}
        # split into chunks
        parallel_args = split_dict_chunks(part_syntenies , self.n_worker_chunks)

        with self.console.status('Add sequences, tax id and species name to flanking genes'):
            dict_list = ParallelTools.parallel_wrapper(parallel_args, self.run_each)
            # combine results
            self.sequences = {k: v for d in dict_list for k, v in d.items()}

    def run_each(self, args: list[dict]) -> dict:
        """
        Run the assignment of sequences to the flanking genes for one target used
        in parallel processing.

        Args:
            args (tuple[dict]): The arguments for the sequence assignment.
                The dictionary with the flanking gene information

        Returns:
            dict: The dictionary with the flanking genes and their sequences.
        """        
        part_syntenies = args

        # get all ncbi codes
        ncbi_codes = [ncbi_code for content_dict in part_syntenies.values() 
                      for ncbi_code in content_dict['flanking_genes']['ncbi_codes']]

        # get from database
        self.find_sequences(ncbi_codes)

        # adapt syntenies
        for target, content_dict in part_syntenies.items():
            # update flanking genes with sequence
            # update flanking genes with sequence
            sequence_list = [self.get_sequence(ncbi_code) for ncbi_code in content_dict['flanking_genes']['ncbi_codes']]
            content_dict['flanking_genes']['sequences'] = sequence_list
            part_syntenies |= {target: content_dict}

        return part_syntenies

    def find_sequences(self, ncbi_codes: list) -> None:
        """
        Find sequences for all flanking genes using the EntrezQuery class.

        Args:
            ncbi_codes (list): The list of NCBI codes to find sequences for.
        """        
        # select from database
        sequences_db = SequenceDBHandler(os.path.join(self.database_path))
        self.sequences_dict = sequences_db.select_as_dict(ncbi_codes)

    def get_sequence(self, ncbi_code: str) -> str:
        """
        Get the sequence for a ncbi code. 
        If the sequence is not found, a fake sequence is returned.

        Args:
            ncbi_code (str): The NCBI code.

        Returns:
            str: The sequence for the NCBI code.
        """        
        entry = self.sequences_dict.get(ncbi_code, {})        
        return entry.get('seq', 'FAKESEQUENCEFAKESEQUENCEFAKESEQUENCEFAKESEQUENCE')

    