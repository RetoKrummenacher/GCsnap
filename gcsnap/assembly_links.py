import os
from datetime import datetime, timedelta
import time
# pip install pypdl
from pypdl import Pypdl

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

class AssemblyLinks:    
    """
    Methods and attributes to download and parse the assembly summaries from NCBI
    for RefSeq and Genbank to retreive the links to the assemblies.

    Attributes:
        console (RichConsole): The RichConsole object to print messages.
        cores (int): The number of CPU cores to use.
        assembly_dir (str): The path to store the assembly summaries.
        db_list (list): The list of databases to download the assembly summaries.
        links (dict): The final dictionary with the assembly code and the link to it.
    """

    def __init__(self, config: Configuration): 
        """
        Initialize the AssemblyLinks object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
        """           
        self.console = RichConsole()

        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value']
        self.age = config.arguments['assemblies_data_update_age']['value']

        parent_path = config.arguments['assemblies_data_folder']['value']
        
        if parent_path is None:
            # set path to store assembly summaries
            parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        self.assembly_dir = os.path.join(parent_path,'data','assembly_summaries')
        
        # database list
        self.db_list = ['genbank','refseq']
        
        # final dictionaire with {assembly code: link}
        self.links = {}
    
    def run(self) -> None:    
        """
        Run the download and parsing of the assembly summaries.
        """        
        self.create_folder()
        for db in self.db_list:
            self.check_and_download(db)        
            self.links.update(self.parse_summaries(db))      

    def get(self) -> dict[str, str]:
        """
        Getter for the links attribute.

        Returns:
            dict[str, str]: The dictionary with the assembly code and the link to it.
        """        
        return self.links              
        
    def create_folder(self) -> None:     
        """
        Create the folder to store the assembly summaries.
        """             
        if not os.path.exists(self.assembly_dir):
            os.makedirs(self.assembly_dir)
            
    def check_and_download(self, db) -> None:
        """
        Check if the assembly summary exists and download it if it does not.
        Using time check to download again if older than 14 days.

        Args:
            db (_type_): The database (genbank or refseq) to download the assembly summary.
        """        
        file = os.path.join(self.assembly_dir,'assembly_summary_{0}.txt'.format(db))
        days=self.age
        if os.path.exists(file):
            # time check, download again if older than days
            file_time = datetime.fromtimestamp(os.path.getmtime(file))
            if datetime.now() - file_time > timedelta(days=days):
                # print('Assembly summary is older than {} days, downloading again.'.format(days))
                self.console.print_info('Assembly summary {} is older than {} days, downloading again.'.format(db, days))
                self.download_summary(db)
            else:
                # print('Assembly summary is not older than {} days, not downloading.'.format(days))
                self.console.print_info('Assembly summary {} is not older than {} days, not downloading.'.format(db, days))
        else:
            # print('Assembly summary does not exist, downloading.')
            self.console.print_info('Assembly summary {} does not exist, downloading.'.format(db))
            self.download_summary(db)

    def download_summary(self, db: str) -> None:
        """
        Download the assembly summary from NCBI in parallel with pypdl Downloader.

        Args:
            db (str): The database (genbank or refseq) to download the assembly summary.
        """        
        url = 'https://ftp.ncbi.nlm.nih.gov/genomes/{0}/assembly_summary_{0}.txt'.format(db) 

        # Download multithreaded with pypdl Downloader
        dl = Pypdl()
        dl.start(url=url,
                file_path=self.assembly_dir, 
                segments=self.cores, 
                multisegment=True, 
                block=False, 
                display=False)
        
        # retreive total size of the file does not work, the file does not have this information
        # total_size = dl.size
        # print(f"Total size: {total_size}")

        # set it manually (lower than the actual size, but it is enough for the progress bar)
        total_size = 1 * 10**9 if db == self.db_list[0] else 0.1 * 10**9

        # Download multithreaded with pypdl Downloader
        with self.console.progress('Downloading {} assembly summary'.format(db), 
                                   total=total_size) as (progress, task_id):
             
            while not dl.completed:
                current_size = dl.current_size
                progress.update(task_id, completed=current_size)
                time.sleep(1)   

        dl.shutdown() 
        
    def parse_summaries(self, db: str) -> dict[str,str]:
        """
        Parse the assembly summary to extract the assembly code and the link to it.

        Args:
            db (str): The database (genbank or refseq) to parse the assembly summary.

        Returns:
            dict[str,str]: The dictionary with the assembly code and the link to it.
        """        
        links = {}
        file = os.path.join(self.assembly_dir,'assembly_summary_{0}.txt'.format(db))
        with open(file, 'r', encoding='utf-8') as f:
            content = f.readlines()

        with self.console.status(f'Parsing {db} assembly summary'):
            for line in content:
                if not line.startswith('#'):
                    data = line.strip().split('\t')
                    # the first column has the Assemlby id, column 19 the url to it
                    if data[19] != 'na':
                        links[data[0]] = data[19]
                    
        return links