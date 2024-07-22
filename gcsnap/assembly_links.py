import os
from datetime import datetime, timedelta
import time
# pip install pypdl
from pypdl import Pypdl

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

class AssemblyLinks:    
    def __init__(self, config: Configuration):    
        self.console = RichConsole()

        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value']
        
        # set path to store assembly summaries
        parent_path = os.path.dirname(os.getcwd())
        self.assembly_dir = os.path.join(parent_path,'data','assembly_summaries')
        
        # database list
        self.db_list = ['genbank','refseq']
        
        # final dictionaire with {assembly code: link}
        self.links = {}
    
    def run(self) -> None:    
        self.create_folder()
        for db in self.db_list:
            self.check_and_download(db)        
            self.links.update(self.parse_summaries(db))      

    def get(self) -> dict[str, str]:
        return self.links              
        
    def create_folder(self) -> None:          
        if not os.path.exists(self.assembly_dir):
            os.makedirs(self.assembly_dir)
            
    def check_and_download(self, db) -> None:
        file = os.path.join(self.assembly_dir,'assembly_summary_{0}.txt'.format(db))
        days=14
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
        links = {}
        file = os.path.join(self.assembly_dir,'assembly_summary_{0}.txt'.format(db))
        with open(file, 'r', encoding='utf-8') as f:
            content = f.readlines()

        with self.console.status(f'Parsing {db} assembly summary'):
            for line in content:
                if not line.startswith('#'):
                    data = line.strip().split('\t')
                    # the first column has the Assemlby id, column 19 the url to it
                    links[data[0]] = data[19]
                    
        return links