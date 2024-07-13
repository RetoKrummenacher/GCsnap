import os
import subprocess
import numpy as np
# pip install scipy
from scipy.cluster import hierarchy
from scipy.spatial import distance

from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.rich_console import RichConsole

class MMseqsCluster:
    def __init__(self, config: Configuration, gc: GenomicContext, out_dir: str):
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.max_evalue = config.arguments['max_evalue']['value']
        self.min_coverage = config.arguments['min_coverage']['value']
        self.num_iterations = config.arguments['num_iterations']['value']
        self.mmseqs_executable = r'{}'.format(config.arguments['mmseqs_executable_path']['value'])
        self.default_base = config.arguments['default_base']['value']

        # set arguments
        self.gc = gc
        self.out_dir = out_dir
  
        self.sensitivity = 7.5

        self.console = RichConsole()

    def run(self) -> None:
        with self.console.status('Prepare data for MMseqs'):
            self.fasta_file = self.gc.write_to_fasta('flanking_sequences.fasta', 
                                            self.out_dir, exclude_pseudogenes = False)  
            self.cluster_order = self.gc.get_fasta_order(exclude_pseudogenes = False) 
        with self.console.status('Running MMseqs'):            
            self.run_mmseqs()
        with self.console.status('Extracting distance matrix'):
            self.extract_distance_matrix()
        with self.console.status('Find clusters'):
            self.find_clusters()
            self.mask_singleton_clusters()        

    def get_distance_matrix(self) -> np.array:
        return self.distance_matrix

    def get_clusters_list(self) -> list[int]:
        return self.cluster_list      

    def get_cluster_order(self) -> list[str]:
        return self.cluster_order       

    def run_mmseqs(self) -> None:
            self.mmseqs_results = os.path.join(self.out_dir, '{}_{}.mmseqs'.format(
                os.path.basename(self.fasta_file)[:-6], self.max_evalue))
            
            # TODO: Why checking if file exists first?
            if not os.path.isfile(self.mmseqs_results):
                try:
                    _, stderr = self.mmseqs_command('mmseqs')
                    if len(stderr) > 0:
                        raise FileNotFoundError
                except FileNotFoundError:
                    try:
                        _, stderr = self.mmseqs_command(self.mmseqs_executable)
                    except:
                        self.console.print_error('No MMseqs installation was found.') 
                        self.console.print_hint('Please install MMseqs or add the path to the executable to config.yaml.')
                        exit(1)                    

    def mmseqs_command(self, mmseqs :str) -> tuple:
        # returns stdout,stderr
        command = [mmseqs, 
                'easy-search', 
                self.fasta_file, 
                self.fasta_file, 
                self.mmseqs_results, 
                self.out_dir, 
                '-e', str(self.max_evalue), 
                '-s', str(self.sensitivity),
                '-c', str(self.min_coverage),
                '--num-iterations', str(self.num_iterations),
                '--threads', str(self.cores),
                '--format-output',
                'query,target,evalue']
        
        result = subprocess.run(command, capture_output=True, text=True, shell=True)        
        return result.stdout, result.stderr       
    
    def extract_distance_matrix(self) -> None:
        # crate base distance matrix
        distance_matrix = [[self.default_base if i!=j else 0 for i in self.cluster_order] 
                        for j in self.cluster_order]
        queries_labels = {query: i for i, query in enumerate(self.cluster_order)}

        # read mmseqs results
        with open(self.mmseqs_results, 'r') as f:
            mmseqs_records = f.readlines()

        for hsp in mmseqs_records:
            hsp = hsp.split()
            if len(hsp) > 0:
                query = hsp[0].split('|')[0]
                query_index = queries_labels[query]
                target = hsp[1].split('|')[0]
                if target != query:
                    target_index = queries_labels[target]
                    distance_matrix[query_index][target_index] = 0
                    distance_matrix[target_index][query_index] = 0

        self.distance_matrix = np.array(distance_matrix)      

    def find_clusters(self, t: int = 0) -> None:
        distance_matrix = distance.squareform(self.distance_matrix)
        linkage = hierarchy.linkage(distance_matrix, method = 'single')
        clusters = hierarchy.fcluster(linkage, t, criterion = 'distance')
        self.cluster_list = [int(i) for i in clusters]

    def mask_singleton_clusters(self, mask: int = 0) -> None:
        self.cluster_list = [mask if list(self.cluster_list).count(value) == 1 
                             else value for value in self.cluster_list]      