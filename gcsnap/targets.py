import os
from typing import Union

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

class Target():
    def __init__(self, configuration: Configuration):
        self.arguments = configuration.arguments
        self.targets = configuration.targets
        self.console = RichConsole()

        # empty dictionary to store the targets
        self.targets_lists = {}

    def run(self) -> None:
        with self.console.status('Parsing targets'):
            self.parse_targets()

    def parse_targets(self, targets: Union[str, list] = None) -> dict:
        if targets is None:
            targets = self.targets

        for target in targets:
            # target is a file
            if os.path.isfile(target):
                if target.endswith('.clans'):
                    # TODO: the cluster is a dictionary, the rest are lists??
                    self.get_clusters_from_clans(target)
                else:
                    self.get_targets_from_file(target)
            else:
                # target is a list of ids
                # the name of this list (the label) by default is the 'default' 
                label = self.arguments['out_label']['value']
                if label not in self.targets_lists:
                    self.targets_lists[label] = []
                
                self.targets_lists[label].append(target)

    def get_targets_from_file(self, target_file: str) -> None:
        curr_label = target_file.split('/')[-1].split('.')[0]
        if curr_label not in self.targets_lists:
            self.targets_lists[curr_label] = []

        is_fasta = False
        with open(target_file, 'r') as f:
            content = f.readlines()

        for line in content:
            if line.startswith('>'):
                is_fasta = True
                curr_target = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()

            elif not is_fasta:
                curr_target = line.strip().split()[0]

            if curr_target not in self.targets_lists[curr_label]:
                self.targets_lists[curr_label].append(curr_target)        


    def get_clusters_from_clans(self, clans_file: str) -> dict:
        ncbis_ordered = self.get_ncbicodes_order_in_clans(clans_file)
        cluster_codes = self.arguments['cluster_patterns']['value']

        # TODO: This what bothers me, but I assume its the target list
        clusters = {}
        
        if cluster_codes != None:
            found_seqgroup = False
            with open(clans_file, 'r') as f:
                content = f.readlines()

            for line in content:
                for cluster_code in cluster_codes:
                    if '<seqgroup' in line:
                        found_seqgroup = True
                        found_allowed_cluster = False

                    elif found_seqgroup:
                        if 'name=' in line and cluster_code in line:
                            current_cluster = line.split('=')[-1].strip()
                            found_allowed_cluster = True
                        
                        elif 'numbers=' in line and found_allowed_cluster:
                            numbers = line.split('=')[-1].split(';')[:-1]
                            numbers = [int(i) for i in numbers]
                            numbers = [ncbis_ordered[i] for i in numbers]
                            self.targets_lists[current_cluster] = numbers

                            found_allowed_cluster = False

        else:
            label = clans_file.split('/')[-1].split('.')[0]
            self.targets_lists[label] = ncbis_ordered
    
    def get_ncbicodes_order_in_clans(self, clans_file: str) -> list:   
        ncbids_ordered = []

        found_seq = False
        with open(clans_file, 'r') as clans:
            for line in clans:
                if '<seq>' in line:
                    found_seq = True
                elif '</seq>' in line:
                    found_seq = False
                elif found_seq and line.startswith('>'):
                    line = line[1:]
                    ncbi_code = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()
                    ncbids_ordered.append(ncbi_code)
                    
        return ncbids_ordered
 