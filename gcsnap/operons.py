import numpy as np
from typing import Union
# pip install pacmap
import pacmap
# pip install scikit-learn
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture
# pip install scipy
from scipy.spatial import distance
from scipy.cluster import hierarchy

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext

import logging
logger = logging.getLogger(__name__) # inherits configuration from main logger

class Operons:
    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str):
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.advanced = config.arguments['operon_cluster_advanced']['value']
        self.max_freq = config.arguments['max_family_freq']['value']
        self.min_freq = config.arguments['min_family_freq']['value']

        # get families and syntenies from the synties object
        self.families = gc.get_families()
        self.syntenies = gc.get_syntenies()
        self.out_label = out_label

        self.console = RichConsole() 

        # data structures created within the class
        # set default values if nothing has to be done in cases only 
        # one target is present in syntenies
        self.sorted_targets = list(self.syntenies.keys())
        self.operon_clusters = [1] * len(self.sorted_targets)           # list

    def get_operons(self) -> None:
        return self.syntenies        

    def run(self) -> None:
        # if syntenies has only one target, the default values are kept
        if len(self.syntenies) > 1:
            if self.advanced:
                with self.console.status('Running advanced operon clustering with PaCMAP'):
                    self.run_advanced()
            else:
                with self.console.status('Running standard operon clustering'):
                    self.run_standard()

        # update syntenies with the new operon clusters
        self.update_syntenies()       

    def run_advanced(self, args: tuple) -> None:
        # get the clusters by excluding the most common families
        res = self.find_operon_clusters_with_PaCMAP(clean = True, coordinates_only = False) 
        self.clean_coordinates, self.operon_clusters, self.sorted_targets = res
        # and now all PacMap coordinates by using all families. 
        # This will be later used for sorting the dendogram             
        res = self.find_operon_clusters_with_PaCMAP(clean = False, coordinates_only = True) 
        self.all_coordinates, self.sorted_targets = res

    def find_operon_clusters_with_PaCMAP(self, clean: bool = True, 
                                         coordinates_only: bool = True,
                                         max_eps_cost: float = 1.3, 
                                         eps_step: float = 0.02
                                         ) -> Union[tuple[np.ndarray,list],
                                                    tuple[np.ndarray,np.ndarray,list]]:
        
        presence_matrix, sorted_targets, sorted_families = self.get_family_presence_matrix(clean)
        paCMAP_embedding = pacmap.PaCMAP(n_components = 2)
        paCMAP_coordinat = paCMAP_embedding.fit_transform(presence_matrix)          

        if coordinates_only:
            return paCMAP_coordinat, sorted_targets
        else:
            # embed into n-D pacMAP space
            n_dims = len(sorted_families)
            if n_dims < 2:
                n_dims = 2
            
            paCMAP_embedding = pacmap.PaCMAP(n_components = n_dims)
            paCMAP_coordinat = paCMAP_embedding.fit_transform(presence_matrix)

		# find clusters in the paCMAP space
		# do this by selecting the best eps based on the number of clusters it creates compared to the number of operons 
		# that are not assigned a clusters (i.e., a given maximum cost)        
        eps = self.calculate_start_eps(paCMAP_coordinat)    

        n_clusters = [0]
        n_singletons = [0]
        cost = 0
        logger.info('Refining EPS for DBSCAN clustering based on Gaussian Mixture Model')
        while cost <= max_eps_cost and eps > 0:
            eps = eps - eps_step

            model = DBSCAN(eps = eps)
            model.fit(paCMAP_coordinat)
            clusters = model.fit_predict(paCMAP_coordinat)

            n = len(set(clusters))
            delta_n_clusters = n - n_clusters[-1]
            delta_singletons = list(clusters).count(-1) - n_singletons[-1]
            if delta_n_clusters > 0:
                cost = delta_singletons/delta_n_clusters
            else:
                cost = 0

        logger.info('Calculated EPS: {}'.format(eps))
        logger.info('Cost: {}'.format(cost))  
        logger.info('Number of clusters: {}'.format(n))              

        return paCMAP_coordinat, clusters, sorted_targets                

    def get_family_presence_matrix(self, clean: bool = True) -> tuple[np.ndarray, list, list]:
        sorted_targets = sorted(list(self.syntenies.keys()))
        sorted_families   = [i for i in sorted(list(self.families.keys())) if (i >= 0 and i < 10000)]

        # select only the protein families that are not very frequenct but also not very rare
        if clean and len(sorted_families) > 10:
            families_frequency = [len(self.families[family]['members']) for family in sorted_families]
            families_frequency = [i*100/len(self.syntenies) for i in families_frequency]
            sorted_families	= [family for i, family in enumerate(sorted_families) 
                               if families_frequency[i] <= self.max_freq 
                               and families_frequency[i] >= self.min_freq]

        presence_matrix = [[0 for i in sorted_families] for i in sorted_targets]
        for i, target_i in enumerate(sorted_targets):
            operon_i = self.syntenies[target_i]['flanking_genes']['families']
            for family in operon_i:
                if family in sorted_families:
                    presence_matrix[i][sorted_families.index(family)] += 1

        return np.array(presence_matrix), sorted_targets, sorted_families   

    def calculate_start_eps(coordinates: list) -> float:
        distances = []
        for i, vector_i in enumerate(coordinates):
            for j, vector_j in enumerate(coordinates):
                if j>i:
                    dist = np.linalg.norm(np.array(vector_i) - np.array(vector_j))
                    distances.append(dist)

        distances = np.array(distances)
        mixture = GaussianMixture(n_components=2).fit(distances.reshape(-1,1))
        means = mixture.means_.flatten()
        sds = np.sqrt(mixture.covariances_).flatten()

        mean = min(means)
        sd = sds[list(means).index(mean)]

        eps = mean - sd        
        return round(eps, 2)    

    def run_standard(self) -> None:
        distance_matrix, self.sorted_targets = self.compute_operon_distance_matrix()
        self.operon_cluster = self.find_clusters_in_distance_matrix(distance_matrix)

    def compute_operon_distance_matrix(self) -> tuple[np.ndarray,list]:
        distance_matrix = [[0 for target in self.syntenies] for target in self.syntenies]
        sorted_targets = sorted(list(self.syntenies.keys()))

        for i, operon_i in enumerate(sorted_targets):
            for j, operon_j in enumerate(sorted_targets):
                if i < j:
                    vector_i = self.syntenies[operon_i]['flanking_genes']['families']
                    reference_family_i = self.syntenies[operon_i]['target_family']
                    
                    vector_j = self.syntenies[operon_j]['flanking_genes']['families']
                    reference_family_j = self.syntenies[operon_j]['target_family']

                    if len(vector_i) >= len(vector_j):
                        reference_vector = vector_i
                        curr_vector = vector_j
                        index_reference_family_reference_vector = (
                            list(reference_vector).index(reference_family_i))
                        index_reference_family_curr_vector = (
                            list(curr_vector).index(reference_family_j))
                    else:
                        reference_vector = vector_j
                        curr_vector = vector_i
                        index_reference_family_reference_vector = (
                            list(reference_vector).index(reference_family_j))
                        index_reference_family_curr_vector = (
                            list(curr_vector).index(reference_family_i))

                    number_fillings = 0
                    if len(curr_vector) < len(reference_vector):
                        if (index_reference_family_curr_vector < 
                            index_reference_family_reference_vector):
                            for a in range(index_reference_family_reference_vector - 
                                           index_reference_family_curr_vector):
                                curr_vector = np.insert(curr_vector, 0, -1)
                                number_fillings += 1
                        
                        if len(curr_vector) < len(reference_vector):
                            for b in range(len(reference_vector) - len(curr_vector)):
                                curr_vector = np.append(curr_vector, -1)
                                number_fillings += 1

                    families_present = set(np.concatenate([reference_vector,curr_vector]))
                    overlapping_families = set(reference_vector).intersection(set(curr_vector))
                    dist = 1-len(overlapping_families)/len(families_present)

                    distance_matrix[i][j] = dist
                    distance_matrix[j][i] = dist

        return np.array(distance_matrix), sorted_targets       

    def find_clusters_in_distance_matrix(self, distance_matrix: np.ndarray, t: int = 0.2) -> list: 
        distance_matrix = distance.squareform(distance_matrix)
        linkage = hierarchy.linkage(distance_matrix, method = 'single')
        clusters = hierarchy.fcluster(linkage, t, criterion = 'distance')
        return [int(i) for i in clusters]
    
    def update_syntenies(self) -> None:
        for i, target in enumerate(self.sorted_targets):
            self.syntenies[target]['operon_type'] = int(self.operon_clusters[i])

            if self.advanced:
                self.syntenies[target]['operon_filtered_PaCMAP'] = (
                    list([float(a) for a in self.clean_coordinates[i]]))
                self.syntenies[target]['operon_PaCMAP'] = (
                    list([float(a) for a in self.all_coordinates[i]]))