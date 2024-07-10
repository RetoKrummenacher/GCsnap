import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import statistics
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial import distance
# pip install pacmap
import pacmap

from collections import Counter

import networkx as nx
from Bio import Phylo

from bokeh.plotting import figure, output_file, gridplot, save
from bokeh.colors import RGB
from bokeh.models import HoverTool, TapTool, LassoSelectTool, Range1d, LinearAxis, WheelZoomTool, Circle, MultiLine, Panel, Tabs
from bokeh.models import ColumnDataSource, DataTable, DateFormatter, TableColumn, Legend, HTMLTemplateFormatter
from bokeh.models.callbacks import OpenURL
from bokeh.models.graphs import from_networkx, NodesAndLinkedEdges
from bokeh.models.widgets import Div
from bokeh.layouts import row, column

from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.figure_genomic_context import GenomicContextFigure
from gcsnap.figure_interactive import InteractiveFigure
from gcsnap.figure_interactive_advanced import AdvancedInteractiveFigure

class Figures:
    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str, starting_directory: str):
        self.config = config
        self.opernon_cluster_advanced = config.arguments['operon_cluster_advanced']['value']
        self.cmap = config.arguments['genomic_context_cmap']['value']

        # set parameters
        self.gc = gc
        self.operons = gc.get_selected_operons()
        self.most_populated_operon = gc.get_most_populated_operon()
        self.synthenise = gc.get_synthenise()
        self.families_summary = gc.get_families()
        self.out_label = out_label
        Figures.starting_directory = starting_directory

        # define the reference family as the one 
        # of the target in the most populated operon type
        self.reference_family = self.synthenise[self.operons[self.most_populated_operon]
                                                ['target_members'][0]]['target_family']
        

    def run(self) -> None:
        # 1. Create genomic context figures
        family_colors = self.define_family_colors(mode = 'matplotlib', cmap = self.cmap)
        gcf = GenomicContextFigure(self.config, self.gc, self.out_label, 
                                   self.reference_family, family_colors)  
        gcf.run()

        # 2. Interactive HTML figures
        family_colors = self.define_family_colors(mode = 'bokeh', cmap = self.cmap)          
        if self.opernon_cluster_advanced:
            aigcf = AdvancedInteractiveFigure(self.config, self.gc, self.out_label,
                                            self.reference_family, family_colors)
            aigcf.run()
        else:         
            igcf = InteractiveFigure(self.config, self.gc, self.out_label,
                                    self.reference_family, family_colors)
            igcf.run()

    def define_family_colors(self, mode: str,  cmap: str):
        families = list(self.families_summary.keys())
        colors = {}

        cmap = matplotlib.cm.get_cmap(cmap)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(families))

        colours = [cmap(norm(i)) for i in range(len(families))]
        #random.shuffle(colours)

        for i, label in enumerate(sorted(families)):
            if label not in colors:
                colors[label] = {}

            if label == self.reference_family:			   # the reference gene
                colors[label]['Color (RGBA)'] = 'grey'
                colors[label]['Color (tuplet)'] = 'grey'
                colors[label]['Line color'] = 'black'
                colors[label]['Line style'] = '-'
            elif label == 0:			# a gene without any other homolog in the figure
                colors[label]['Color (RGBA)'] = 'lightgrey'
                colors[label]['Color (tuplet)'] = 'lightgrey'
                colors[label]['Line color'] = 'lightgrey'
                colors[label]['Line style'] = '-'
            elif type(label) == int and label == 10000:  # a pseudogene
                colors[label]['Color (RGBA)'] = 'white'
                colors[label]['Color (tuplet)'] = 'white'
                colors[label]['Line color'] = 'grey'
                colors[label]['Line style'] = ':'
            else:						# real genes that occur many times in the figure
                if mode == 'matplotlib':
                    colors[label]['Color (RGBA)'] = [int(255*j) for j in colours[i]]
                    colors[label]['Color (tuplet)'] = colours[i]
                    colors[label]['Line color'] = 'black'
                    colors[label]['Line style'] = '-'
                elif mode == 'bokeh':
                    colors[label]['Color (RGBA)'] = RGB(int(255*list(colours[i])[0]), 
                                                        int(255*list(colours[i])[1]), 
                                                        int(255*list(colours[i])[2]))
                    colors[label]['Color (tuplet)'] = colours[i]
                    colors[label]['Line color'] = 'black'
                    colors[label]['Line style'] = '-'

        return colors   

    @staticmethod
    def set_class_attributes(cls, **kwargs):
        """Static method to set attributes on the given class."""
        for key, value in kwargs.items():
            setattr(cls, key, value)  
    
    @staticmethod
    def make_dendogram_figure(**kwargs) -> tuple[figure, dict]:
        # set attributes from kwargs, it updates any keywords that are passed to method
        Figures.set_class_attributes(Figures, **kwargs)

        
        dendogram, leaf_labels = Figures.compute_dendogram(**kwargs)
        den_tooltips, den_data = Figures.create_dendogram_features(dendogram = dendogram, 
                                                                leaf_labels = leaf_labels, **kwargs)

        y_range = [0, max(den_data['y'])+min(den_data['y'])+(den_data['y'][1]-den_data['y'][0])]

        if len(leaf_labels) < 5:
            height = int(len(leaf_labels)*Figures.height_factor)*3
        elif len(leaf_labels) < 10:
            height = int(len(leaf_labels)*Figures.height_factor)*2
        else:
            height = int(len(leaf_labels)*Figures.height_factor)

        if Figures.in_tree != None:
            title = 'Input phylogeny/hierarchy'
        elif Figures.sort_mode == 'taxonomy':
            title = 'Taxonomy hierarchy'
        elif Figures.sort_mode == 'operons':
            title = 'Genomic contexts similarity hierarchy'
        else:
            title = 'Genomic contexts clusters hierarchy'

        if Figures.show_leafs:
            den = figure(title = title, height=height, width=500, 
                         x_range=[-max(max(dendogram['dcoord']))-
                                  (max(max(dendogram['dcoord']))-min(min(dendogram['dcoord'])))*0.2, 40], 
                         y_range = y_range, toolbar_location="left")
        else:
            den = figure(title = title, height=height, width=250, 
                         x_range=[-max(max(dendogram['dcoord']))-
                                  (max(max(dendogram['dcoord']))-min(min(dendogram['dcoord'])))*0.2, 0.2], 
                         y_range = y_range, toolbar_location="left")

        for i, d in zip(dendogram['icoord'], dendogram['dcoord']):
            d = list(map(lambda x: -x, d))
            den.line(x=d, y=i, line_color='black')

        if Figures.show_leafs:
            leafs = den.text(x='x', y='y', text = 'leaf_label', text_baseline='middle', 
                             text_font_size='8pt', source = den_data)
        else:
            leafs = den.circle(x=0.1, y='y', line_color = 'black', fill_color = 'color', 
                               source = den_data, size = (den_data['y'][1] - den_data['y'][0])*0.7)
        den.add_tools(HoverTool(tooltips=den_tooltips, renderers=[leafs]))

        den.title.align = "right"
        den.axis.major_tick_line_color = None
        den.axis.minor_tick_line_color = None
        den.axis.major_label_text_color = None
        den.axis.major_label_text_font_size = '0pt'
        den.axis.axis_line_color = None
        den.grid.grid_line_color = None
        den.outline_line_color = None

        return den, den_data        
    
    @staticmethod
    def compute_dendogram(**kwargs) -> tuple[dict, list]:
        Figures.set_class_attributes(Figures, **kwargs)
        
        if distance_matrix is None:
            if Figures.in_tree != None:
                distance_matrix, labels = Figures.get_phylogeny_distance_matrix(print_tree = True, **kwargs)
            elif Figures.sort_mode == 'operon':
                distance_matrix, labels = Figures.get_operons_distance_matrix(**kwargs)
            elif Figures.sort_mode == 'operon clusters':
                distance_matrix, labels = Figures.get_avgoperons_distance_matrix(**kwargs)		  
            else:
                distance_matrix, labels = Figures.get_taxonomy_distance_matrix(**kwargs)

        distance_matrix = distance.squareform(distance_matrix)

        Z = hierarchy.linkage(distance_matrix, method = 'single')
        results = hierarchy.dendrogram(Z, no_plot=True, count_sort = 'descending')

        if Figures.sort_mode == 'taxonomy' or Figures.in_tree != None or 'operon' in Figures.sort_mode:
            labels_indeces = list(map(int, results['ivl']))
            labels = [labels[i] for i in labels_indeces]

        return results, labels    
    
    @staticmethod
    def get_phylogeny_distance_matrix(**kwargs) -> tuple[np.ndarray, list]:
        Figures.set_class_attributes(Figures, **kwargs)

        if Figures.starting_directory not in Figures.in_tree:
            Figures.in_tree = os.path.join(os.getcwd(),'{}'.format(Figures.in_tree))

        tree = Phylo.read(Figures.in_tree, Figures.in_tree_format)  

        if Figures.print_tree:
            Phylo.draw_ascii(tree)

        tree.ladderize()
        tree_depths = tree.depths()	
        if not max(tree_depths.values()): 
            tree_depths = tree.depths(unit_branch_lengths=True) 

        max_depth = max([tree_depths[clade] for clade in tree_depths if clade.name != None])

        # align everything to the right by extending the branches
        all_leafs = []
        for clade in tree.get_terminals():
            curr_depth = tree_depths[clade]
            clade.branch_length += max_depth - curr_depth
            all_leafs.append(clade.name)

        # remove branches that are not in our dataset
        for leaf in all_leafs:
            if leaf not in Figures.input_targets:
                tree.prune(leaf)

        labels = [leaf.name for leaf in tree.get_terminals()]

        distance_matrix = [[0 for leaf in labels] for leaf in labels]
        for i, leaf_i in enumerate(labels):
            for j, leaf_j in enumerate(labels):
                distance_matrix[i][j] = tree.distance(leaf_i, leaf_j)

        distance_matrix = np.array(distance_matrix)

        Phylo.write(tree, 'out_tree.nwk', 'newick')

        return distance_matrix, labels
    
    @staticmethod
    def get_operons_distance_matrix(**kwargs) -> tuple[np.ndarray, list]:
        # set attributes from kwargs, it updates any keywords that are passed to method
        Figures.set_class_attributes(Figures, **kwargs)       

        input_targets = kwargs.get('input_targets', None)

        labels = []
        if 'operon_filtered_PaCMAP' in Figures.operons[list(Figures.operons.keys())[0]]:
            paCMAP_coordinat = []
            for operon_type in Figures.operons:
                for i, target in enumerate(Figures.operons[operon_type]['target_members']):
                    if input_targets is None or target in input_targets:
                        paCMAP_coordinat.append(Figures.operons[operon_type]['operon_filtered_PaCMAP'][i])
                        labels.append(target)

            distance_matrix = [[0 for i in paCMAP_coordinat] for j in paCMAP_coordinat]
            for i, vector_i in enumerate(paCMAP_coordinat):
                for j, vector_j in enumerate(paCMAP_coordinat):
                    if i < j:
                        dist = np.linalg.norm(np.array(vector_i) - np.array(vector_j))
                        distance_matrix[i][j] = dist
                        distance_matrix[j][i] = dist

        else:
            genomic_contexts = []
            for operon_type in Figures.operons:
                for i, target in enumerate(Figures.operons[operon_type]['target_members']):
                    if input_targets is None or target in input_targets:
                        genomic_contexts.append(Figures.operons[operon_type]['operon_protein_families_structure'][i])
                        labels.append(target)

            distance_matrix = [[0 for i in genomic_contexts] for j in genomic_contexts]
            for i, vector_i in enumerate(genomic_contexts):
                for j, vector_j in enumerate(genomic_contexts):
                    if i < j:
                        families_present = set(vector_i+vector_j)
                        overlapping_families = set(vector_i).intersection(set(vector_j))
                        dist = 1-len(overlapping_families)/len(families_present)
                        distance_matrix[i][j] = dist
                        distance_matrix[j][i] = dist

        return distance_matrix, labels    
    
    @staticmethod
    def get_avgoperons_distance_matrix(**kwargs) -> tuple[np.ndarray, list]:  
        Figures.set_class_attributes(Figures, **kwargs) 

        matrix = [[0 for i in Figures.operons if '-' not in i] for i in Figures.operons if '-' not in i]
        
        selected_operons = [i for i in Figures.operons if '-' not in i]
        selected_operons = sorted(selected_operons)
        for i, i_operon_type in enumerate(selected_operons):
            for j, j_operon_type in enumerate(selected_operons):
                if i > j:
                    dists = []
                    for i_member_pacmap in Figures.operons[i_operon_type]['operon_filtered_PaCMAP']:
                        for j_member_pacmap in Figures.operons[j_operon_type]['operon_filtered_PaCMAP']:
                            dist = np.linalg.norm(np.array(i_member_pacmap)-np.array(j_member_pacmap))
                            dists.append(dist)
                    dist = min(dists)
                    matrix[i][j] = dist
                    matrix[j][i] = dist        
        return np.array(matrix), selected_operons    
    
    @staticmethod
    def get_taxonomy_distance_matrix(**kwargs) -> tuple[np.ndarray, list]:
        Figures.set_class_attributes(Figures, **kwargs)
        
        targets_taxonomy_vector = {}

        for superkingdom in Figures.taxonomy:
            for phylum in Figures.taxonomy[superkingdom]:
                for taxclass in Figures.taxonomy[superkingdom][phylum]:
                    for order in Figures.taxonomy[superkingdom][phylum][taxclass]:
                        for genus in Figures.taxonomy[superkingdom][phylum][taxclass][order]:
                            for species in Figures.taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                for target in Figures.taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
                                    # find if this target is in the operons
                                    operon = [operon for operon in Figures.operons if target in Figures.operons[operon]['target_members']]
                                    if len(operon) > 0:
                                        targets_taxonomy_vector[target] = [superkingdom, phylum, taxclass, order, genus, species]

        if Figures.sort_mode == 'taxonomy':
            labels = sorted(list(targets_taxonomy_vector.keys()))
            taxonomy_distance_matrix = [[0 for label in targets_taxonomy_vector] 
                                        for label in targets_taxonomy_vector]

            for i, target_i in enumerate(labels):
                vector_i = targets_taxonomy_vector[target_i]
                for j, target_j in enumerate(labels):
                    vector_j = targets_taxonomy_vector[target_j]
                    if i >= j:
                        dist = 0
                        for k, level in enumerate(vector_i):
                            if level != vector_j[k]:
                                dist = 6-k
                                break
                        
                        if dist == 1: # it means they correspond to different species. check if they are just not different strains and fuse them if so
                            if len(list(set(vector_i[-1].split()) & set(vector_j[-1].split()))) >= 2: # it means they share the genus and the species at least
                                dist = 0						
                            
                        taxonomy_distance_matrix[i][j] = dist
                        taxonomy_distance_matrix[j][i] = dist

        elif Figures.sort_mode == 'as_input':
            labels = list(targets_taxonomy_vector.keys())
            taxonomy_distance_matrix = [[0 for label in targets_taxonomy_vector] 
                                        for label in targets_taxonomy_vector]

            if Figures.input_targets != None:
                labels = [in_target for in_target in Figures.input_targets if in_target in labels]
        
        taxonomy_distance_matrix = np.array(taxonomy_distance_matrix)        
        return taxonomy_distance_matrix, labels 

    @staticmethod
    def create_dendogram_features(**kwargs) -> tuple[list, dict]:
        Figures.set_class_attributes(Figures, **kwargs)

        if Figures.sort_mode == 'operon clusters':
            data = {'cluster size': [],
                    'x': [],
                    'y': [],
                    'color': [],
                    'leaf_label': Figures.leaf_labels}
            
            tooltips = [('Cluster type', '@leaf_label'),
                        ('Cluster size', '@cluster_size')]
            
        else:
            data = {'superkingdom': [],
                    'phylum': [],
                    'class': [],
                    'order': [],
                    'genus': [],
                    'species': [],
                    'x': [],
                    'y': [],
                    'color': [],
                    'leaf_label': Figures.leaf_labels}
            
            tooltips = [('Superkingdom', '@superkingdom'),
                        ('Phylum', '@phylum'),
                        ('Class', '@class'),
                        ('Order', '@order'),
                        ('Genus', '@genus'),
                        ('Species', '@species')]

        icoord, dcoord = Figures.dendogram['icoord'], Figures.dendogram['dcoord']

        data['y'] = list(np.linspace(min([num for sublist in icoord for num in sublist]), 
                                     max([num for sublist in icoord for num in sublist]), len(Figures.leaf_labels)))
        data['x'] = [1 for y in data['y']] 

        for label in Figures.leaf_labels:
            if Figures.colors is not None and label in Figures.colors:
                color = Figures.colors[label]['Color (RGBA)']
            else:
                color = 'darkgrey'
            
            data['color'].append(color)
            
            if Figures.sort_mode == 'operon clusters':
                data['cluster size'].append(str(len(Figures.operons[label]['target_members'])))
            else:
                found_taxonomy = False
                for superkingdom in Figures.taxonomy:
                    for phylum in Figures.taxonomy[superkingdom]:
                        for taxclass in Figures.taxonomy[superkingdom][phylum]:
                            for order in Figures.taxonomy[superkingdom][phylum][taxclass]:
                                for genus in Figures.taxonomy[superkingdom][phylum][taxclass][order]:
                                    for species in Figures.taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                        if species == label:
                                            data['superkingdom'].append(superkingdom)
                                            data['phylum'].append(phylum)
                                            data['class'].append(taxclass)
                                            data['order'].append(order)
                                            data['genus'].append(genus)
                                            found_taxonomy = True
                                        else:
                                            for target in Figures.taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
                                                if target == label:
                                                    data['superkingdom'].append(superkingdom)
                                                    data['phylum'].append(phylum)
                                                    data['class'].append(taxclass)
                                                    data['order'].append(order)
                                                    data['genus'].append(genus)
                                                    data['species'].append(species)
                                                    found_taxonomy = True
                if not found_taxonomy:
                    data['superkingdom'].append('na')
                    data['phylum'].append('na')
                    data['class'].append('na')
                    data['order'].append('na')
                    data['genus'].append('na')
                    data['species'].append('na')

        return tooltips, data            