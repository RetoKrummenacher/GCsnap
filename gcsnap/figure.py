import os
import numpy as np
import statistics
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial import distance

from collections import Counter

from Bio import Phylo

from bokeh.plotting import figure
from bokeh.models import HoverTool, TapTool
from bokeh.models.callbacks import OpenURL

class Figure:
    """ 
    Methods to create the figures used by AdvancedInteractiveFigure and InteractiveFigure.
    """

    @staticmethod
    def set_class_attributes(cls, **kwargs):
        """
        Set class attributes from kwargs.
        """        
        for key, value in kwargs.items():
            setattr(cls, key, value)  
    
    @staticmethod
    def make_dendogram_figure(**kwargs) -> tuple[figure, dict]:
        """
        Create the dendogram figure.

        Returns:
            tuple[figure, dict]: The figure and the data.
        """        
        # set attributes from kwargs, it updates any keywords that are passed to method
        Figure.set_class_attributes(Figure, **kwargs)

        
        dendogram, leaf_labels = Figure.compute_dendogram(**kwargs)
        den_tooltips, den_data = Figure.create_dendogram_features(dendogram = dendogram, 
                                                                leaf_labels = leaf_labels, **kwargs)

        y_range = [0, max(den_data['y'])+min(den_data['y'])+(den_data['y'][1]-den_data['y'][0])]

        if len(leaf_labels) < 5:
            height = int(len(leaf_labels)*Figure.height_factor)*3
        elif len(leaf_labels) < 10:
            height = int(len(leaf_labels)*Figure.height_factor)*2
        else:
            height = int(len(leaf_labels)*Figure.height_factor)

        if Figure.in_tree != None:
            title = 'Input phylogeny/hierarchy'
        elif Figure.sort_mode == 'taxonomy':
            title = 'Taxonomy hierarchy'
        elif Figure.sort_mode == 'operons':
            title = 'Genomic contexts similarity hierarchy'
        else:
            title = 'Genomic contexts clusters hierarchy'

        if Figure.show_leafs:
            den = figure(title = title, plot_height=height, plot_width=500, 
                         x_range=[-max(max(dendogram['dcoord']))-
                                  (max(max(dendogram['dcoord']))-min(min(dendogram['dcoord'])))*0.2, 40], 
                         y_range = y_range, toolbar_location="left")
        else:
            den = figure(title = title, plot_height=height, plot_width=250, 
                         x_range=[-max(max(dendogram['dcoord']))-
                                  (max(max(dendogram['dcoord']))-min(min(dendogram['dcoord'])))*0.2, 0.2], 
                         y_range = y_range, toolbar_location="left")

        for i, d in zip(dendogram['icoord'], dendogram['dcoord']):
            d = list(map(lambda x: -x, d))
            den.line(x=d, y=i, line_color='black')

        if Figure.show_leafs:
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
        """
        Compute the dendogram.

        Returns:
            tuple[dict, list]: The dendogram and the labels.
        """        
        Figure.set_class_attributes(Figure, **kwargs)

        distance_matrix = kwargs.get('distance_matrix')
        labels = kwargs.get('labels')
        
        if distance_matrix is None:
            if Figure.in_tree != None:
                distance_matrix, labels = Figure.get_phylogeny_distance_matrix(print_tree = True, **kwargs)
            elif Figure.sort_mode == 'operon':
                distance_matrix, labels = Figure.get_operons_distance_matrix(**kwargs)
            elif Figure.sort_mode == 'operon clusters':
                distance_matrix, labels = Figure.get_avgoperons_distance_matrix(**kwargs)		  
            else:
                distance_matrix, labels = Figure.get_taxonomy_distance_matrix(**kwargs)

        distance_matrix = distance.squareform(distance_matrix)

        Z = hierarchy.linkage(distance_matrix, method = 'single')
        results = hierarchy.dendrogram(Z, no_plot=True, count_sort = 'descending')

        if Figure.sort_mode == 'taxonomy' or Figure.in_tree != None or 'operon' in Figure.sort_mode:
            labels_indeces = list(map(int, results['ivl']))
            labels = [labels[i] for i in labels_indeces]

        return results, labels    
    
    @staticmethod
    def get_phylogeny_distance_matrix(**kwargs) -> tuple[np.ndarray, list]:
        """
        Get the phylogeny distance matrix.

        Returns:
            tuple[np.ndarray, list]: The distance matrix and the labels.
        """        
        Figure.set_class_attributes(Figure, **kwargs)

        if Figure.starting_directory not in Figure.in_tree:
            Figure.in_tree = os.path.join(Figure.starting_directory,'{}'.format(Figure.in_tree))

        tree = Phylo.read(Figure.in_tree, Figure.in_tree_format)  

        if Figure.print_tree:
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
            if leaf not in Figure.input_targets:
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
        """
        Get the operons distance matrix.

        Returns:
            tuple[np.ndarray, list]: The distance matrix and the labels.
        """        
        # set attributes from kwargs, it updates any keywords that are passed to method
        Figure.set_class_attributes(Figure, **kwargs)       

        input_targets = kwargs.get('input_targets', None)

        labels = []
        if 'operon_filtered_PaCMAP' in Figure.operons[list(Figure.operons.keys())[0]]:
            paCMAP_coordinat = []
            for operon_type in Figure.operons:
                for i, target in enumerate(Figure.operons[operon_type]['target_members']):
                    if input_targets is None or target in input_targets:
                        paCMAP_coordinat.append(Figure.operons[operon_type]['operon_filtered_PaCMAP'][i])
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
            for operon_type in Figure.operons:
                for i, target in enumerate(Figure.operons[operon_type]['target_members']):
                    if input_targets is None or target in input_targets:
                        genomic_contexts.append(Figure.operons[operon_type]['operon_protein_families_structure'][i])
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
        """
        Get the average operons distance matrix.

        Returns:
            tuple[np.ndarray, list]: The distance matrix and the labels.
        """        
        Figure.set_class_attributes(Figure, **kwargs) 

        matrix = [[0 for i in Figure.operons if '-' not in i] for i in Figure.operons if '-' not in i]
        
        selected_operons = [i for i in Figure.operons if '-' not in i]
        selected_operons = sorted(selected_operons)
        for i, i_operon_type in enumerate(selected_operons):
            for j, j_operon_type in enumerate(selected_operons):
                if i > j:
                    dists = []
                    for i_member_pacmap in Figure.operons[i_operon_type]['operon_filtered_PaCMAP']:
                        for j_member_pacmap in Figure.operons[j_operon_type]['operon_filtered_PaCMAP']:
                            dist = np.linalg.norm(np.array(i_member_pacmap)-np.array(j_member_pacmap))
                            dists.append(dist)
                    dist = min(dists)
                    matrix[i][j] = dist
                    matrix[j][i] = dist        
        return np.array(matrix), selected_operons    
    
    @staticmethod
    def get_taxonomy_distance_matrix(**kwargs) -> tuple[np.ndarray, list]:
        """
        Get the taxonomy distance matrix.

        Returns:
            tuple[np.ndarray, list]: The distance matrix and the labels.
        """        
        Figure.set_class_attributes(Figure, **kwargs)
        
        targets_taxonomy_vector = {}

        for superkingdom in Figure.taxonomy:
            for phylum in Figure.taxonomy[superkingdom]:
                for taxclass in Figure.taxonomy[superkingdom][phylum]:
                    for order in Figure.taxonomy[superkingdom][phylum][taxclass]:
                        for genus in Figure.taxonomy[superkingdom][phylum][taxclass][order]:
                            for species in Figure.taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                for target in Figure.taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
                                    # find if this target is in the operons
                                    operon = [operon for operon in Figure.operons if target in Figure.operons[operon]['target_members']]
                                    if len(operon) > 0:
                                        targets_taxonomy_vector[target] = [superkingdom, phylum, taxclass, order, genus, species]

        if Figure.sort_mode == 'taxonomy':
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

        elif Figure.sort_mode == 'as_input':
            labels = list(targets_taxonomy_vector.keys())
            taxonomy_distance_matrix = [[0 for label in targets_taxonomy_vector] 
                                        for label in targets_taxonomy_vector]

            if Figure.input_targets != None:
                labels = [in_target for in_target in Figure.input_targets if in_target in labels]
        
        taxonomy_distance_matrix = np.array(taxonomy_distance_matrix)        
        return taxonomy_distance_matrix, labels 

    @staticmethod
    def create_dendogram_features(**kwargs) -> tuple[list, dict]:
        """
        Create the dendogram features.

        Returns:
            tuple[list, dict]: The tooltips and the data.
        """        
        Figure.set_class_attributes(Figure, **kwargs)

        if Figure.sort_mode == 'operon clusters':
            data = {'cluster size': [],
                    'x': [],
                    'y': [],
                    'color': [],
                    'leaf_label': Figure.leaf_labels}
            
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
                    'leaf_label': Figure.leaf_labels}
            
            tooltips = [('Superkingdom', '@superkingdom'),
                        ('Phylum', '@phylum'),
                        ('Class', '@class'),
                        ('Order', '@order'),
                        ('Genus', '@genus'),
                        ('Species', '@species')]

        icoord, dcoord = Figure.dendogram['icoord'], Figure.dendogram['dcoord']

        data['y'] = list(np.linspace(min([num for sublist in icoord for num in sublist]), 
                                     max([num for sublist in icoord for num in sublist]), len(Figure.leaf_labels)))
        data['x'] = [1 for y in data['y']] 

        for label in Figure.leaf_labels:
            if Figure.colors is not None and label in Figure.colors:
                color = Figure.colors[label]['Color (RGBA)']
            else:
                color = 'darkgrey'
            
            data['color'].append(color)
            
            if Figure.sort_mode == 'operon clusters':
                data['cluster size'].append(str(len(Figure.operons[label]['target_members'])))
            else:
                found_taxonomy = False
                for superkingdom in Figure.taxonomy:
                    for phylum in Figure.taxonomy[superkingdom]:
                        for taxclass in Figure.taxonomy[superkingdom][phylum]:
                            for order in Figure.taxonomy[superkingdom][phylum][taxclass]:
                                for genus in Figure.taxonomy[superkingdom][phylum][taxclass][order]:
                                    for species in Figure.taxonomy[superkingdom][phylum][taxclass][order][genus]:
                                        if species == label:
                                            data['superkingdom'].append(superkingdom)
                                            data['phylum'].append(phylum)
                                            data['class'].append(taxclass)
                                            data['order'].append(order)
                                            data['genus'].append(genus)
                                            found_taxonomy = True
                                        else:
                                            for target in Figure.taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
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

    @staticmethod
    def create_most_common_genomic_features_figure(**kwargs) -> figure:
        """
        Create the most common genomic features figure.

        Returns:
            figure: The figure.
        """        
        Figure.set_class_attributes(Figure, **kwargs)

        most_common_context = Figure.find_most_common_genomic_context(**kwargs)
        gc_tooltips, gc_data = Figure.create_most_common_genomic_context_features(
                            most_common_context = most_common_context, **kwargs)

        gc = figure(plot_width=2000, plot_height=200, y_range = [0, 4], 
                    title = 'Most conserved gene per position', toolbar_location="left")

        for i, xs in enumerate(gc_data['xs']):
            gc.patch(xs, gc_data['ys'][i], 
                     fill_color = gc_data['facecolor'][i], 
                     line_color = gc_data['edgecolor'][i], 
                     fill_alpha = gc_data['transparency'][i], 
                     line_alpha = gc_data['transparency'][i], line_width = 1)	
        
        gc.patches('xs', 'ys', fill_color = None, line_color = None, line_width = 0, source = gc_data, 
                hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
                selection_fill_color='facecolor', selection_line_color='edgecolor',
                nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', 
                nonselection_fill_alpha=0.2)

        gc.text('text_x', 'text_y', text = 'family', text_baseline="bottom", 
                text_align="center", text_font_size = {'value': '6pt'}, source = gc_data)
        gc.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", 
                text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, 
                source = gc_data)
        
        gc.yaxis.ticker = [1]
        gc.yaxis.major_label_overrides = {1: max([Figure.syntenies[target]['species'] 
                                                  for target in Figure.syntenies], key=len)}
        # gc.yaxis.major_label_text_font_size = {'value': '8pt'}
        
        gc.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
        gc.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
        gc.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
        gc.yaxis.axis_line_width = 0

        # define general features
        gc.grid.visible = False
        gc.outline_line_width = 0
        
        # define xticks
        gc.xaxis.axis_label = "Position relative to target (bp)"
        
        gc.add_tools(HoverTool(tooltips=gc_tooltips))
        gc.add_tools(TapTool(callback = OpenURL(url='@model_links')))
        
        return gc
        
    @staticmethod
    def find_most_common_genomic_context(**kwargs) -> dict:
        """
        Find the most common genomic context.

        Returns:
            dict: The most common genomic context.
        """        
        Figure.set_class_attributes(Figure, **kwargs)

        # will use only the complete genomic contexts and ignore the partial ones	
        operon_matrix = []
        for operon in Figure.operons:
            for curr_context in Figure.operons[operon]['operon_protein_families_structure']:
                if len(curr_context) == Figure.n_flanking5 + Figure.n_flanking3 + 1:
                    operon_matrix.append(curr_context)
        
        operon_matrix = np.array(operon_matrix).T
        
        most_common_context = {'selected_context': [],
                            'families_frequency': [],
                            'average_starts': [],
                            'average_ends': [],
                            'average_size': [],
                            'stdev_size': [],
                            'directions': [],
                            'tm_annotations': []}
        
        for i, column in enumerate(operon_matrix):
            occurence_count = Counter(column) 
            most_common_family = occurence_count.most_common(1)[0][0]

            most_common_context['selected_context'].append(most_common_family)
            most_common_context['families_frequency'].append(round(
                occurence_count.most_common(1)[0][1]*100/len(column), 1))

            all_starts_of_most_common = []
            all_ends_of_most_common = []
            all_orientations = []
            all_sizes = []
            all_tm_annotations = []

            for operon in Figure.operons:
                for j, curr_context in enumerate(Figure.operons[operon]['operon_protein_families_structure']):
                    curr_target = Figure.operons[operon]['target_members'][j]
                
                    if len(curr_context) == Figure.n_flanking5 + Figure.n_flanking3 + 1:			
                        if Figure.operons[operon]['operon_protein_families_structure'][j][i] == most_common_family:
                            all_starts_of_most_common.append(Figure.syntenies[curr_target]['flanking_genes']
                                                             ['relative_starts'][i])
                            all_ends_of_most_common.append(Figure.syntenies[curr_target]['flanking_genes']
                                                           ['relative_ends'][i])
                            all_sizes.append((Figure.syntenies[curr_target]['flanking_genes']
                                              ['relative_ends'][i] - Figure.syntenies[curr_target]
                                              ['flanking_genes']['relative_starts'][i])/3)
                            all_orientations.append(Figure.syntenies[curr_target]['flanking_genes']['directions'][i])
                            
                            if 'TM_annotations' in Figure.syntenies[curr_target]['flanking_genes']:
                                all_tm_annotations.append(Figure.syntenies[curr_target]['flanking_genes']
                                                          ['TM_annotations'][i])
            
            most_common_context['average_starts'].append(int(statistics.median(all_starts_of_most_common)))
            most_common_context['average_ends'].append(int(statistics.median(all_ends_of_most_common)))
            most_common_context['average_size'].append(int(statistics.median(all_sizes)))
            most_common_context['stdev_size'].append(int(stats.median_abs_deviation(all_sizes)))
            
            try:
                most_common_context['directions'].append(statistics.mode(all_orientations))
            except:
                most_common_context['directions'].append('+')
            
            try:
                most_common_context['tm_annotations'].append(statistics.mode(all_tm_annotations))
            except:
                most_common_context['tm_annotations'].append('')
            
        return most_common_context
            
    @staticmethod
    def create_most_common_genomic_context_features(**kwargs) -> tuple:
        """
        Create the most common genomic context features.

        Returns:
            tuple: The tooltips and the data.
        """        
        Figure.set_class_attributes(Figure, **kwargs)

        data = {'xs': [],
                    'ys': [],
                    'edgecolor': [],
                    'facecolor': [],
                    'text_x': [],
                    'text_y': [],
                    'family': [],
                    'tm_text_x': [],
                    'tm_text_y': [],
                    'tm_text': [],
                    'tm_pred_text': [],
                    'protein_name': [],
                    'protein_size': [],
                    'family_frequency': [],
                    'transparency': [],
                    'relative_start': [],
                    'relative_end': [],
                    'found_models': [],
                    'model_links': []}
            
        for i, family in enumerate(Figure.most_common_context['selected_context']):
            gene_dx = Figure.most_common_context['average_ends'][i] - \
                        Figure.most_common_context['average_starts'][i]+1
            gene_direction = Figure.most_common_context['directions'][i]

            if gene_direction == '-':
                gene_x_tail = Figure.most_common_context['average_ends'][i]
                gene_dx = gene_dx*(-1)
                gene_x_head = gene_x_tail + gene_dx
                gene_x_head_start = gene_x_head+100
                text_x = gene_x_tail - (gene_x_tail-gene_x_head_start)/2
            else:
                gene_x_tail = Figure.most_common_context['average_starts'][i]
                gene_x_head = gene_x_tail + gene_dx
                gene_x_head_start = gene_x_head-100
                text_x = gene_x_tail + (gene_x_head_start-gene_x_tail)/2

            if family == 0:
                facecolor = Figure.family_colors[family]['Color (RGBA)']
                edgecolor = Figure.family_colors[family]['Line color']
                linestyle = Figure.family_colors[family]['Line style']
            else:
                facecolor = Figure.family_colors[family]['Color (RGBA)']
                edgecolor = Figure.family_colors[family]['Line color']
                linestyle = Figure.family_colors[family]['Line style'] 
            
            if family == 0:
                relative_start = 'n.a.'
                relative_end = 'n.a.'
                family_frequency = 'n.a.'
                protein_size = 'n.a.' 
                protein_name = 'n.a.'
                transparency = 0.2
                tm_annotation = ''
                tm_pred_text = 'n.a.'
            else:
                relative_start = format(gene_x_tail, ',d')
                relative_end = format(gene_x_head, ',d')
                family_frequency = '{}%'.format(Figure.most_common_context['families_frequency'][i])
                protein_size = r'{} ({})'.format(Figure.most_common_context['average_size'][i], 
                                                 Figure.most_common_context['stdev_size'][i])
                protein_name = Figure.families_summary[family]['name']
                transparency = Figure.most_common_context['families_frequency'][i]/100  
                
                if 'tm_annotations' in Figure.most_common_context:
                    tm_annotation = Figure.most_common_context['tm_annotations'][i]
                    if tm_annotation == 'TM':
                        tm_pred_text = 'Yes'
                    elif tm_annotation == 'SP':
                        tm_pred_text = 'Contains signal peptide'
                    else:
                        tm_pred_text = 'No'
                else:
                    tm_annotation = ''
                    tm_pred_text = 'n.a.' 
                    
            data['relative_start'].append(relative_start)
            data['relative_end'].append(relative_end)
            data['facecolor'].append(facecolor)
            data['edgecolor'].append(edgecolor)
            data['family_frequency'].append(family_frequency)
            data['protein_size'].append(protein_size)
            data['protein_name'].append(protein_name)
            data['xs'].append([gene_x_tail, gene_x_tail, gene_x_head_start, gene_x_head, gene_x_head_start])
            data['ys'].append([1-0.25, 1+0.25, 1+0.25, 1, 1-0.25])
            data['text_x'].append(text_x)
            data['text_y'].append(1+0.25)
            data['transparency'].append(transparency)
            
            data['tm_text'].append(tm_annotation)
            data['tm_text_x'].append(text_x)
            data['tm_text_y'].append(1)
            data['tm_pred_text'].append(tm_pred_text)

            if family > 0 and family != Figure.reference_family:
                data['family'].append(family)
            else:
                data['family'].append(str(''))
            
            if 'model_state' in Figure.families_summary[family]:
                model_state = Figure.families_summary[family]['model_state']

                if model_state == 'Model exists':
                    model_state = 'Yes (click to view in Swiss-Model repository)'
                elif model_state == 'Model does not exist':
                    model_state = 'No (click to model with Swiss-Model)'
                else:
                    if family > 0:
                        model_state = 'Not possible to find'
                    else:
                        model_state = ''
                
                structure = Figure.families_summary[family]['structure']
                if structure == '':
                    uniprot_code = Figure.families_summary[family]['uniprot_code']
                    structure = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
                
            else:
                model_state = 'n.a.'
                structure = 'n.a.'
            
            data['found_models'].append(model_state)
            data['model_links'].append(structure)
            
        tooltips = [('Protein name', "@protein_name"),
                    ("Predicted membrane protein", "@tm_pred_text"),
                    ('Structural model found', '@found_models'),
                    ('Frequency in position', '@family_frequency'),
                    ('Median protein size', '@protein_size'),
                    ('Median starting position', '@relative_start'),
                    ('Median end position', '@relative_end')] 
        
        return tooltips, data      

    @staticmethod
    def create_genomic_context_figure(**kwargs) -> figure:
        """
        Create the genomic context figure.

        Returns:
            figure: The figure.
        """        
        Figure.set_class_attributes(Figure, **kwargs)

        p_tooltips, p_data, p_yyticklabels = Figure.create_genomic_context_features(**kwargs)

        # the genomic_context figure
        p = figure(plot_width=Figure.most_common_gc_figure.plot_width, 
                   plot_height=Figure.syn_dendogram.plot_height, 
                   x_range = Figure.most_common_gc_figure.x_range, 
                   y_range = Figure.syn_dendogram.y_range, 
                   toolbar_location="left", 
                   title = 'Representative genomic contexts (hover to get more information)')    
       
        for i, xs in enumerate(p_data['xs']):
            p.patch(xs, p_data['ys'][i], fill_color = p_data['facecolor'][i], 
                    line_color = p_data['edgecolor'][i], line_width = 1)

        p.patches('xs', 'ys', fill_color = None, line_color = None, line_width = 0, source = p_data, 
                    hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
                    selection_fill_color='facecolor', selection_line_color='edgecolor',
                    nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', 
                    nonselection_fill_alpha=0.2)

        p.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", 
               text_font_size = {'value': '6pt'}, source = p_data)
        p.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", text_baseline="middle", 
               text_align="center", text_font_size = {'value': '6pt'}, source = p_data)
        
        # define yticks on the left
        p.yaxis.ticker = [int(n) for n in list(p_yyticklabels.keys())]
        p.yaxis.major_tick_line_color = None
        p.yaxis.major_label_overrides = {int(i): p_yyticklabels[i] for i in [int(n) for n in p_yyticklabels.keys()]}
        #	 p.yaxis.major_label_text_font_size = {'value': '8pt'}
        #	 p.yaxis.major_label_text_font_style = 'italic'
        p.yaxis.axis_line_width = 0

        # define xticks
        p.xaxis.axis_label = "Position relative to target (bp)"

        # define general features
        p.grid.visible = False
        p.outline_line_width = 0

        p.add_tools(HoverTool(tooltips=p_tooltips))
        p.add_tools(TapTool(callback = OpenURL(url='@ncbi_link')))
            
        return p        
    
    @staticmethod
    def create_genomic_context_features(**kwargs) -> tuple[list, dict, dict]:
        """
        Create the genomic context features.

        Returns:
            tuple[list, dict, dict]: The tooltips, the data, and the yyticklabels.
        """        
        Figure.set_class_attributes(Figure, **kwargs)
        
        data = {'operon':[],
                'target_id':[],
                'ncbi_code':[],
                'ncbi_link':[],
                'assembly': [],
                'name':[],
                'family': [],
                'superkingdom': [],
                'phylum': [],
                'class': [],
                'order': [],
                'genus': [],
                'species': [],
                'relative_start': [],
                'relative_end': [],
                'facecolor': [],
                'edgecolor': [],
                'linestyle': [],
                'xs':[],
                'ys':[],
                'half_heights': [],
                'text_x': [],
                'text_y': [],
                'tm_text_x': [],
                'tm_text_y': [],
                'tm_text': [],
                'tm_pred_text': []}
        
        yyticklabels = []
        yys = []
        y_step = Figure.syn_den_data['y'][1] - Figure.syn_den_data['y'][0]
        y_half_height = y_step/4

        for i, current_target in enumerate(Figure.syn_den_data['leaf_label']):
            operon = [operon for operon in Figure.operons if current_target in Figure.operons[operon]['target_members']][0]
            current_assembly = Figure.syntenies[current_target]['assembly_id'][1]
            current_species = Figure.syntenies[current_target]['species']
            current_genomic_context_block = Figure.syntenies[current_target]['flanking_genes']
            current_species = Figure.syntenies[current_target]['species']
            current_reference_family = Figure.syntenies[current_target]['target_family']
            
            curr_y = Figure.syn_den_data['y'][i]

            for j, flanking_gene in enumerate(current_genomic_context_block['ncbi_codes']):
                family = current_genomic_context_block['families'][j]
                name = current_genomic_context_block['names'][j]
                gene_dx = current_genomic_context_block['relative_ends'][j] - \
                        current_genomic_context_block['relative_starts'][j]+1
                gene_direction = current_genomic_context_block['directions'][j]

                if gene_direction == '-':
                    gene_x_tail = current_genomic_context_block['relative_ends'][j]
                    gene_dx = gene_dx*(-1)
                    gene_x_head = gene_x_tail + gene_dx
                    gene_x_head_start = gene_x_head+100
                    text_x = gene_x_tail - (gene_x_tail-gene_x_head_start)/2
                else:
                    gene_x_tail = current_genomic_context_block['relative_starts'][j]
                    gene_x_head = gene_x_tail + gene_dx
                    gene_x_head_start = gene_x_head-100
                    text_x = gene_x_tail + (gene_x_head_start-gene_x_tail)/2

                if family == 0:
                    facecolor = Figure.family_colors[family]['Color (RGBA)']
                    edgecolor = Figure.family_colors[family]['Line color']
                    linestyle = Figure.family_colors[family]['Line style']
                else:
                    facecolor = Figure.family_colors[family]['Color (RGBA)']
                    edgecolor = Figure.family_colors[family]['Line color']
                    linestyle = Figure.family_colors[family]['Line style'] 
                    
                data['operon'].append(operon)
                data['target_id'].append(current_target)
                data['ncbi_code'].append(flanking_gene)
                data['assembly'].append(current_assembly)
                data['name'].append(name)
                data['relative_start'].append(format(gene_x_tail, ',d'))
                data['relative_end'].append(format(gene_x_head, ',d'))
                data['facecolor'].append(facecolor)
                data['edgecolor'].append(edgecolor)
                data['linestyle'].append(linestyle)
                data['xs'].append([gene_x_tail, gene_x_tail, gene_x_head_start, gene_x_head, gene_x_head_start])
                data['ys'].append([curr_y-y_half_height, curr_y+y_half_height, curr_y+y_half_height, curr_y, curr_y-y_half_height])
                data['half_heights'].append(y_half_height)
                data['text_x'].append(text_x)
                data['text_y'].append(curr_y+y_half_height)
                data['ncbi_link'].append('https://www.ncbi.nlm.nih.gov/protein/{}'.format(flanking_gene))
                data['tm_text_x'].append(text_x)
                data['tm_text_y'].append(curr_y)
                
                if 'TM_annotations' in current_genomic_context_block:
                    tm_annotation = current_genomic_context_block['TM_annotations'][j]
                    if tm_annotation == 'TM':
                        data['tm_pred_text'].append('Yes')
                    elif tm_annotation == 'SP':
                        data['tm_pred_text'].append('Contains signal peptide')
                    else:
                        data['tm_pred_text'].append('No')
                else:
                    tm_annotation = ''
                    data['tm_pred_text'].append('n.a.')
                    
                data['tm_text'].append(tm_annotation)
                                    
                if family > 0 and family != Figure.reference_family:
                    data['family'].append(family)
                else:
                    data['family'].append(str(''))
                    
                data['species'].append(current_species)
                for level in ['superkingdom', 'phylum', 'class', 'order', 'genus']:
                    data[level].append(Figure.syn_den_data[level][i])
            
            if Figure.gc_legend_mode in ['superkingdom', 'phylum', 'class', 'order', 'genus', 'species']:
                label = data[Figure.gc_legend_mode][-1]
                if label == []:
                    label = 'n.a'
                yyticklabels.append(label)
            else:
                yyticklabels.append('{} | {}'.format(current_target, operon))
                
            yys.append(curr_y)
            
        yyticklabels = {int(yys[i]): yyticklabels[i] for i in range(len(yyticklabels))}

        tooltips = [('GC type', "@operon"),
                    ('InputID', "@target_id"),
                    ('EntrezID', '@ncbi_code'),
                    ('Genome assembly', '@assembly'),
                    ('Gene relative start', '@relative_start'),
                    ('Gene relative end', '@relative_end'),
                    ("Protein name", "@name"),
                    ("Protein family code", "@family"),
                    ("Predicted membrane protein", "@tm_pred_text"),
                    ('Superkingdom', '@superkingdom'),
                    ('Phylum', '@phylum'),
                    ('Class', '@class'),
                    ('Order', '@order'),
                    ('Genus', '@genus'),
                    ("Species", "@species")] 
            
        return tooltips, data, yyticklabels                  