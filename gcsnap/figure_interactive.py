import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import statistics
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial import distance

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
from gcsnap.figures import Figures
from gcsnap.rich_console import RichConsole

class InteractiveFigure:
    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str, ref_family: str, 
                 family_colors: dict, input_targets: list):

        # Extract parameters from config and gc
        kwargs = {k: v['value'] for k, v in config.arguments.items()}
        kwargs.update({
            'config': config,
            'out_label': out_label,
            'gc': gc,
            'operons': gc.get_selected_operons(),
            'most_populated_operon': gc.get_most_populated_operon(),
            'synthenise': gc.get_synthenise(),
            'families_summary': gc.get_families(),
            'taxonomy': gc.get_taxonomy(),
            'reference_family': ref_family,
            'family_colors': family_colors,
            'input_targets': input_targets,
        })
        # change sort_mode and input targets if needed
        if kwargs['in_tree'] is not None:
            kwargs['sort_mode'] = 'tree'    
            kwargs['input_targets'] = [target for operon in self.operons 
                                       for target in self.operons[operon]['target_members']]                

        # Set all attributes
        self._set_attributes(**kwargs)

        self.console = RichConsole()

    def _set_attributes(self, **kwargs):
        """Helper method to set attributes from kwargs."""
        for key, value in kwargs.items():
            # Only set the attribute if it does not already exist
            if not hasattr(self, key):
                setattr(self, key, value)   

    def run(self) -> None:
        with self.console.status('Creating interactive genomic context figures'):
            # Make a copy of the current object's attributes (all contained in instance __dict__)
            kwargs = self.__dict__.copy()
            # output to static HTML file
            output_file(os.path.join(os.getcwd,'{}_interactive_output.html'.format(self.out_label)))
            # Work on most conserved genomic context figure        
            most_common_gc_figure = self.create_most_common_genomic_features_figure(**kwargs) 
            # Work on gene co-occurence figure        
            coocurrence_figure, graph_coord = self.create_graph_figure(mgc = most_common_gc_figure, 
                                                                    graph_coord = {},
                                                                    mode = 'coocurrence',
                                                                    previous_net = '',
                                                                    **kwargs) 
            adjacency_figure, graph_coord = self.create_graph_figure(mode = 'adjacency', 
                                                            graph_coord=graph_coord, 
                                                            previous_net=coocurrence_figure, **kwargs)
            # Work on dendogram for the genomic context block
            syn_dendogram, syn_den_data = Figures.make_dendogram_figure(show_leafs = False,                                                                        
                                                                height_factor = 25*1.2, 
                                                                distance_matrix = None, 
                                                                labels = None, 
                                                                colors = None,
                                                                **kwargs)
            
            genomic_context_figure = self.create_genomic_context_figure(syn_dendogram = syn_dendogram,
                                                                    syn_den_data = syn_den_data,
                                                                    **kwargs)

            legend_figure = self.create_legend_figure(grid = genomic_context_figure,
                                                    rescale_height = False,
                                                    **kwargs)
            
            # Make the tabs for the network and most common genomic context
            tab1 = Panel(child=coocurrence_figure, title='Gene co-occurrence network')
            tab2 = Panel(child=adjacency_figure, title='Gene adjacency network')
            tab3 = Panel(child=most_common_gc_figure, title='Most common genomic features')
            tabs = Tabs(tabs=[tab1, tab2, tab3])

            # Make a grid out of them
            grid = gridplot([[None, tabs, None], 
                                [syn_dendogram, genomic_context_figure, legend_figure]], merge_tools = True)

            save(grid)

    def create_most_common_genomic_features_figure(self, **kwargs) -> figure:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        most_common_context = self.find_most_common_genomic_context(**kwargs)
        gc_tooltips, gc_data = self.create_most_common_genomic_context_features(most_common_context, **kwargs)

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
        gc.yaxis.major_label_overrides = {1: max([self.syntenies[target]['species'] 
                                                  for target in self.syntenies], key=len)}
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
        
    def find_most_common_genomic_context(self, **kwargs) -> dict:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)        

        # will use only the complete genomic contexts and ignore the partial ones	
        operon_matrix = []
        for operon in self.operons:
            for curr_context in self.operons[operon]['operon_protein_families_structure']:
                if len(curr_context) == self.n_flanking5 + self.n_flanking3 + 1:
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

            for operon in self.operons:
                for j, curr_context in enumerate(self.operons[operon]['operon_protein_families_structure']):
                    curr_target = self.operons[operon]['target_members'][j]
                
                    if len(curr_context) == self.n_flanking5 + self.n_flanking3 + 1:			
                        if self.operons[operon]['operon_protein_families_structure'][j][i] == most_common_family:
                            all_starts_of_most_common.append(self.syntenies[curr_target]['flanking_genes']
                                                             ['relative_starts'][i])
                            all_ends_of_most_common.append(self.syntenies[curr_target]['flanking_genes']
                                                           ['relative_ends'][i])
                            all_sizes.append((self.syntenies[curr_target]['flanking_genes']
                                              ['relative_ends'][i] - self.syntenies[curr_target]
                                              ['flanking_genes']['relative_starts'][i])/3)
                            all_orientations.append(self.syntenies[curr_target]['flanking_genes']['directions'][i])
                            
                            if 'TM_annotations' in self.syntenies[curr_target]['flanking_genes']:
                                all_tm_annotations.append(self.syntenies[curr_target]['flanking_genes']
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
            
    def create_most_common_genomic_context_features(self, **kwargs) -> tuple:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)
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
            
        for i, family in enumerate(self.most_common_context['selected_context']):
            gene_dx = self.most_common_context['average_ends'][i] - \
                        self.most_common_context['average_starts'][i]+1
            gene_direction = self.most_common_context['directions'][i]

            if gene_direction == '-':
                gene_x_tail = self.most_common_context['average_ends'][i]
                gene_dx = gene_dx*(-1)
                gene_x_head = gene_x_tail + gene_dx
                gene_x_head_start = gene_x_head+100
                text_x = gene_x_tail - (gene_x_tail-gene_x_head_start)/2
            else:
                gene_x_tail = self.most_common_context['average_starts'][i]
                gene_x_head = gene_x_tail + gene_dx
                gene_x_head_start = gene_x_head-100
                text_x = gene_x_tail + (gene_x_head_start-gene_x_tail)/2

            if family == 0:
                facecolor = self.family_colors[family]['Color (RGBA)']
                edgecolor = self.family_colors[family]['Line color']
                linestyle = self.family_colors[family]['Line style']
            else:
                facecolor = self.family_colors[family]['Color (RGBA)']
                edgecolor = self.family_colors[family]['Line color']
                linestyle = self.family_colors[family]['Line style'] 
            
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
                family_frequency = '{}%'.format(self.most_common_context['families_frequency'][i])
                protein_size = r'{} ({})'.format(self.most_common_context['average_size'][i], 
                                                 self.most_common_context['stdev_size'][i])
                protein_name = self.families_summary[family]['name']
                transparency = self.most_common_context['families_frequency'][i]/100  
                
                if 'tm_annotations' in self.most_common_context:
                    tm_annotation = self.most_common_context['tm_annotations'][i]
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

            if family != 0 and family != self.reference_family and family < 10000:
                data['family'].append(family)
            else:
                data['family'].append(str(''))
            
            if 'model_state' in self.families_summary[family]:
                model_state = self.families_summary[family]['model_state']

                if model_state == 'Model exists':
                    model_state = 'Yes (click to view in Swiss-Model repository)'
                elif model_state == 'Model does not exist':
                    model_state = 'No (click to model with Swiss-Model)'
                else:
                    if family > 0 and family < 10000:
                        model_state = 'Not possible to find'
                    else:
                        model_state = ''
                
                structure = self.families_summary[family]['structure']
                if structure == '':
                    uniprot_code = self.families_summary[family]['uniprot_code']
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

    def create_graph_figure(self, **kwargs) -> tuple[nx.Graph, dict]:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        matrix, selected_families_summary = self.get_coocurrence_matrix(**kwargs)
        graph = self.get_graph_from_matrix(matrix, selected_families_summary, **kwargs)
        
        if self.mode == 'coocurrence':
            title = 'Gene co-occurrence network'
        elif self.mode == 'adjacency':
            graph = self.remove_non_adjacent_edges(graph=graph, 
                                    families_present = sorted(list(selected_families_summary.keys())),
                                    **kwargs)
            title = 'Gene adjcency network'
        
        if len(graph_coord) == 0:
            graph_coord = nx.spring_layout(graph)
        
        node_data, node_tooltips = self.create_node_features(node_graph_coord=graph_coord, graph=graph, **kwargs)
            
        if self.previous_net != '':
            x_range = self.previous_net.x_range
            y_range = self.previous_net.y_range
            
        else:
            x_range = (min([graph_coord[node][0] for node in graph_coord])-0.5, 
                       max([graph_coord[node][0] for node in graph_coord])+0.5)
            y_range = (min([graph_coord[node][1] for node in graph_coord])-0.5, 
                       max([graph_coord[node][1] for node in graph_coord])+0.5)

        g = figure(width = self.mgc.plot_width, height = self.mgc.plot_height, x_range=x_range, y_range=y_range, title = title)

        graph_renderer = from_networkx(graph, graph_coord, scale=1, center=(0, 0))
        graph_renderer.edge_renderer.glyph = MultiLine(line_width="line_width", line_color = "edge_color")
        graph_renderer.node_renderer.glyph = Circle(size=22, fill_color = "node_color")
        
        g.renderers.append(graph_renderer)
        
        g.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", text_font_size = {'value': '6pt'}, source = node_data)
        g.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, source = node_data)
        g.circle('tm_text_x', 'tm_text_y', color = None, size = 22, source = node_data)
        
        g.add_tools(HoverTool(tooltips=node_tooltips))
        g.add_tools(TapTool(callback = OpenURL(url='@model_links')))

        g.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
        g.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
        g.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
        g.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
        g.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
        g.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
        g.yaxis.axis_line_width = 0
        g.xaxis.axis_line_width = 0
        # define general features
        g.grid.visible = False
        g.outline_line_width = 0
        
        return g, graph_coord        
    
    def get_coocurrence_matrix(self, **kwargs) -> tuple[np.ndarray, dict]:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        family_labels = sorted([family for family in self.families_summary.keys() 
                                if family > 0 and family < 10000])
        matrix = [[0 for family in family_labels] for family in family_labels]

        context_count = 0
        for operon in self.operons:
            for genomic_context in self.operons[operon]['operon_protein_families_structure']:
                for i in range(len(set(genomic_context))):
                    for j in range(len(set(genomic_context))):
                        if i > j:
                            out_family = list(set(genomic_context))[i]
                            in_family = list(set(genomic_context))[j]

                            if all(family > 0 and family < 10000 for family in [out_family, in_family]):
                                matrix[family_labels.index(out_family)][family_labels.index(in_family)] += 1
                                matrix[family_labels.index(in_family)][family_labels.index(out_family)] += 1


        matrix = np.array(matrix)
        matrix = matrix/np.amax(matrix)
        matrix = np.where(matrix < self.min_coocc, 0, matrix)
        matrix = np.where(matrix!=0, (np.exp(matrix*2)-1)*1, matrix)

        # remove all columns and lines with all zero
        indices_with_all_zeros = np.all(matrix == 0, axis=1)

        matrix = matrix[~indices_with_all_zeros]
        matrix = matrix.T[~indices_with_all_zeros]

        selected_families = np.array(family_labels)[~indices_with_all_zeros]
        selected_families_summary = {family: self.families_summary[family] for family in selected_families}
        
        return matrix, selected_families_summary        
    
    def get_graph_from_matrix(self, **kwargs) -> nx.Graph:  
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        G=nx.from_numpy_array(self.matrix)

        # take care of the edges
        edge_params = {'color': {}, 'weight': {}, 'line_width': {}}
        edge_cmap = matplotlib.cm.get_cmap('Greys')
        edge_norm = matplotlib.colors.Normalize(vmin = 0, vmax = 4)

        for start_node, end_node, params in G.edges(data=True):
            edge_color = edge_cmap(edge_norm(round(params['weight'])))
            edge_params['line_width'][(start_node, end_node)] = params['weight']
            edge_params['color'][(start_node, end_node)] = RGB(int(255*list(edge_color)[0]), 
                                                               int(255*list(edge_color)[1]), 
                                                               int(255*list(edge_color)[2]))
            edge_params['weight'][(start_node, end_node)] = 1
            
        # take care of the nodes
        node_params = {'color': {}, 'label': {}}
        for node, _ in G.nodes(data=True):
            node_label = sorted(list(self.selected_families_summary.keys()))[node]
            node_params['label'][node] = node_label

            node_color = self.family_colors[node_label]['Color (RGBA)']
            node_params['color'][node] = node_color

        nx.set_node_attributes(G, node_params['color'], "node_color")
        nx.set_node_attributes(G, node_params['label'], "node_label")
        nx.set_edge_attributes(G, edge_params['color'], "edge_color")
        nx.set_edge_attributes(G, edge_params['line_width'], "line_width")
        nx.set_edge_attributes(G, edge_params['weight'], "weight")
        
        return G    
    
    def remove_non_adjacent_edges(self, **kwargs) -> nx.Graph:  
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        adjacency_matrix, family_labels = self.get_adjcency_matrix(**kwargs)        
        edges_to_remove = []
        for start_node, end_node, params in self.G.edges(data=True):
            start_node_index = family_labels.index(self.families_present[start_node])
            end_node_index = family_labels.index(self.families_present[end_node])
            
            if adjacency_matrix[start_node_index][end_node_index] == 0:
                edges_to_remove.append((start_node, end_node))
        
        for edge in edges_to_remove:
            self.G.remove_edge(*edge)
        
        return self.G 

    def get_adjcency_matrix(self, **kwargs) -> tuple[np.ndarray, list]: 
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        family_labels = sorted([family for family in self.families_summary.keys() 
                                if family > 0 and family < 10000])
        matrix = [[0 for family in family_labels] for family in family_labels]
        
        for operon in self.operons:
            for genomic_context in self.operons[operon]['operon_protein_families_structure']:
                for i in range(len(genomic_context)-1):
                    out_family = genomic_context[i]
                    in_family = genomic_context[i+1]

                    if all(family > 0 and family < 10000 for family in [out_family, in_family]):
                        matrix[family_labels.index(out_family)][family_labels.index(in_family)] += 1
                        matrix[family_labels.index(in_family)][family_labels.index(out_family)] += 1

        matrix = np.array(matrix)        
        return matrix, family_labels      
    
    def create_node_features(self, **kwargs) -> tuple[dict, list]:   
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        data = {'text_x': [],
                'text_y': [],
                'family': [],
                'tm_text_x': [],
                'tm_text_y': [],
                'tm_text': [],
                'tm_type': [],
                'tm_pred_text': [],
                'protein_name': [],
                'found_models': [],
                'model_links': []}
        
        for node in self.node_graph_coords:
            
            family = self.G.nodes[node]['node_label']
            coords = self.node_graph_coords[node]
            protein_name = self.families_summary[family]['name']
            
            if family > 0 and family < 10000 and 'function' in self.families_summary[family]:
                if 'TM_topology' in self.families_summary[family]['function']:
                    tm_type = self.families_summary[family]['function']["TM_topology"]
                    
                    if len(tm_type) > 0:
                        tm_text = 'TM'
                        tm_mode = 'Yes -> type:'
                    else:
                        tm_text = ''
                        tm_mode = 'No'
                else:
                    tm_type = ''
                    tm_text = ''
                    tm_mode = ''
            else:
                tm_type = 'n.a.'
                tm_text = ''
                tm_mode = 'n.a.'
            
            if 'model_state' in self.families_summary[family]:
                model_state = self.families_summary[family]['model_state']

                if model_state == 'Model exists':
                    model_state = 'Yes (click to view in Swiss-Model repository)'
                elif model_state == 'Model does not exist':
                    model_state = 'No (click to model with Swiss-Model)'
                else:
                    if family > 0 and family < 10000:
                        model_state = 'Not possible to find'
                    else:
                        model_state = ''
                
                structure = self.families_summary[family]['structure']
                if structure == '':
                    uniprot_code = self.families_summary[family]['uniprot_code']
                    structure = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
                
            else:
                model_state = 'n.a.'
                structure = 'n.a.'
            
            if family != 0 and family != self.reference_family and family < 10000:
                data['family'].append(family)
            else:
                data['family'].append(str(''))
                
            data['found_models'].append(model_state)
            data['model_links'].append(structure)
            
            y_range = (min([self.node_graph_coords[node][1] for node in self.node_graph_coords])-0.5, 
                       max([self.node_graph_coords[node][1] for node in self.node_graph_coords])+0.5)
            y_step = (y_range[1]-y_range[0])*0.09

            data['text_x'].append(coords[0])
            data['text_y'].append(coords[1]+y_step)
            
            data['tm_text_x'].append(coords[0])
            data['tm_text_y'].append(coords[1])
            data['tm_text'].append(tm_text)
            data['tm_pred_text'].append(tm_mode)
            data['tm_type'].append(tm_type)
            
            data['protein_name'].append(protein_name)
            
        
        tooltips = [('Protein name', "@protein_name"),
                    ('Protein family code', '@family'),
                    ("Predicted membrane protein", "@tm_pred_text @tm_type"),
                    ('Structural model found', '@found_models')] 
        
        return data, tooltips    
    
    def create_genomic_context_figure(self, **kwargs) -> figure:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        p_tooltips, p_data, p_yyticklabels = self.create_genomic_context_features(**kwargs)

        # the genomic_context figure
        p = figure(plot_width=self.most_common_gc_figure.plot_width, 
                   plot_height=self.syn_dendogram.height, 
                   x_range = self.most_common_gc_figure.x_range, 
                   y_range = self.syn_dendogram.y_range, 
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
    
    def create_genomic_context_features(self, **kwargs) -> tuple[list, dict, dict]:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)
        
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
        y_step = self.syn_den_data['y'][1] - self.syn_den_data['y'][0]
        y_half_height = y_step/4
        
        for i, current_target in enumerate(self.syn_den_data['leaf_label']):
            operon = [operon for operon in self.operons if current_target in self.operons[operon]['target_members']][0]
            current_assembly = self.all_syntenies[current_target]['assembly_id'][1]
            current_species = self.all_syntenies[current_target]['species']
            current_genomic_context_block = self.all_syntenies[current_target]['flanking_genes']
            current_species = self.all_syntenies[current_target]['species']
            current_reference_family = self.all_syntenies[current_target]['target_family']
            
            curr_y = self.syn_den_data['y'][i]

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
                    facecolor = self.family_colors[family]['Color (RGBA)']
                    edgecolor = self.family_colors[family]['Line color']
                    linestyle = self.family_colors[family]['Line style']
                else:
                    facecolor = self.family_colors[family]['Color (RGBA)']
                    edgecolor = self.family_colors[family]['Line color']
                    linestyle = self.family_colors[family]['Line style'] 
                    
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
                                    
                if family != 0 and family != self.reference_family and family < 10000:
                    data['family'].append(family)
                else:
                    data['family'].append(str(''))
                    
                data['species'].append(current_species)
                for level in ['superkingdom', 'phylum', 'class', 'order', 'genus']:
                    data[level].append(self.syn_den_data[level][i])
            
            if self.legend_mode in ['superkingdom', 'phylum', 'class', 'order', 'genus', 'species']:
                label = data[self.legend_mode][-1]
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
    
    def create_legend_figure(self, **kwargs) -> figure:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)        

        l_tooltips, l_data = self.create_legend_features(dx = 50, **kwargs)
        
        height = int(len(self.family_colors)*25*1.2)
        
        if len(l_data['ys']) > len(self.operons) or self.rescale_height:
            height = self.grid.height
            
        l = figure(plot_width=500, plot_height=height, x_range = [0,500], 
                   title = 'Protein families (click to view in Swiss-Model Repository)', 
                   toolbar_location="right")  # the legend figure
        
        for i, xs in enumerate(l_data['xs']):
            l.patch(xs, l_data['ys'][i], fill_color = l_data['facecolor'][i], 
                    line_color = l_data['edgecolor'][i], line_width = 1)	
        
        l.patches('xs', 'ys', fill_color = None, line_color = None, line_width = 0, source = l_data, 
                hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
                selection_fill_color='facecolor', selection_line_color='edgecolor',
                nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', 
                nonselection_fill_alpha=0.2)

        l.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", 
               text_font_size = {'value': '6pt'}, source = l_data)
        l.text('label_x', 'label_y', text = 'text', text_baseline="middle", text_align="left", 
               text_font_size = {'value': '8pt'}, source = l_data)
        l.text('tm_text_x', 'tm_text_y', text = 'tm_text',  text_color = "white", 
               text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, source = l_data)
        
        l.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
        l.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
        l.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
        l.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
        l.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
        l.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
        l.yaxis.axis_line_width = 0
        l.xaxis.axis_line_width = 0
        # define general features
        l.grid.visible = False
        l.outline_line_width = 0
        
        l.add_tools(HoverTool(tooltips=l_tooltips))
        l.add_tools(TapTool(callback = OpenURL(url='@model_links')))
                
        return l    
    
    def create_legend_features(self, **kwargs) -> tuple[list, dict]:
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)     

        data = {'xs': [],
                'ys': [],
                'edgecolor': [],
                'facecolor': [],
                'text': [],
                'label_x': [],
                'label_y': [],
                'text_x': [],
                'text_y': [],
                'tm_text_x': [],
                'tm_text_y': [],
                'tm_text': [],
                'tm_type': [],
                'tm_mode': [],
                'family': [],
                'found_models': [],
                'model_links': [],
                'keywords': [],
                'go_terms': [],
                'function': []}
        
        
        curr_y = len(self.family_colors.keys())
        
        for family in sorted(list(self.families_summary.keys())):
            if family == 0:
                data['text'].append('Non-conserved gene')
            elif family == self.reference_family:
                data['text'].append('Target protein: {}'.format(self.families_summary[family]['name']))
            elif family == 10000:
                data['text'].append('Pseudogene')
            elif family in self.families_summary:
                data['text'].append(self.families_summary[family]['name'])
            
            if family != 0 and family != self.reference_family and family < 10000:
                data['family'].append(family)
            else:
                data['family'].append(str(''))
            
            if 'model_state' in self.families_summary[family]:
                model_state = self.families_summary[family]['model_state']

                if model_state == 'Model exists':
                    model_state = self.families_summary[family]['structure']
                elif model_state == 'Model does not exist':
                    model_state = 'click to model/view with Swiss-Model'
                else:
                    if family > 0 and family < 10000:
                        model_state = 'Not possible to find'
                    else:
                        model_state = 'n.a.'
            else:
                model_state = ''
                
            data['found_models'].append(model_state)
            
            if family > 0 and family < 10000 and 'function' in self.families_summary[family]:
                if 'TM_topology' in self.families_summary[family]['function']:
                    tm_type = self.families_summary[family]['function']["TM_topology"]
                    keywords = ', '.join(sorted(self.families_summary[family]['function']['Keywords']))
                    go_terms = '; '.join(sorted(self.families_summary[family]['function']['GO_terms']))
                    function = self.families_summary[family]['function']['Function_description']

                    if len(tm_type) > 0:
                        tm_text = 'TM'
                        tm_mode = 'Yes -> type:'
                    else:
                        tm_text = ''
                        tm_mode = 'No'
                else:
                    tm_type = ''
                    tm_text = ''
                    tm_mode = ''
                    keywords = ''
                    go_terms = ''   
                    function = ''
            else:
                tm_type = 'n.a.'
                tm_text = ''
                tm_mode = 'n.a.'
                keywords = 'n.a.'
                go_terms = 'n.a.'   
                function = 'n.a.'
            
            if 'structure' in self.families_summary[family]:
                structure = self.families_summary[family]['structure']
                if structure == '':
                    uniprot_code = self.families_summary[family]['uniprot_code']
                    structure = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
            else:
                structure = 'n.a.'
                
            data['model_links'].append(structure)
            
            data['facecolor'].append(self.family_colors[family]['Color (RGBA)'])
            data['edgecolor'].append(self.family_colors[family]['Line color'])
            data['xs'].append([5, 5, self.dx-5, self.dx, self.dx-5])
            data['ys'].append([curr_y-0.25, curr_y+0.25, curr_y+0.25, curr_y, curr_y-0.25])
            data['text_x'].append(((self.dx-10)/2)+5)
            data['text_y'].append(curr_y+0.25)
            
            data['label_x'].append(self.dx + 10) 
            data['label_y'].append(curr_y) 
            
            data['tm_text_x'].append(((self.dx-10)/2)+5)
            data['tm_text_y'].append(curr_y)
            data['tm_text'].append(tm_text)
            data['tm_type'].append(tm_type)
            data['tm_mode'].append(tm_mode)
            
            data['go_terms'].append(go_terms)
            data['keywords'].append(keywords)
            data['function'].append(function)
            
            curr_y -= 1
            
        tooltips = [('Protein family', '@text'),
                    ('Protein family code', '@family'),
                    ('Structure', '@found_models'),
                    ('Predicted membrane protein', '@tm_mode @tm_type'),
                    ('Keywords', '@keywords'),
                    ('GO terms', '@go_terms'),
                    ('Function', '@function')]
        
        return tooltips, data    