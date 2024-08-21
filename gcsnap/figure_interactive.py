import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import statistics
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial import distance

from collections import Counter

import networkx as nx
from Bio import Phylo

# pip install bokeh
from bokeh.plotting import figure, output_file, gridplot, save
from bokeh.colors import RGB
from bokeh.models import HoverTool, TapTool, LassoSelectTool, Range1d, LinearAxis, WheelZoomTool, Circle, MultiLine, Panel, Tabs
from bokeh.models import ColumnDataSource, DataTable, DateFormatter, TableColumn, Legend, HTMLTemplateFormatter
from bokeh.models.callbacks import OpenURL
from bokeh.models.widgets import Div
from bokeh.models.graphs import from_networkx
from bokeh.layouts import row, column

from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.figure import Figure
from gcsnap.rich_console import RichConsole

import logging
logger = logging.getLogger('iteration')

class InteractiveFigure:
    """
    Methods to create the interactive genomic context figures.
    """
    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str, ref_family: str, 
                 family_colors: dict, starting_directory: str):
        """
        Initialize the InteractiveFigure object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_label (str): The label of the output.
            ref_family (str): The reference family.
            family_colors (dict): The dictionary with the family colors.
            starting_directory (str): The path to the starting directory.        """             
        
        # Extract parameters from config and gc
        kwargs = {k: v['value'] for k, v in config.arguments.items()}
        kwargs.update({
            'config': config,
            'out_label': out_label,
            'gc': gc,
            'operons': gc.get_selected_operons(),
            'most_populated_operon': gc.get_most_populated_operon(),
            'syntenies': gc.get_syntenies(),
            'families_summary': gc.get_families(),
            'taxonomy': gc.get_taxonomy(),
            'input_targets': gc.get_curr_targets(),            
            'reference_family': ref_family,
            'family_colors': family_colors,
            'starting_directory': starting_directory            
        })

        # Set all attributes
        self._set_attributes(**kwargs)

        # change sort_mode and input targets if needed
        if kwargs['in_tree'] is not None:
            kwargs['sort_mode'] = 'tree'    
            kwargs['input_targets'] = [target for operon in self.operons 
                                       for target in self.operons[operon]['target_members']]                

        self.console = RichConsole()

    def _set_attributes(self, **kwargs):
        """
        Set attributes from kwargs.
        """        
        for key, value in kwargs.items():
            # Only set the attribute if it does not already exist
            if not hasattr(self, key):
                setattr(self, key, value)   

    def run(self) -> None:
        """
        Run the creation of the interactive genomic context figures
            - Create the most common genomic context figure.
            - Create the gene co-occurrence network figure.
            - Create the gene adjacency network figure.
            - Create the dendogram for the genomic context block.
            - Create the genomic context figure.
            - Create the legend figure.
        """        
        with self.console.status('Creating interactive genomic context figures'):
            # Make a copy of the current object's attributes (all contained in instance __dict__)
            kwargs = self.__dict__.copy()

            # Work on most conserved genomic context figure        
            most_common_gc_figure = Figure.create_most_common_genomic_features_figure(**kwargs) 
            # Work on gene co-occurence figure    
            try:    
                coocurrence_figure, graph_coord = self.create_graph_figure(
                                                    most_common_gc_figure = most_common_gc_figure, 
                                                    graph_coord = {},
                                                    mode = 'coocurrence',
                                                    previous_net = '',
                                                    **kwargs) 
                adjacency_figure, graph_coord = self.create_graph_figure(mode = 'adjacency', 
                                                            graph_coord=graph_coord, 
                                                            previous_net=coocurrence_figure, **kwargs)
            except ValueError as e:
                self.console.print_warning('No co-occurrence detected returning empty figure.')
                message = 'No co-occurrence detected'
                coocurrence_figure = self.create_empty_graph_figure(**kwargs, 
                                                                    mode = 'coocurrence', fig_message = message)
                adjacency_figure = self.create_empty_graph_figure(**kwargs, 
                                                                  mode = 'adjacency',fig_message = message)
            except Exception as e:
                self.console.print_warning('Error creating gene co-occurrence or adjacency network figure.')
                logger.warning(e) 
                message = 'Error creating gene co-occurrence or adjacency network figure'
                coocurrence_figure = self.create_empty_graph_figure(**kwargs, 
                                                                    mode = 'coocurrence', fig_message = message)
                adjacency_figure = self.create_empty_graph_figure(**kwargs, 
                                                                  mode = 'adjacency',fig_message = message)                


            # Work on dendogram for the genomic context block
            syn_dendogram, syn_den_data = Figure.make_dendogram_figure(show_leafs = False,                                                                        
                                                                height_factor = 25*1.2, 
                                                                distance_matrix = None, 
                                                                labels = None, 
                                                                colors = None,
                                                                **kwargs)
            
            # height_factor already present in kwargs as added in previous method
            # same true for most_common_gc_figure
            genomic_context_figure = Figure.create_genomic_context_figure(syn_dendogram = syn_dendogram,
                                                                most_common_gc_figure = most_common_gc_figure,      
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

            # output to static HTML file
            f_name = '{}_interactive_output.html'.format(self.out_label)
            output_file(os.path.join(os.getcwd(),f_name))
            save(grid)
        self.console.print_info('Genomic context visualization created in {}'.format(f_name))

    def create_empty_graph_figure(self, **kwargs) -> figure:
        """
        Create an empty gene co-occurrence or adjacency network figure.

        Returns:
            tuple[figure,dict]: The Bokeh figure and the dictionary with the node coordinates.
        """        
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        if self.mode == 'coocurrence':
            title = 'Gene co-occurrence network'
        elif self.mode == 'adjacency':
            title = 'Gene adjcency network'

        g = figure(plot_width = self.most_common_gc_figure.plot_width, 
                   plot_height = self.most_common_gc_figure.plot_height, 
                   title = title)
        g.text(x=[0.5], y=[0.5], text=[self.fig_message], text_align='center', text_baseline='middle', 
                text_font_size='12pt')
        g.xaxis.visible = False
        g.yaxis.visible = False
        g.grid.visible = False        
        
        return g

    def create_graph_figure(self, **kwargs) -> tuple[figure, dict]:
        """
        Create the gene co-occurrence or adjacency network figure.

        Returns:
            tuple[nx.Graph, dict]: The networkx graph and the dictionary with the node coordinates.
        """        
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        matrix, selected_families_summary = self.get_coocurrence_matrix(**kwargs)
        graph = self.get_graph_from_matrix(matrix = matrix, 
                                           selected_families_summary = selected_families_summary, 
                                           **kwargs)
        
        if self.mode == 'coocurrence':
            title = 'Gene co-occurrence network'
        elif self.mode == 'adjacency':
            graph = self.remove_non_adjacent_edges(graph=graph, 
                                    families_present = sorted(list(selected_families_summary.keys())),
                                    **kwargs)
            title = 'Gene adjcency network'
        
        if len(self.graph_coord) == 0:
            self.graph_coord = nx.spring_layout(graph)         
        
        node_data, node_tooltips = self.create_node_features(node_graph_coords=self.graph_coord, 
                                                             graph=graph, **kwargs)
            
        if self.previous_net != '':
            x_range = self.previous_net.x_range
            y_range = self.previous_net.y_range
            
        else:
            x_range = (min([self.graph_coord[node][0] for node in self.graph_coord])-0.5, 
                       max([self.graph_coord[node][0] for node in self.graph_coord])+0.5)
            y_range = (min([self.graph_coord[node][1] for node in self.graph_coord])-0.5, 
                       max([self.graph_coord[node][1] for node in self.graph_coord])+0.5)

        g = figure(plot_width = self.most_common_gc_figure.plot_width, 
                   plot_height = self.most_common_gc_figure.plot_height, 
                   x_range = x_range, y_range = y_range, title = title)

        graph_renderer = from_networkx(graph, self.graph_coord, scale=1, center=(0, 0))
        graph_renderer.edge_renderer.glyph = MultiLine(line_width="line_width", line_color = "edge_color")
        graph_renderer.node_renderer.glyph = Circle(size=22, fill_color = "node_color")
        
        g.renderers.append(graph_renderer)
        
        g.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", 
               text_font_size = {'value': '6pt'}, source = node_data)
        g.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", text_baseline="middle", 
               text_align="center", text_font_size = {'value': '6pt'}, source = node_data)
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
        
        return g, self.graph_coord        
    
    def get_coocurrence_matrix(self, **kwargs) -> tuple[np.ndarray, dict]:
        """
        Get the co-occurrence matrix.

        Returns:
            tuple[np.ndarray, dict]: The co-occurrence matrix and the dictionary with the selected families summary.
        """        
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        family_labels = sorted([family for family in self.families_summary.keys() 
                                if family > 0])
        matrix = [[0 for family in family_labels] for family in family_labels]

        context_count = 0
        for operon in self.operons:
            for genomic_context in self.operons[operon]['operon_protein_families_structure']:
                for i in range(len(set(genomic_context))):
                    for j in range(len(set(genomic_context))):
                        if i > j:
                            out_family = list(set(genomic_context))[i]
                            in_family = list(set(genomic_context))[j]

                            if all(family > 0 for family in [out_family, in_family]):
                                matrix[family_labels.index(out_family)][family_labels.index(in_family)] += 1
                                matrix[family_labels.index(in_family)][family_labels.index(out_family)] += 1


        matrix = np.array(matrix)
        # normalize the matrix if its not empty
        if np.amax(matrix) == 0:
            raise ValueError('Maximum value in the co-occurrence matrix is zero')
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
        """
        Get the graph from the matrix.

        Returns:
            nx.Graph: The networkx graph.
        """        
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        G = nx.from_numpy_array(self.matrix)

        # take care of the edges
        edge_params = {'color': {}, 'weight': {}, 'line_width': {}}
        edge_cmap = plt.get_cmap('Greys')
        edge_norm = mcolors.Normalize(vmin = 0, vmax = 4)

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
        """
        Remove non-adjacent edges from the graph.

        Returns:
            nx.Graph: The networkx graph.
        """         
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
        """
        Get the adjacency matrix.

        Returns:
            tuple[np.ndarray, list]: The adjacency matrix and the list with the family labels.
        """        
        # set attributes from kwargs, it updates any keywords that are passed to method
        self._set_attributes(**kwargs)

        family_labels = sorted([family for family in self.families_summary.keys() 
                                if family > 0])
        matrix = [[0 for family in family_labels] for family in family_labels]
        
        for operon in self.operons:
            for genomic_context in self.operons[operon]['operon_protein_families_structure']:
                for i in range(len(genomic_context)-1):
                    out_family = genomic_context[i]
                    in_family = genomic_context[i+1]

                    if all(family > 0 for family in [out_family, in_family]):
                        matrix[family_labels.index(out_family)][family_labels.index(in_family)] += 1
                        matrix[family_labels.index(in_family)][family_labels.index(out_family)] += 1

        matrix = np.array(matrix)        
        return matrix, family_labels      
    
    def create_node_features(self, **kwargs) -> tuple[dict, list]:  
        """
        Create the node features.

        Returns:
            tuple[dict, list]: The dictionary with the node features and the list with the tooltips.
        """         
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
            family = self.graph.nodes[node]['node_label']
            coords = self.node_graph_coords[node]
            protein_name = self.families_summary[family]['name']
            
            if family > 0 and 'function' in self.families_summary[family]:
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
                    if family > 0:
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
            
            if family > 0 and family != self.reference_family:
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
        
    def create_legend_figure(self, **kwargs) -> figure:
        """
        Create the legend figure.

        Returns:
            figure: The Bokeh figure.
        """        
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
        """
        Create the legend features.

        Returns:
            tuple[list, dict]: The list with the tooltips and the dictionary with the legend features.
        """        
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
            elif family == -1:
                data['text'].append('Pseudogene')
            elif family in self.families_summary:
                data['text'].append(self.families_summary[family]['name'])
            
            if family > 0 and family != self.reference_family:
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
                    if family > 0:
                        model_state = 'Not possible to find'
                    else:
                        model_state = 'n.a.'
            else:
                model_state = ''
                
            data['found_models'].append(model_state)
            
            if family > 0 and 'function' in self.families_summary[family]:
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