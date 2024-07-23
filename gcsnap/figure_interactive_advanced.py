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
from bokeh.models.widgets import Div
from bokeh.layouts import row, column

from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.figure import Figure
from gcsnap.mmseqs_cluster import MMseqsCluster
from gcsnap.rich_console import RichConsole

class AdvancedInteractiveFigure:
    """ 
    Methods to create advanced interactive figures.
    """
    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str, 
                 ref_family: str, family_colors: dict, starting_directory: str):
        """
        Initialize the AdvancedInteractiveFigure object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_label (str): The label of the output.
            ref_family (str): The reference family.
            family_colors (dict): The dictionary with the colors of the families.
            starting_directory (str): The path to the starting directory.
        """
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
            'reference_family': ref_family,
            'family_colors': family_colors,
            'starting_directory': starting_directory
        })

        # store configuration arguments needed with different values at certain points
        # to avoid multiple keyword arguments in the methods
        self.input_targets_start =  gc.get_curr_targets()
        self.sort_mode_start = kwargs['sort_mode']
        self.gc_legend_mode_start = kwargs['gc_legend_mode']
        # remove from kwargs
        kwargs.pop('sort_mode')
        kwargs.pop('gc_legend_mode')

        # handling of mode in case in_tree argument was set
        if kwargs['in_tree'] is not None:
            self.sort_mode_start = 'tree'   
            self.input_targets_start = [target for operon in self.operons 
                                       for target in self.operons[operon]['target_members']]                 

        # Set all attributes
        self._set_attributes(**kwargs)

        self.console = RichConsole()

    def _set_attributes(self, **kwargs) -> None:
        """
        Helper method to set attributes from kwargs.
        """        
        for key, value in kwargs.items():
            # Only set the attribute if it does not already exist
            if not hasattr(self, key):
                setattr(self, key, value)   

    def run(self) -> None:
        """
        Run the creation of the advanced interactive figures.
            - Summary page
            - Per operon type page
        """        
        # Make a copy of the current object's attributes (all contained in instance __dict__)
        kwargs = self.__dict__.copy()

        self.console.print_working_on('Making advanced interactive figure summary page.')  
        f_name = self.create_interactive_output_html(**kwargs)
        self.console.print_info('Summary visualization created in {}'.format(f_name))       

        with self.console.status('Making advanced interactive figure per operon type page.'):           
            self.create_advanced_operon_interactive_output(**kwargs)
        f_name = '{}_advanced_operons_interactive_output_??.html'.format(self.out_label)                    
        self.console.print_info('Operon visualization created in {}'.format(f_name))    

    def create_interactive_output_html(self, **kwargs) -> str:
        """
        Create the interactive output html file.

        Returns:
            str: The name of the output file.
        """        
        self._set_attributes(**kwargs)

        # Create the dataframe with all data used to plot. it allows for cross interativity between plots
        cluster_colors = self.define_operons_colors(clusters = self.operons.keys(), 
                                                    mode = 'bokeh', cmap = 'gist_rainbow')
        scatter_data = self.create_data_structure(cluster_colors = cluster_colors,
                                                **kwargs)
        
        # Plot the scatter of the operons PaCMAP coordinates
        operons_scatter = self.create_operons_clusters_scatter(scatter_data = scatter_data, **kwargs)
        # Create network of operon types similarities (based on the centroids in PaCMAP space and the minimum distance between cluster members)
        operons_network, operons_distance_matrix = self.create_avg_operons_clusters_network(
                            operons_scatter = operons_scatter, scatter_data = scatter_data, 
                            **kwargs)
        
        # Plot the CLANS map of the input target sequences if given
        clans_scatter = self.create_clans_map_scatter(scatter_data = scatter_data, **kwargs)
        scatter_row = gridplot([[operons_scatter, operons_network, clans_scatter]], merge_tools = True)

        # Now create two tabs
        all_tabs = []

        # 1. The tab of species (i.e. a table where each line corresponds 
        # to an input target and includes the taxonomy information as well as the cluster type it belongs to)
        targets_table = self.create_targets_summary_table(**kwargs)
        all_tabs.append(Panel(child = targets_table, title = 'Input targets summary'))

        # 2. The tab the shows the gene composition of each cluster type, 
        # sorted by similarity and connected by a cladogram build from the 
        # distance matrix used to build the network above

        # Make the dendogram
        # sort mode for this call is changed in kwargs
        oprn_dendogram, oprn_den_data = Figure.make_dendogram_figure(show_leafs = False, 
                                                input_targets = None,                                                                        
                                                height_factor = 25*1.2, 
                                                sort_mode = 'operon clusters', 
                                                distance_matrix = (1-operons_distance_matrix), 
                                                labels = sorted([i for i in self.operons if '-' not in i]), 
                                                colors = cluster_colors,
                                                **kwargs)
    
        family_freq_figure = self.create_family_frequency_per_operon_figure(oprn_dendogram = oprn_dendogram, 
                                                                    oprn_den_data = oprn_den_data, 
                                                                    height_factor = 25*1.2, 
                                                                    min_freq = self.min_family_freq_accross_contexts/100,
                                                                    **kwargs)
        
        gc_row = gridplot([[oprn_dendogram, family_freq_figure]], merge_tools = True)
        all_tabs.append(Panel(child = gc_row, title = 'Genomic contexts clusters hierarchy'))

        all_tabs = Tabs(tabs=all_tabs)
        grid = gridplot([[scatter_row],[all_tabs]], merge_tools = False)

        # output to static HTML file
        f_name = '{}_advanced_operons_output_summary.html'.format(self.out_label)
        output_file(os.path.join(os.getcwd(),f_name))
        save(grid)           
        return f_name
    
    def create_advanced_operon_interactive_output(self, **kwargs) -> None:
        """
        Create the advanced operon interactive output.
        """        
        self._set_attributes(**kwargs)

        for operon_type in sorted(self.operons.keys()):
            curr_operon = {operon_type: self.operons[operon_type]}

            if len(curr_operon[operon_type]['target_members']) > 1:                
                all_tabs = []

                div = Div(text="""<b>Detailed view of individual genomic context types:</b></br></br>
                            Below you find multiple <b>tabs</b> coresponding to each individual
                            cluster you see on the scatter plot above.</br>
                            <b>Click</b> on the tab to have a detailed view of the different genomic
                            context clusters, clustered based on the similarity of their family 
                            composition.</br></br>
                            The depiction is interactive, <b>hover</b> and <b>click</b> to get more information!</br></br>  """) 

                # Work on most conserved genomic context figure
                # remove operons to not have multiple values in keyqord arguments
                operons = self.operons
                kwargs.pop('operons', None)
                most_common_gc_figure = Figure.create_most_common_genomic_features_figure(
                                            operons = curr_operon, **kwargs) 
                kwargs['operons'] = operons
                # Work on dendogram for the genomic context block
                syn_dendogram, syn_den_data = Figure.make_dendogram_figure(show_leafs = False, 
                                                                    input_targets = self.input_targets_start,
                                                                    sort_mode = self.sort_mode_start,                                                                       
                                                                    height_factor = 25*1.2, 
                                                                    distance_matrix = None, 
                                                                    labels = None, 
                                                                    colors = None,
                                                                    **kwargs)
                # Work on the genomic context block
                genomic_context_figure = Figure.create_genomic_context_figure(
                                                most_common_gc_figure = most_common_gc_figure,
                                                syn_dendogram = syn_dendogram,
                                                syn_den_data = syn_den_data,
                                                gc_legend_mode = 'species',
                                                **kwargs)
                
                # Make the table of family frequencies
                family_table, table_div = self.create_families_frequency_table(**kwargs)

                gc_row = row(syn_dendogram, genomic_context_figure)

                grid = gridplot([[gc_row],[table_div],[family_table]], merge_tools = True)
                curr_tab = Panel(child = grid, title = operon_type)
                all_tabs.append(curr_tab)

                all_tabs = Tabs(tabs=all_tabs)
                grid = gridplot([[div], [all_tabs]], merge_tools = True)

                f_name = '{}_advanced_operons_interactive_output_{}.html'.format(self.out_label, operon_type)
                output_file(os.path.join(os.getcwd(),f_name))
                save(grid)

    def define_operons_colors(self, **kwargs) -> dict:
        """
        Define the colors of the operons.

        Returns:
            dict: The dictionary with the colors of the operons.
        """        
        self._set_attributes(**kwargs)

        colors = {}
        cmap = plt.get_cmap(self.cmap)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(self.clusters))
        colours = [cmap(norm(i)) for i in range(len(self.clusters))]
        #random.shuffle(colours)

        for i, label in enumerate(sorted(self.clusters)):
            if label not in colors:
                colors[label] = {}

            if '-' in label:   # singleton clusters
                colors[label]['Color (RGBA)'] = 'grey'
                colors[label]['Color (tuplet)'] = 'grey'
                colors[label]['Line color'] = 'black'
                colors[label]['Line style'] = '-'
                colors[label]['Size'] = '2'

            else:  
                if self.mode == 'matplotlib':
                    colors[label]['Color (RGBA)'] = [int(255*j) for j in colours[i]]
                    colors[label]['Color (tuplet)'] = colours[i]
                    colors[label]['Line color'] = 'black'
                    colors[label]['Line style'] = '-'
                elif self.mode == 'bokeh':
                    colors[label]['Color (RGBA)'] = RGB(int(255*list(colours[i])[0]), 
                                                        int(255*list(colours[i])[1]), 
                                                        int(255*list(colours[i])[2]))
                    colors[label]['Color (tuplet)'] = colours[i]
                    colors[label]['Line color'] = 'black'
                    colors[label]['Line style'] = '-'
                colors[label]['Size'] = '5'

        return colors         
    
    def create_data_structure(self, **kwargs) -> tuple:
        """
        Create the data structure for the scatter plots.

        Returns:
            tuple: The tooltips and the ColumnDataSource.
        """        
        self._set_attributes(**kwargs)

        if self.clans_file is not None:
            clans_coords = self.parse_coordinates_from_clans(**kwargs)
            seq_in_clans_without_operon = list(set(clans_coords.keys()).difference(list(self.syntenies.keys())))
        else:
            clans_coords = self.generate_coordinates_for_clans(**kwargs)
            seq_in_clans_without_operon = []
        
        data = {'x': [],		  # pacmap x coord for individual operon
                'y': [],		  # pacmap y coord for individual operon
                'avg_x': [],	  # average x coord of the corresponding operon type
                'avg_y': [],	  # average y coord of the corresponding operon type
                'clans_x': [],	# clans x coord for the input target
                'clans_y': [],	# clans y coord for the input target
                'edgecolor': [],
                'facecolor': [],
                'size': [],
                'node_size': [],
                'text': [],
                'type': [],
                'species': [],
                'target': []}
        
        # first go through the sequences in the clans map that were not assigned any gene cluster type 
        # (eg, the clans map is for a much larger set of sequences, our set of sequences is a subset of the clans map
        # or sequences were filteresd out because they were in a partial genome context)
        for member in seq_in_clans_without_operon:
            clans_x, clans_y = clans_coords[member]['xy']
            text = 'Not analysed'
            operon_type = 'n.a.'
            species	 = clans_coords[member]['species']
            facecolor   = 'grey'
            edgecolor   = 'grey'
            size		= 1
        
            data['clans_x'].append(clans_x)
            data['clans_y'].append(-clans_y)
            data['x'].append(np.nan)
            data['y'].append(np.nan)
            data['avg_x'].append(np.nan)
            data['avg_y'].append(np.nan)
            data['text'].append(text)
            data['edgecolor'].append(edgecolor)
            data['facecolor'].append(facecolor)
            data['size'].append(size)
            data['node_size'].append(np.nan)
            data['type'].append(operon_type)
            data['species'].append(species)
            data['target'].append(member)
        
        # now go through the sequences in the gene clusters and collect their xy clans coordinates if they
        # are in the clans map
        for operon_type in sorted(list(self.operons.keys())):
            facecolor = self.cluster_colors[operon_type]['Color (RGBA)']
            edgecolor = self.cluster_colors[operon_type]['Line color']
            size = self.cluster_colors[operon_type]['Size']
            
            for i, member in enumerate(self.operons[operon_type]['target_members']):
                if member in clans_coords:
                    clans_x, clans_y = clans_coords[member]['xy']
                    clans_y = -clans_y
                else:
                    clans_x, clans_y = np.nan, np.nan
                
                x, y = self.operons[operon_type]['operon_filtered_PaCMAP'][i]
                text = self.operons[operon_type]['operon_protein_families_structure'][i]

                data['clans_x'].append(clans_x)
                data['clans_y'].append(clans_y)
                data['x'].append(x)
                data['y'].append(y)
                data['avg_x'].append(np.nan)
                data['avg_y'].append(np.nan)
                data['text'].append(text)
                data['edgecolor'].append(edgecolor)
                data['facecolor'].append(facecolor)
                data['size'].append(size)
                data['node_size'].append(np.nan)
                data['type'].append(operon_type.split()[2])
                data['species'].append(self.syntenies[member]['species'])
                data['target'].append(member)
            
            # add the data for the average position of the GC in the pacmap space if it is not a -00001 type
            if '-' not in operon_type:
                avg_x, avg_y = self.operons[operon_type]['operon_centroid_PaCMAP']
                
                data['clans_x'].append(np.nan)
                data['clans_y'].append(np.nan)
                data['x'].append(avg_x)
                data['y'].append(avg_y)
                data['avg_x'].append(avg_x)
                data['avg_y'].append(avg_y)
                data['text'].append(np.nan)
                data['edgecolor'].append(edgecolor)
                data['facecolor'].append(facecolor)
                data['size'].append('1')
                data['node_size'].append('10')
                data['type'].append(operon_type.split()[2])
                data['species'].append(np.nan)
                data['target'].append(np.nan)
                    
        tooltips = [('Operon/GC type', '@type'),
                    ('EntrezID', '@target'),
                    ('Species', '@species'),]
        
        return tooltips, ColumnDataSource(data)  

    def parse_coordinates_from_clans(self, **kwargs) -> dict:
        """
        Parse the coordinates from the CLANS file.

        Returns:
            dict: The dictionary with the coordinates.
        """        
        self._set_attributes(**kwargs)
        
        clans_coords = {}
        seq_map = {}
        
        found_seq_block = False
        found_coords = False
        seq_count = 0
        with open(self.clans_file, 'r') as inclans:
            for line in inclans:
                if '<seq>' in line:
                    found_seq_block = True
                elif '</seq>' in line:
                    found_seq_block = False
                elif found_seq_block and line.startswith('>'):
                    line = line[1:]
                    ncbi_code = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()
                    seq_map[seq_count] = ncbi_code
                    
                    if '[' in line and ']' in line:
                        species = line.split('[')[1].split(']')[0]
                    else:
                        species = 'n.a.'
                    
                    clans_coords[ncbi_code]={'species': species, 'xy': None}                    
                    seq_count += 1
                
                elif '<pos>' in line:
                    found_coords = True
                elif '</pos>' in line:
                    found_coords = False
                elif not found_seq_block and found_coords:
                    coords = [float(i) for i in line.strip().split(' ')]
                    ncbi_code = seq_map[int(coords[0])]
                    clans_coords[ncbi_code]['xy'] = coords[1:3]
                
        return clans_coords  
    
    def generate_coordinates_for_clans(self, **kwargs) -> dict:
        """
        Generate the coordinates for the CLANS map.

        Returns:
            dict: The dictionary with the coordinates.
        """        
        self._set_attributes(**kwargs)

        in_syntenies = {}
        for target in self.syntenies:
            in_syntenies[target] = {'flanking_genes': {}}

            assembly_targetid = self.syntenies[target]['assembly_id'][0]
            context_idx = self.syntenies[target]['flanking_genes']['ncbi_codes'].index(assembly_targetid)

            for key in self.syntenies[target]['flanking_genes']:
                if type(self.syntenies[target]['flanking_genes'][key]) == list:
                    if key == 'ncbi_codes':
                        in_syntenies[target]['flanking_genes'][key] = [target]
                    else:
                        in_syntenies[target]['flanking_genes'][key] = [self.syntenies[target]
                                                                       ['flanking_genes'][key][context_idx]]
        
        # Here MMSeqs only with targets and not flanking genes
        # hence we use a copy of the gc object and replaye the syntenies with the in_syntenies
        gc_temp = self.gc.copy()
        gc_temp.syntenise = in_syntenies
        out_label = '{}_targets'.format(self.out_label)
        out_dir = os.path.join(os.getcwd(), out_label)
        mmseqs = MMseqsCluster(self.config, gc_temp, out_dir)
        mmseqs.run()
        distance_matrix = mmseqs.get_distance_matrix()
        ordered_ncbi_codes = mmseqs.get_cluster_order()

        paCMAP_embedding = pacmap.PaCMAP(n_components = 2)
        paCMAP_coordinat = paCMAP_embedding.fit_transform(distance_matrix)

        clans_coords = {ordered_ncbi_codes[i]: {'xy': paCMAP_coordinat[i]} 
                        for i in range(len(paCMAP_coordinat))}

        return clans_coords    
    
    def create_operons_clusters_scatter(self, **kwargs) -> figure:
        """
        Create the scatter plot of the operons clusters.

        Returns:
            figure: The Bokeh figure.
        """        
        self._set_attributes(**kwargs)

        p_tooltips, p_data = self.scatter_data

        p = figure(title = 'Genomic context types/clusters', width=500, height=500)
        p.add_layout(Legend(orientation="horizontal"), 'above')

        p.circle('x', 'y', size='size', line_color='edgecolor', fill_color='facecolor', 
                 alpha=1, source = p_data)

        p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
        p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
        p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
        p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
        p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
        p.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
        p.yaxis.axis_line_width = 0
        p.xaxis.axis_line_width = 0
        # define general features
        p.grid.visible = False
    #	 p.outline_line_width = 0

        p.xaxis.axis_label = ""
        p.yaxis.axis_label = ""

        p.add_tools(HoverTool(tooltips=p_tooltips))
        p.add_tools(LassoSelectTool())

        p.background_fill_color = "lightgrey"
        p.background_fill_alpha = 0.2

        p.legend.click_policy="hide"

        return p
    
    def create_avg_operons_clusters_network(self, **kwargs) -> tuple:
        """
        Create the network of the operons clusters.

        Returns:
            tuple: The Bokeh figure and the distance matrix.
        """        
        self._set_attributes(**kwargs)

        p_tooltips, p_data = self.scatter_data

        similarity_matrix, operons_labels = self.get_avgoperons_distance_matrix(**kwargs)
        similarity_matrix = self.normalize_matrix(similarity_matrix = similarity_matrix, 
                                                  power = 30, **kwargs)

        edge_tooltips, edge_data = self.create_edges_data(operons_labels = operons_labels,
                                                          similarity_matrix = similarity_matrix,
                                                          **kwargs)

        p = figure(title = 'Genomic context types/clusters similarity network', 
                   plot_width=self.operons_scatter.width,
                   plot_height=self.operons_scatter.height, 
                   x_range = self.operons_scatter.x_range, 
                   y_range = self.operons_scatter.y_range)
        p.add_layout(Legend(orientation="horizontal"), 'above')

        p.multi_line('x', 'y', color='color', alpha='alpha', source=edge_data, name='edges')
        p.circle('avg_x', 'avg_y', size='node_size', line_color='edgecolor', 
                 fill_color='facecolor', alpha=1, source = p_data, name='nodes')

        p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
        p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
        p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
        p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
        p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
        p.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
        p.yaxis.axis_line_width = 0
        p.xaxis.axis_line_width = 0
        # define general features
        p.grid.visible = False
    #	 p.outline_line_width = 0

        p.xaxis.axis_label = ""
        p.yaxis.axis_label = ""

        p.add_tools(HoverTool(tooltips=p_tooltips, names=['nodes']))
        p.add_tools(HoverTool(tooltips=edge_tooltips, names=['edges']))
        p.add_tools(LassoSelectTool())

        p.background_fill_color = "lightgrey"
        p.background_fill_alpha = 0.2

        p.legend.click_policy="hide"

        return p, similarity_matrix    
    
    def get_avgoperons_distance_matrix(self,**kwargs) -> tuple:
        """
        Get the average operons distance matrix.

        Returns:
            tuple: The distance matrix and the operons labels.
        """        
        self._set_attributes(**kwargs)
            
        matrix = [[0 for i in self.operons if '-' not in i] for i in self.operons if '-' not in i]        
        selected_operons = [i for i in self.operons if '-' not in i]
        selected_operons = sorted(selected_operons)
        for i, i_operon_type in enumerate(selected_operons):
            for j, j_operon_type in enumerate(selected_operons):
                if i > j:
                    dists = []
                    for i_member_pacmap in self.operons[i_operon_type]['operon_filtered_PaCMAP']:
                        for j_member_pacmap in self.operons[j_operon_type]['operon_filtered_PaCMAP']:
                            dist = np.linalg.norm(np.array(i_member_pacmap)-np.array(j_member_pacmap))
                            dists.append(dist)
                    dist = min(dists)
                    matrix[i][j] = dist
                    matrix[j][i] = dist
        
        return np.array(matrix), selected_operons    
    
    def normalize_matrix(self, **kwargs) -> np.array:
        """
        Normalize the matrix.

        Returns:
            np.array: The normalized matrix.
        """        
        self._set_attributes(**kwargs)

        min_dist = pd.DataFrame(self.similarity_matrix).min().min()
        max_dist = pd.DataFrame(self.similarity_matrix).max().max()
        
        normalised_matrix = 1-(np.array(self.similarity_matrix)-min_dist)/max_dist
        normalised_matrix = np.power(normalised_matrix, self.power)
        
        return normalised_matrix    
    
    def create_edges_data(self, **kwargs) -> tuple:
        """
        Create the edges data.

        Returns:
            tuple: The tooltips and the ColumnDataSource.
        """        
        self._set_attributes(**kwargs)

        data = {'x': [],
                'y': [],
                'color': [],
                'alpha':[]}

        for i, i_operon_type in enumerate(self.operons_labels):
            for j, j_operon_type in enumerate(self.operons_labels):
                if i > j:
                    x_start, y_start = self.operons[i_operon_type]['operon_centroid_PaCMAP']
                    x_end, y_end = self.operons[j_operon_type]['operon_centroid_PaCMAP']

                    data['x'].append([x_start, x_end])
                    data['y'].append([y_start, y_end])
                    data['color'].append('black')
                    data['alpha'].append(round(self.similarity_matrix[i][j], 1))

        tooltips = [('Relative distance/alpha', '@alpha')]

        return tooltips, data    
    
    def create_clans_map_scatter(self, **kwargs) -> figure:
        """
        Create the CLANS map scatter plot.

        Returns:
            figure: The Bokeh figure.
        """        
        self._set_attributes(**kwargs)

        p_tooltips, p_data = self.scatter_data

        p = figure(title = 'Sequence similarity cluster (CLANS) map', width=500, height=500)
        p.add_layout(Legend(orientation="horizontal"), 'above')

        p.circle('clans_x', 'clans_y', size='size', line_color='edgecolor', 
                 fill_color='facecolor', alpha=1, source = p_data)

        p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
        p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
        p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
        p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
        p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
        p.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
        p.yaxis.axis_line_width = 0
        p.xaxis.axis_line_width = 0
        # define general features
        p.grid.visible = False
    #	 p.outline_line_width = 0

        p.xaxis.axis_label = ""
        p.yaxis.axis_label = ""

        p.add_tools(HoverTool(tooltips=p_tooltips))
        p.add_tools(LassoSelectTool())

        p.background_fill_color = "lightgrey"
        p.background_fill_alpha = 0.2

        return p    
    
    def create_targets_summary_table(self, **kwargs) -> DataTable:
        """
        Create the targets summary table.

        Returns:
            DataTable: The Bokeh DataTable.
        """        
        self._set_attributes(**kwargs)

        t_data, t_columns = self.create_targets_table_data(**kwargs)	
        t = DataTable(source=ColumnDataSource(t_data), columns=t_columns, width=1500, height=500)
        return t
    
    def create_targets_table_data(self, **kwargs) -> tuple:
        """
        Create the targets table data.

        Returns:
            tuple: The data and the columns.
        """        
        self._set_attributes(**kwargs)

        data = {'Target EntrezID': [],
                'Superkingdom': [],
                'Phylum': [],
                'Class': [],
                'Order': [],
                'Genus': [],
                'Species': [],
                'Genomic context type': [],
                'color': []}

        for superkingdom in self.taxonomy.keys():
            for phylum in self.taxonomy[superkingdom].keys():
                for taxclass in self.taxonomy[superkingdom][phylum].keys():
                    for order in self.taxonomy[superkingdom][phylum][taxclass].keys():
                        for genus in self.taxonomy[superkingdom][phylum][taxclass][order].keys():
                            for species in self.taxonomy[superkingdom][phylum][taxclass][order][genus].keys():
                                for target in self.taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
                                    operon_type = self.syntenies[target]['operon_type']
                                    operon_type = 'GC Type {:05d}'.format(operon_type)

                                    if '-' not in operon_type:
                                        operon_color = self.cluster_colors[operon_type]['Color (RGBA)']

                                        data['Target EntrezID'].append(target)
                                        data['Superkingdom'].append(superkingdom)
                                        data['Phylum'].append(phylum)
                                        data['Class'].append(taxclass)
                                        data['Order'].append(order)
                                        data['Genus'].append(genus)
                                        data['Species'].append(species)
                                        data['Genomic context type'].append(operon_type.split()[-1])
                                        data['color'].append(operon_color)

        # columns = [TableColumn(field=i, title=i) for i in data.keys()]

        columns = [TableColumn(field=i, title=i) if i not in ['color']
            else TableColumn(field=i, title='Genomic context color', 
                             formatter=HTMLTemplateFormatter(template='<span style="color:<%= value %>;font-size:18pt;text-shadow: 1px 1px 2px #000000;">&#9632;</span>'))
            for i in data.keys()]																			 

        return data, columns    
    
    def create_family_frequency_per_operon_figure(self, **kwargs) -> figure:
        """
        Create the family frequency per operon figure.

        Returns:
            figure: The Bokeh figure.
        """        
        self._set_attributes(**kwargs)

        families_frequencies = self.compute_families_frequencies_per_operon_cluster(**kwargs)
        families_frequencies_matrix = self.clean_uncommon_families(families_frequencies = families_frequencies,
                                                                   operons_labels = self.oprn_den_data['leaf_label'], 
                                                                   **kwargs)

        p_tooltips, p_data, p_yyticklabels = self.create_family_spectrum_data(
                                families_frequencies_matrix = families_frequencies_matrix,
                                dx = 1,
                                **kwargs)
            
        p = figure(plot_width=1250, plot_height=self.oprn_dendogram.height, 
                   x_range = [0, len(families_frequencies_matrix.columns)], 
                   y_range = self.oprn_dendogram.y_range, 
                   toolbar_location="left", 
                   title = 'Genomic context family spectrum (hover to get more information and click to model with SWISS-MODEL)') 
        
        p.patches('xs', 'ys', fill_color = 'facecolor', line_color = 'edgecolor', line_width = 1, source = p_data, 
                hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
                selection_fill_color='facecolor', selection_line_color='edgecolor',
                nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', nonselection_fill_alpha=0.2)
            
        p.yaxis.ticker = list(p_yyticklabels.keys())
        p.yaxis.major_label_overrides = {int(i): p_yyticklabels[i] for i in p_yyticklabels.keys()}
        
        p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
        p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
        p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
        p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
        p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
        p.yaxis.axis_line_width = 0
        p.xaxis.axis_line_width = 0
        # define general features
        p.grid.visible = False
        p.outline_line_width = 0

        p.add_tools(HoverTool(tooltips=p_tooltips))
        p.add_tools(TapTool(callback = OpenURL(url='@model_links')))

        return p    
    
    def compute_families_frequencies_per_operon_cluster(self, **kwargs) -> dict:
        """
        Compute the families frequencies per operon cluster.

        Returns:
            dict: The families frequencies per operon cluster.
        """        
        self._set_attributes(**kwargs)

        families_frequencies = {}
        for operon_type in self.operons:
            curr_frequencies, curr_number_of_operons = self.compute_families_frequencies(
                oper = {operon_type: self.operons[operon_type]})
            families_frequencies[operon_type] = curr_frequencies

            for family in families_frequencies[operon_type]:
                families_frequencies[operon_type][family] = families_frequencies[operon_type][family]/curr_number_of_operons

        return families_frequencies    
    
    def compute_families_frequencies(self, **kwargs) -> tuple:
        """
        Compute the families frequencies.

        Returns:
            tuple: The families frequencies and the number of operons.
        """        
        self._set_attributes(**kwargs)

        families_frequencies = {}
        number_of_operons = 0
        for operon_type in self.oper:
            for i, target in enumerate(self.oper[operon_type]['target_members']):
                number_of_operons += 1
                for family in self.oper[operon_type]['operon_protein_families_structure'][i]:
                    if family not in families_frequencies:
                        families_frequencies[family] = 1
                    else:
                        families_frequencies[family] += 1

        return families_frequencies, number_of_operons    

    def clean_uncommon_families(self, **kwargs) -> pd.DataFrame:
        """
        Clean the uncommon families.

        Returns:
            pd.DataFrame: The cleaned uncommon families.
        """        
        self._set_attributes(**kwargs)

        families_labels = [i for i in sorted(self.families_summary.keys()) if i < 0]
        matrix = [[0 for family in families_labels] for operon_type in self.operons_labels]

        for i, operon_type in enumerate(self.operons_labels):
            for j, family in enumerate(families_labels):
                if family in self.families_frequencies[operon_type]:
                    matrix[i][j] = self.families_frequencies[operon_type][family]

        matrix = pd.DataFrame(matrix, index=self.operons_labels, columns=families_labels).T
        matrix = matrix.loc[matrix.max(axis=1) > self.min_freq]
        matrix['sum'] = matrix.sum(axis=1)
        matrix = matrix.sort_values(by='sum', ascending=False)
        matrix = matrix.drop('sum', axis=1).T    

        return matrix

    def create_family_spectrum_data(self, **kwargs) -> tuple:
        """
        Create the family spectrum data.

        Returns:
            tuple: The tooltips, the data, and the y tick labels.
        """        
        self._set_attributes(**kwargs)

        data = {'xs': [],
                'ys': [],
                'edgecolor': [],
                'facecolor': [],
                'text': [],
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

        yyticklabels = self.oprn_den_data['leaf_label']
        yys = []
        y_step = self.oprn_den_data['y'][1] - self.oprn_den_data['y'][0]

        for i, operon_type in enumerate(self.families_frequencies_matrix.index):
            curr_y = self.oprn_den_data['y'][i]-(y_step/2)
            
            for j, family in enumerate(self.families_frequencies_matrix.columns):
                dy = self.families_frequencies_matrix.at[operon_type,family]
                if dy > 1:
                    dy = 1
                    
                dy = dy*y_step
                curr_x = j
                
                if family == self.reference_family:
                    data['text'].append('Target protein: {}'.format(self.families_summary[family]['name']))
                elif family in self.families_summary:
                    data['text'].append(self.families_summary[family]['name'])

                data['family'].append(family)

                if 'model_state' in self.families_summary[family]:
                    model_state = self.families_summary[family]['model_state']

                    if model_state == 'Model exists':
                        model_state = self.families_summary[family]['structure']
                    elif model_state == 'Model does not exist':
                        model_state = 'click to model/view with Swiss-Model'
                    else:
                        model_state = 'Not possible to find'
                else:
                    model_state = ''

                data['found_models'].append(model_state)

                if 'function' in self.families_summary[family]:
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
                data['xs'].append([curr_x, curr_x, curr_x+self.dx,curr_x+self.dx])
                data['ys'].append([curr_y, curr_y+dy, curr_y+dy, curr_y])
                data['text_x'].append(((curr_x+self.dx)/2)+curr_x)
                data['text_y'].append(curr_y+0.5)

                data['tm_text_x'].append(((curr_x+self.dx)/2)+curr_x)
                data['tm_text_y'].append(curr_y+0.25)
                data['tm_text'].append(tm_text)
                data['tm_type'].append(tm_type)
                data['tm_mode'].append(tm_mode)

                data['go_terms'].append(go_terms)
                data['keywords'].append(keywords)
                data['function'].append(function)

            curr_y -= 1
            
            yys.append(curr_y+(y_step/2))

        yyticklabels = {yys[i]: yyticklabels[i] for i in range(len(yyticklabels))}

        tooltips = [('Protein family', '@text'),
                    ('Protein family code', '@family'),
                    ('Structure', '@found_models'),
                    ('Predicted membrane protein', '@tm_mode @tm_type'),
                    ('Keywords', '@keywords'),
                    ('GO terms', '@go_terms'),
                    ('Function', '@function')]

        return tooltips, data, yyticklabels     

    def create_families_frequency_table(self, **kwargs) -> tuple:
        """
        Create the families frequency table.

        Returns:
            tuple: The Bokeh DataTable and the Div.
        """        
        self._set_attributes(**kwargs)
        
        t_data, t_columns = self.create_family_table_data(**kwargs)        
        t = DataTable(source=ColumnDataSource(t_data), columns=t_columns, width=2250, height=300)
        div = Div(text="""<b>Represented families:</b>""") 
        
        return t, div       

    def create_family_table_data(self, **kwargs) -> tuple:
        """
        Create the family table data.

        Returns:
            tuple: The data and the columns.
        """        
        self._set_attributes(**kwargs)
            
        data = {'Family': [],
                'Name': [],
                'UniprotKB': [],
                'Frequency': [],
                'Keywords': [],
                'GO terms: Molecular function (F)': [],
                'GO terms: Cellular component (C)': [],
                'GO terms: Biological process (P)': [],
                'Function': [],
                'Predicted membrane protein': [],
                'Structure': []}
        
        families_frequencies, number_of_operons = self.compute_families_frequencies(**kwargs)
        
        for family in families_frequencies:  
            if family > 0:
                data['Family'].append(family)
                data['Frequency'].append(round(families_frequencies[family]*100/number_of_operons, 1))
                data['Name'].append(self.families_summary[family]['name'])
                uniprot_code = 'nan'
                
                if 'function' in self.families_summary[family]:
                    uniprot_code = self.families_summary[family]['uniprot_code']
                    if 'TM_topology' in self.families_summary[family]['function']:
                        tm_type = self.families_summary[family]['function']["TM_topology"]
                        keywords = ', '.join(sorted(self.families_summary[family]['function']['Keywords']))
                        go_terms = sorted(self.families_summary[family]['function']['GO_terms'])
                        function = self.families_summary[family]['function']['Function_description']

                        if len(tm_type) > 0:
                            tm_mode = 'Yes'
                        else:
                            tm_mode = 'No'
                            
                        f_go = '; '.join([go for go in go_terms if go.startswith('F:')])
                        c_go = '; '.join([go for go in go_terms if go.startswith('C:')])
                        p_go = '; '.join([go for go in go_terms if go.startswith('P:')])
                
                    else:
                        tm_mode = ''
                        keywords = ''
                        f_go = '' 
                        c_go = '' 
                        p_go = '' 
                        function = ''
                else:
                    tm_mode = 'n.a.'
                    keywords = 'n.a.'
                    f_go = 'n.a.' 
                    c_go = 'n.a.' 
                    p_go = 'n.a.'   
                    function = 'n.a.'
                
                if 'structure' in self.families_summary[family]:
                    structure = self.families_summary[family]['structure']
                    if structure == '':
                        structure = 'Click to model'
                else:
                    structure = 'n.a.'
                
                data['UniprotKB'].append(uniprot_code)
                data['Keywords'].append(keywords)
                data['GO terms: Molecular function (F)'].append(f_go)
                data['GO terms: Cellular component (C)'].append(c_go)
                data['GO terms: Biological process (P)'].append(p_go)
                data['Function'].append(function)
                data['Predicted membrane protein'].append(tm_mode)
                data['Structure'].append(structure)
        
        columns = [TableColumn(field=i, title=i) if i not in ['Structure', 'UniprotKB']
                    else TableColumn(field=i, title=i, formatter=HTMLTemplateFormatter(
                        template='<a href="https://swissmodel.expasy.org/repository/uniprot/<%= UniprotKB %>"><%= value %></a>')) if i=='Structure'
                    else TableColumn(field=i, title='Representative UniprotKB', 
                             formatter=HTMLTemplateFormatter(
                                 template='<a href="https://www.uniprot.org/uniprot/<%= value %>"><%= value %></a>'))  
                    for i in data.keys()]
        
        df = pd.DataFrame(data).sort_values(by='Frequency', ascending=False)
        # ColumnDataSource takes a dictionary
        data = df.to_dict(orient='list')
        
        return data, columns
    
    def compute_families_frequencies(self, **kwargs) -> tuple:
        """
        Compute the families frequencies.

        Returns:
            tuple: The families frequencies and the number of operons.
        """        
        self._set_attributes(**kwargs)

        families_frequencies = {}
        number_of_operons = 0
        for operon_type in self.operons:
            for i, target in enumerate(self.operons[operon_type]['target_members']):
                number_of_operons += 1
                for family in self.operons[operon_type]['operon_protein_families_structure'][i]:
                    if family not in families_frequencies:
                        families_frequencies[family] = 1
                    else:
                        families_frequencies[family] += 1

        return families_frequencies, number_of_operons    