import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import statistics
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial import distance

from collections import Counter

from Bio import Phylo

from bokeh.plotting import figure
from bokeh.colors import RGB
from bokeh.models import HoverTool, TapTool
from bokeh.models.callbacks import OpenURL

from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.figure_genomic_context import GenomicContextFigure
from gcsnap.figure_interactive import InteractiveFigure
from gcsnap.figure_interactive_advanced import AdvancedInteractiveFigure

class Figures:
    """ 
    Methods and attributes to create the genomic context figures.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        operon_cluster_advanced (bool): The boolean to use the advanced interactive figure.
        cmap (str): The name of the colormap.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        operons (dict): The dictionary with the operons.
        most_populated_operon (str): The most populated operon.
        syntenies (dict): The dictionary with the syntenies.
        families_summary (dict): The dictionary with the families.
        out_label (str): The label of the output.
        starting_directory (str): The path to the starting directory.
        reference_family (str): The reference family
    """

    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str, starting_directory: str):
        """
        Initialize the Figures object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_label (str): The label of the output.
            starting_directory (str): The path to the starting directory.
        """        
        self.config = config
        self.operon_cluster_advanced = config.arguments['operon_cluster_advanced']['value']
        self.cmap = config.arguments['genomic_context_cmap']['value']

        # set parameters
        self.gc = gc
        self.operons = gc.get_selected_operons()
        self.most_populated_operon = gc.get_most_populated_operon()
        self.syntenies = gc.get_syntenies()
        self.families_summary = gc.get_families()
        self.out_label = out_label
        self.starting_directory = starting_directory

        # define the reference family as the one 
        # of the target in the most populated operon type
        self.reference_family = self.syntenies[self.operons[self.most_populated_operon]
                                                ['target_members'][0]]['target_family']
        
    def run(self) -> None:
        """
        Run the creation of the genomic context figures:
            - Create genomic context figure with GenomicContextFigure.
            - Create interactive figure with InteractiveFigure or AdvancedInteractiveFigure.
        """        
        # 1. Create genomic context figures
        family_colors = self.define_family_colors(mode = 'matplotlib', cmap = self.cmap)
        gcf = GenomicContextFigure(self.config, self.gc, self.out_label, 
                                   self.reference_family, family_colors)  
        gcf.run()

        # 2. Interactive HTML figures
        if self.config.arguments['interactive']['value']:
            family_colors = self.define_family_colors(mode = 'bokeh', cmap = self.cmap)          
            if self.operon_cluster_advanced:
                aigcf = AdvancedInteractiveFigure(self.config, self.gc, self.out_label,
                                                self.reference_family, family_colors,
                                                self.starting_directory)
                aigcf.run()
            else:         
                igcf = InteractiveFigure(self.config, self.gc, self.out_label,
                                        self.reference_family, family_colors,
                                        self.starting_directory)
                igcf.run()

    def define_family_colors(self, mode: str,  cmap: str) -> dict:
        """
        Define the colors for the families.

        Args:
            mode (str): The mode to define the colors (matplotlib or bokeh).
            cmap (str): The name of the colormap.

        Returns:
            dict: The dictionary with the colors for the families.
        """        
        families = list(self.families_summary.keys())
        colors = {}

        cmap = plt.get_cmap(cmap)
        norm = mcolors.Normalize(vmin=0, vmax=len(families))

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
            elif type(label) == int and label == -1:  # a pseudogene
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