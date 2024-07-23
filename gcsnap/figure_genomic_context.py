import matplotlib.pyplot as plt
import numpy as np

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext

import logging
logger = logging.getLogger(__name__) # inherits configuration from main logger

class GenomicContextFigure:
    """ 
    Methods to create genomic context figures.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cmap (str): The colormap to use.
        out_format (str): The output format of the figures.
        out_label (str): The label of the output.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        operons (dict): The dictionary with the operons.
        most_populated_operon (str): The most populated operon.
        syntenies (dict): The dictionary with the syntenies.
        families_summary (dict): The dictionary with the families summary.
        reference_family (str): The reference family.
        family_colors (dict): The dictionary with the family colors.
        console (RichConsole): The RichConsole object to print messages.
    """
    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str, ref_family: str, family_colors: dict):
        """
        Initialize the GenomicContextFigure object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_label (str): The label of the output.
            ref_family (str): The reference family.
            family_colors (dict): The dictionary with the family colors.
        """        
        self.config = config
        self.cmap = config.arguments['genomic_context_cmap']['value']
        self.out_format = config.arguments['out_format']['value']

        # set parameters
        self.out_label = out_label
        self.gc = gc
        self.operons = gc.get_selected_operons()
        self.most_populated_operon = gc.get_most_populated_operon()
        self.syntenies = gc.get_syntenies()
        self.families_summary = gc.get_families()
        self.reference_family = ref_family
        self.family_colors = family_colors

        self.console = RichConsole()

    def run(self) -> None:
        """
        Run the creation of genomic context figures.
        """        
        try:
            with self.console.status('Creating genomic context figures'):                
                self.draw_genomic_context()
                self.draw_genomic_context_legend()
            f_name = '{}_genomic_context.{}'.format(self.out_label, self.out_format)
            self.console.print_info('Genomic context figures created in {}'.format(f_name))
        except Exception as e:
            self.console.print_warning('Images not created due to error.')
            logger.warning('Error message: {}'.format(e))

    def draw_genomic_context(self) -> None:
        """
        Draw the genomic context figure.
        """        
        curr_y_level = len(self.operons.keys())
        all_xs = []
        all_populations = []
        yticklabels = []

        plt.clf()
        if len(self.operons) == 1:
            fig, ax = plt.subplots(1, 2, figsize=(20, 2), gridspec_kw={'width_ratios': [4, 1]})
        elif len(self.operons) < 5:
            fig, ax = plt.subplots(1, 2, figsize=(20, len(self.operons)), gridspec_kw={'width_ratios': [4, 1]})
        else:
            fig, ax = plt.subplots(1, 2, figsize=(20, int(len(self.operons)/1.5)), gridspec_kw={'width_ratios': [4, 1]})

        for operon in sorted(list(self.operons.keys())):
            current_target = self.operons[operon]['target_members'][0]
            current_genomic_context_block = self.syntenies[current_target]['flanking_genes']
            current_species = self.syntenies[current_target]['species']

            operon_population = len(self.operons[operon]['target_members'])*100/len(self.syntenies)

            for i, flanking_gene in enumerate(current_genomic_context_block['ncbi_codes']):
                family = current_genomic_context_block['families'][i]
                gene_dx = current_genomic_context_block['relative_ends'][i] - \
                            current_genomic_context_block['relative_starts'][i]+1
                gene_direction = current_genomic_context_block['directions'][i]

                if gene_direction == '-':
                    gene_x_tail = current_genomic_context_block['relative_ends'][i]
                    gene_dx = gene_dx*(-1)
                else:
                    gene_x_tail = current_genomic_context_block['relative_starts'][i]

                # make genomic_context side			
                if family == 0:
                    zorder = family
                    facecolor = self.family_colors[family]['Color (tuplet)']
                    edgecolor = self.family_colors[family]['Line color']
                    linestyle = self.family_colors[family]['Line style']
                else:
                    zorder = len(current_genomic_context_block['ncbi_codes']) - i + 1
                    facecolor = self.family_colors[family]['Color (tuplet)']
                    edgecolor = self.family_colors[family]['Line color']
                    linestyle = self.family_colors[family]['Line style']

                ax[0].arrow(gene_x_tail, curr_y_level, gene_dx, 0, width=0.5, head_width=0.5, 
                            length_includes_head = True, head_length = 100, zorder = zorder, 
                            facecolor = facecolor, edgecolor = edgecolor, linestyle = linestyle)
                    
                text_x = gene_x_tail + (gene_dx/2)
                if family < 0 and family != self.reference_family:
                    ax[0].text(text_x, curr_y_level+0.3, str(family), horizontalalignment='center')

                # make histogram side
                ax[1].arrow(0, curr_y_level, operon_population, 0, width=0.5, head_width=0.5,
                             length_includes_head = True, head_length = 0, facecolor = 'black', 
                             edgecolor = 'black')
                ax[1].text(operon_population+2.5, curr_y_level, '{}'.format(
                    round(operon_population, 1)), verticalalignment = 'center')
        
                all_xs.append(current_genomic_context_block['relative_ends'][i])
                all_xs.append(current_genomic_context_block['relative_starts'][i])
            
            all_populations.append(operon_population)
            curr_y_level -= 1
            yticklabels.append('{}: {} ({})'.format(operon, current_target, current_species))
            
        yticklabels.append('')
        yticklabels.reverse()

        # format genomic_context side
        ax[0].set_xlim(min(all_xs)-100, max(all_xs)+100)
        ax[0].set_ylim(0, len(self.operons.keys())+1)

        # TODO: Here 
        #   File "C:\MT\GCsnap\gcsnap\figure_genomic_context.py", line 115, in draw_genomic_context
        #     ax[0].set_yticklabels(yticklabels, fontsize = 10, horizontalalignment='left')
        #   File "C:\Users\retok\miniconda3\envs\gcsnap\Lib\site-packages\matplotlib\axes\_base.py", line 74, in wrapper
        #     return get_method(self)(*args, **kwargs)
        #            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #   File "C:\Users\retok\miniconda3\envs\gcsnap\Lib\site-packages\matplotlib\axis.py", line 2060, in set_ticklabels
        #     raise ValueError(
        # ValueError: The number of FixedLocator locations (3), usually from a call to set_ticks, does not match the number of labels (2).        

        ax[0].set_yticks(np.arange(0, len(yticklabels), 1.0))
        ax[0].set_yticklabels(yticklabels, fontsize = 10, horizontalalignment='left')
        ax[0].spines['right'].set_visible(False)
        ax[0].spines['left'].set_visible(False)
        yax = ax[0].get_yaxis()
        pad = max(len(label) for label in yticklabels)*6
        yax.set_tick_params(pad=pad)
        
        ax[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

        # format histogram side
        ax[1].set_xlim(0, max(all_populations)+5)
        ax[1].set_ylim(0, len(self.operons.keys())+1)
        ax[1].spines['right'].set_visible(False)
        ax[1].tick_params(axis='y', which='both', left = False, right = False, labelleft=False, labelright=False)
        ax[1].set_xlabel('Operon frequency (%)')
        ax[1].set_title('{}% total operons represented'.format(round(sum(all_populations), 1)))

        try:
            plt.tight_layout()
        except:
            pass
        plt.savefig('{}_genomic_context.{}'.format(self.out_label, self.out_format), format = self.out_format)
        plt.close('all')        

    def draw_genomic_context_legend(self) -> None:
        """
        Draw the genomic context legend.
        """        
        curr_y_level = len(self.family_colors.keys())
        x_tail = 0
        dx = 5

        plt.clf()
        fig, ax = plt.subplots(figsize=(10, int(len(self.family_colors)/1.5)))
        for family in sorted(list(self.family_colors.keys())):
            plt.arrow(x_tail, curr_y_level, dx, 0, width=0.5, head_width=0.5, 
                      length_includes_head = True, head_length = 0.5, 
                      facecolor = self.family_colors[family]['Color (tuplet)'], 
                      edgecolor = self.family_colors[family]['Line color'], 
                      linestyle = self.family_colors[family]['Line style'])
            if family == 0:
                plt.text(dx + 2, curr_y_level, 'Non-conserved gene')
            elif family == self.reference_family:
                plt.text(dx + 2, curr_y_level, 'Target protein: {}'.format(
                    self.families_summary[family]['name']) )
            elif family == -1:
                plt.text(dx + 2, curr_y_level, 'Pseudogene')
            elif family in self.families_summary:
                plt.text(2.25, curr_y_level+0.3, str(family), horizontalalignment='center')
                plt.text(dx + 2, curr_y_level, self.families_summary[family]['name'])

            curr_y_level -= 1

        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.yticks(fontsize = 10)
        plt.tick_params(axis='y', which='both', left = False, right = False, 
                        labelleft=False, labelright=False)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        
        plt.xlim(0, 50)
        plt.ylim(0, len(self.family_colors.keys())+1)
        
        plt.savefig('{}_genomic_context_legend.{}'.format(self.out_label, self.out_format), 
                    format = self.out_format)
        plt.close('all')