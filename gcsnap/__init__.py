"""
GCsnap Package

This package contains the main classes and functions to run the GCsnap pipeline.
The following modules are included:

Modules used in the pipeline:
    - rich_console: RichConsole class to print messages in color and style.
    - configuration: Configuration class to parse the configuration file and arguments.
    - timing: Timing and Timer class to measure the time of each step.
    - targets: Target class to parse the target list.
    - sequence_mapping: SequenceMapping class to map sequences to the genome.
    - assemblies: Assemblies class to finde, download and parse the assemblies.
    - genomic_context: GenomicContext class to store the genomic context.
    - sequences: Sequences class to find and extract sequences.
    - families: Families class to find and assign families to genes.
    - families_functions_structures: FamiliesFunctionsStructures class to find and
        assign functions and structures to families.
    - operons: Operons class to find operons to genes.
    - taxonomy: Taxonomy class to find and add taxonomy information.
    - tm_segments: TMsegments class to find and add transmembrane segments.
    - figures: Figures class to create figures.

Helper modules than are usable standalone as well:
    - utils: Helper functions to split lists and dictionaries and parallelize functions.
    - entrez_query: EntrezQuery class to query the NCBI Entrez database.
    - mmseqs_cluster: MMseqsCluster class to cluster sequences with MMseqs2.
    - uniprot_api: UniprotAPI class to query the Uniprot API.
    - uniprot_dbs_dict: UniprotDBsDict class to store the Uniprot databases.
    - apis: SwissProtAPI, AlphaFoldAPI and UniProtAPI classes to query functional annotation.
    - assembly_links: AssemblyLinks class to find and store assembly links.

"""

# # Import classes and functions from various modules of the pipeline in __main__.py
# from .rich_console import RichConsole 
# from .configuration import Configuration 
# from .timing import Timing
# from .targets import Target 
# from .sequence_mapping import SequenceMapping
# from .assemblies import Assemblies
# from .genomic_context import GenomicContext
# from .sequences import Sequences
# from .families import Families
# from .families_functions_structures import FamiliesFunctionsStructures
# from .operons import Operons
# from .taxonomy import Taxonomy
# from .tm_segments import TMsegments
# from .figures import Figures

# # import helper classes and function
# from .utils import processpool_wrapper, split_dict_chunks, split_list_chunks
# from .entrez_query import EntrezQuery
# from .mmseqs_cluster import MMseqsCluster
# from .uniprot_api import submit_id_mapping
# from .uniprot_api import check_id_mapping_results_ready
# from .uniprot_api import get_id_mapping_results_link
# from .uniprot_api import get_id_mapping_results_search
# from .uniprot_dbs_dict import UniprotDict
# from .apis import SwissProtAPI, AlphaFoldAPI, EbiAPI
# from .assembly_links import AssemblyLinks


# # using from GCsnap import *, only the names included in __all__ will be imported.
# __all__ = [
#     'RichConsole', 'Configuration', 'Timing', 'Target', 'SequenceMapping', 'Assemblies',
#     'GenomicContext', 'Sequences', 'Families', 'FamiliesFunctionsStructures',
#     'Operons', 'Taxonomy', 'TMsegments', 'Figures', 'processpool_wrapper', 'split_dict_chunks',
#     'split_list_chunks', 'EntrezQuery', 'MMseqsCluster', 'UniprotDict',
#     'SwissProtAPI', 'AlphaFoldAPI', 'EbiAPI', 'AssemblyLinks', 'submit_id_mapping',
#     'check_id_mapping_results_ready', 'get_id_mapping_results_link', 'get_id_mapping_results_search'
# ]