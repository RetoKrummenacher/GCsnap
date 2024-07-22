

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.entrez_query import EntrezQuery

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import split_dict_chunks

class Taxonomy:
    def __init__(self, config: Configuration, gc: GenomicContext):
        self.config = config
        self.cores = config.arguments['n_cpu']['value']

        if config.arguments['get_taxonomy']['value']:
            self.mode = 'taxonomy'
            self.msg = 'Find and map taxonomy'
        else:
            self.mode = 'as_input'
            self.msg = 'Create input taxonomy dictionary. No taxomomies are searched'

        # set parameters
        self.gc = gc

        self.console = RichConsole()

    def get_taxonomy(self) -> dict:
        return self.taxonomy

    def run(self) -> None:
        self.cleaned_taxonomy = {}
        if self.mode == 'taxonomy':
            # find all taxonomies via entrez
            self.find_taxonomies(self.gc.get_all_taxids())

        # here we parallellize over chunks, so as many chunks as 
        # there are cores
        parallel_args = split_dict_chunks(self.gc.get_syntenies(), self.cores)

        with self.console.status(self.msg):
            dict_list = processpool_wrapper(self.cores, parallel_args, self.run_each)
            # combine results
            # as this is a heavily nested dictionary, we need some recursive functionality
            self.taxonomy = self.merge_all_dicts(dict_list)

    def run_each(self, arg: dict) -> dict:
        content_dict = arg
        taxonomy = {}
        for target in content_dict.keys():
            ncbi_code = content_dict[target]['assembly_id'][0]
            taxID = content_dict[target]['flanking_genes']['taxID']
            species = content_dict[target]['flanking_genes']['species']
  
            genus = 'na'
            if self.mode == 'taxonomy':
                if species.split()[0] == 'uncultured':
                    genus = species.split()[1]
                else:
                    genus = species.split()[0]

            # take what is available from entrez taxonomies requests
            # if mode is as_input, we will have 'na' for all values as 
            # self.clean_taxonomy is empty
            superkingdom = self.clean_taxonomy.get(taxID,{}).get('superkingdom','na')
            phylum = self.clean_taxonomy.get(taxID,{}).get('phylum','na')
            taxclass = self.clean_taxonomy.get(taxID,{}).get('class','na')
            order = self.clean_taxonomy.get(taxID,{}).get('order','na')

            if self.mode == 'taxonomy':
                # now all are either correct string or 'na'
                if phylum == 'na':
                    phylum = '{}_na'.format(superkingdom)
                if taxclass == 'na':
                    taxclass = '{}_na'.format(phylum)
                if order == 'na':
                    order = '{}_na'.format(taxclass)

            # create the nested dictionaries
            taxonomy.setdefault(superkingdom, {}).setdefault(phylum, {}).setdefault(taxclass, {}) \
                .setdefault(order, {}).setdefault(genus, {}).setdefault(
                    species, {'ncbi_codes': [], 'target_members': []}
                )            

            taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members'].append(target)
            taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes'].append(ncbi_code)            

        return taxonomy

    def find_taxonomies(self, taxids: list) -> None:
        # get the information for all ncbi codes
        entrez = EntrezQuery(self.config, taxids, db='taxonomy', 
                             retmode='xml', logging=True)
        self.entrez_taxonomy = entrez.run()
        # remove all entries that are None
        # we do this here to have keep entrez result as is (to realy see what was not there)
        self.clean_taxonomy = {
            target: {key: value for key, value in sub_dict.items() if value is not None}
            for target, sub_dict in self.entrez_taxonomy.items()
        }

    def merge_nested_dicts(self, dict1: dict, dict2: dict) -> dict: 
        for key, value in dict2.items():
            if key in dict1:
                if isinstance(value, dict) and isinstance(dict1[key], dict):
                    self.merge_nested_dicts(dict1[key], value)
                elif key == "ncbi_codes" or key == "target_members":
                    if dict1[key] == value:
                        dict1[key] = list(set(dict1[key] + value))
                    else:
                        dict1[key].extend(value)
                        dict1[key] = list(set(dict1[key]))
                else:
                    dict1[key] = value
            else:
                dict1[key] = value
        return dict1

    def merge_all_dicts(self, dict_list: list[dict]) -> dict:
        merged_dict = dict_list[0]
        for d in dict_list[1:]:
            merged_dict = self.merge_nested_dicts(merged_dict, d)
        return merged_dict

