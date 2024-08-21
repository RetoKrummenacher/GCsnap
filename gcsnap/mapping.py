import pandas as pd
import numpy as np
import os
import glob
# pip install pypdl
from pypdl import Pypdl
import gzip
import shutil

import dask.dataframe as dd
from dask.distributed import LocalCluster, Client
from dask_jobqueue import SLURMCluster

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.utils import process_wrapper
from gcsnap.uniprot_dbs_dict import UniprotDict
from gcsnap.db_handler_uniprot_mappings import UniprotMappingsDBHandler


import logging
logger = logging.getLogger(__name__) # inherits configuration from main logger

class SequenceMapping:
    """ 
    Methods and attributes to extract the target sequences for the overall raw mapping file.
    The resulting DataFrame is the base to map targets to other types and also
    retrieve the Taxonomy ID and the PDB ID.
    The raw mapping file si donwloaded in advance and read with Dask.
    The raw files and further information can be found here:
    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

    The idmapping_selected.tab raw file is stored in parquet partitions for further processing.
    This is done in advance with the create_parquet_partitions method.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        target_list (list): The list of target sequences.
        cores (int): The number of CPU cores to use.
        console (RichConsole): The RichConsole object to print messages.
        uniprot_dict (dict): The dictionary with the identifiers standard information and the target sequences.
        mapping_df (pd.DataFrame): The DataFrame with the mapping results.
        block_size (int): The block size of the partition files to create.
        mapping_file (str): The path to the raw mapping file including the file name.
        parquet_path (str): The path to the parquet partition files.
    """

    def __init__(self, config: Configuration, path_to_raw_mapping: str, target_list: list = None,
                 block_size: int = 128e6):
        """
        Initialize the SequenceMapping object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            path_to_raw_mapping (str): The path to the raw mapping file including the file name.
            target_list (list, optional): The list of target sequences. Defaults to None
            block_size (int, optional): The block size of the partition files to create. Defaults to 128e6.
        """        
        self.target_list = target_list
        self.config = config

        # file with the data
        self.mapping_file = path_to_raw_mapping
        self.mapping_path = os.path.dirname(path_to_raw_mapping)
        self.parquet_path = os.path.join(os.path.dirname(path_to_raw_mapping),
                                         'idmapping.parquet.' + str(round(block_size/10**6)))  
        # block size of the partition files to create
        self.block_size = block_size

        # extract needed values from config
        self.cores = config.arguments['n_cpu']['value']
        self.console = RichConsole()

    def run(self):
        """
        Run the mapping of target sequences to the specified id standard:
            - Create the UniProt dictionary with the target sequences.
            - Read the parquet partitions with Dask.
            - Finalize the mapping.
        Uses Dask DataFrames within a Dask Cluster for parallel processing.
        """        
        with self.console.status('Mapping'):
            self.create_uniprot_dbs_with_targets()
   
            # Generate filtered data (not yet computed)
            ddf = self.generate_data_filter()

            # compute the filtered data
            cluster = SLURMcontrol(self.cores).start_cluster()
            with Client(cluster) as client:
                mapping_df = ddf.compute()
            cluster.close()

            # make them unique
            self.mapping_df = self.make_unique(mapping_df, self.uniprot_dict.get_target_types())
            self.finalize()

    def download_and_decrompress_mapping_file(self) -> None:
        """
        Download the mapping file from UniProt in parallel with pypdl Downloader from
        https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz.

        """        
        url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'
        target_file = os.path.join(self.mapping_path, os.path.basename(url))

        # Download multithreaded with pypdl Downloader
        dl = Pypdl()
        dl.start(url=url,
                file_path=target_file, 
                segments=self.cores, 
                multisegment=True, 
                block=False, 
                display=False)
        
        # retreive total size of the file does not work, the file does not have this information
        # total_size = dl.size
        # print(f"Total size: {total_size}")

        # set it manually (lower than the actual size, but it is enough for the progress bar)
        total_size = 11 * 10**9

        # Download multithreaded with pypdl Downloader
        with self.console.progress('Downloading mapping raw file {}'.format(os.path.basename(url)), 
                                   total=total_size) as (progress, task_id):
             
            while not dl.completed:
                current_size = dl.current_size
                progress.update(task_id, completed=current_size)

        dl.shutdown()       

        # Decompress the file
        with self.console.status('Decompressing mapping file to {}'.format(os.path.basename(self.mapping_file))):
            with gzip.open(target_file, 'rb') as f_in:
                with open(self.mapping_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)   

    def create_parquet_partitions(self) -> None:
        """
        Create the parquet partitions of the raw mapping file.
        This is done in advance to speed up the mapping process.
        """        
        with self.console.status('Create parquet partitions of raw mapping file'):
            # empty uniprot dictionary
            self.uniprot_dict = UniprotDict()

            # read mapping file with dask
            mapping_ddf = self.read_uniprot_mappings()

            # create folder if not existing
            if not os.path.exists(self.parquet_path):
                os.makedirs(self.parquet_path)
            # empty folder
            for files in glob.glob(os.path.join(self.parquet_path, '*')):
                os.remove(files)

            # write to parquet partitions
            mapping_ddf.to_parquet(self.parquet_path)

    def read_uniprot_mappings(self) -> dd.DataFrame:
        """
        Read the mapping file with the UniProt mappings with Dask.

        Returns:
            dd.dataframe: The Dask DataFrame with the mapping file.
        """        
        mapping_ddf = dd.read_table(self.mapping_file, 
                                sep = '\t', 
                                header = None, 
                                names = list(self.uniprot_dict.get_uniprot_dict().keys()), 
                                dtype = self.uniprot_dict.get_all_types(), 
                                usecols = self.uniprot_dict.get_columns_to_keep(),
                                blocksize = self.block_size)
        
        return mapping_ddf
    
    def read_parquet_partitions(self) -> dd.DataFrame:
        """
        Read the parquet partitions with Dask.

        Returns:
            dd.dataframe: The Dask DataFrame read from the parquet partitions.
        """        
        return dd.read_parquet(os.path.join(self.parquet_path), ignore_metadata_file=True)
    
    def filter_data(self) -> dd.DataFrame:
        """
        Filter the data with the data filter.

        Returns:
            dd.dataframe: The filtered Dask DataFrame.
        """        
        mappings_ddf = self.read_parquet_partitions()

        filter_ddf = mappings_ddf[
                    (mappings_ddf['UniProtKB-AC'].isin(self.uniprot_dbs['UniProtKB-AC']['targets'])) |
                    (mappings_ddf['UniProtKB-ID'].isin(self.uniprot_dbs['UniProtKB-ID']['targets'])) |
                    (mappings_ddf['RefSeq'].isin(self.uniprot_dbs['RefSeq']['targets'])) |
                    (mappings_ddf['GeneID'].isin(self.uniprot_dbs['GeneID']['targets'])) |
                    (mappings_ddf['Ensembl'].isin(self.uniprot_dbs['Ensembl']['targets'])) |
                    (mappings_ddf['UniParc'].isin(self.uniprot_dbs['UniParc']['targets']))
                        ]  
        
        return filter_ddf

    def finalize(self) -> None: 
        """
        Finalize the mapping of target sequences to the specified database:
            - Add the NCBI column.
            - Log the failed mappings.
            - Save the mapping DataFrame to a CSV file.        
        """        
        # possible existing target column is overwritten
        self.add_target_column()
        self.add_ncbi_column()
        self.log_failed_finalize()

        # write dataframe to csv file
        self.mapping_df.to_csv('mapping.csv', index=False) 
        self.console.print_done('All mapping done. Table saved to mapping.csv')
    
    def get_codes(self, id_type: str = None) -> list:
        """
        Get the list of codes from the mapping DataFrame.

        Args:
            id_type (str, optional): The id standard of the codes to get. Defaults to None using to_type attribute.

        Returns:
            list: The list of codes.
        """             
        return self.mapping_df[id_type].dropna().tolist()
    
    def get_target_to_result_dict(self) -> dict:
        """
        Get a dictionary with the target sequences as keys and the mapped results as values.

        Returns:
            dict: The dictionary with the target sequences and the mapped results.
        """        
        df = self.mapping_df[(self.mapping_df['target'].notna()) &
                              (self.mapping_df[self.to_type].notna())] 
        return df.set_index('target')[self.to_type].to_dict()
    
    def get_targets_and_ncbi_codes(self) -> list[tuple]:
        """
        Get a list of tuples with the target sequences and the NCBI codes.

        Returns:
            list[tuple]: The list of tuples with the target sequences and the NCBI codes.
        """        
        # return list of tuples with target and ncbi_code
        df = self.mapping_df[(self.mapping_df['target'].notna()) &
                              (self.mapping_df['ncbi_code'].notna())]        
        return list(zip(df['target'], df['ncbi_code']))
           
    def create_uniprot_dbs_with_targets(self) -> None:
        """
        Create the UniProt dictionary and assign the target sequences to the different type standards.
        """        
        self.uniprot_dict = UniprotDict()
        # fill uniprot_dbs with the targets        
        self.uniprot_dict.assign_target_to_type(self.target_list)       

    def make_unique(self, mapping_df: pd.DataFrame, from_type: list) -> pd.DataFrame:
        """
        Make the mapping DataFrame unique by removing duplicates. If tie, select the
        one with the shorte 'UniProtKB-AC' id.

        Args:
            mapping_df (pd.DataFrame): The DataFrame with the mapping results.
            from_type (list): The list of database id standards to map from.
            to_type (str): The database id standard to map to.

        Returns:
            pd.DataFrame: The DataFrame with the unique mapping results.
        """        
        # make copy before handling data frame
        df = mapping_df.copy()

        # add length colum for sorting
        # Rule of thumb, just take the first one for UniParc.            
        # The longer UniProt-KB are currated automatically. Take the samller one
        # pd.NA is a float, hence len() directly is not possible
        df['Length'] = df['UniProtKB-AC'].apply(lambda x: len(x) if pd.notna(x) else 0)
        # sort, default ascending
        df = df.sort_values(by='Length')         

        # do it for every column in from_type
        for col in from_type: 
            # is duplicated detects Nan as duplicated, avoid this
            mask = df[col].notna()
            # Checking for duplicates only on non-NaN values while preserving the original index
            is_duplicated = df.loc[mask, col].duplicated(keep=False).reindex(df.index, fill_value=False)
            # DataFrame of non-unique entries, selecting the first occurrence of each duplicated entry
            # sorting is not changed by group by
            first_non_unique_df = df[is_duplicated].groupby(col, sort=False).head(1)
            # DataFrame of unique entries (i.e., entries that are not duplicated)
            unique_df = df[~is_duplicated]
            # Combine the two DataFrames
            filtered_unique_df = pd.concat([unique_df, first_non_unique_df], ignore_index=True)
        
        # drop sorting column
        filtered_unique_df.drop(columns='Length', inplace=True)   
        return filtered_unique_df
    
    def add_target_column(self) -> None:
        """
        Add the target column to the mapping DataFrame and fill it with the target sequences.
        """        
        # add target column with NA or fill it with NA if present
        self.mapping_df['target'] = pd.NA
        for key in self.target_types:
            targets = self.uniprot_dict[key]['targets']
            # Update the 'target' column with the value from 'key' column where 'key' is in targets
            self.mapping_df.loc[self.mapping_df[key].isin(targets), 'target'] = self.mapping_df[key]

    def add_ncbi_column(self) -> None:
        """
        Add the NCBI column to the mapping.
        """        
        self.mapping_df['ncbi_code'] = pd.NA

        # Update 'ncbi_code' with values from 'RefSeq' if not NA, otherwise from 'EMBL-CDS'
        # np.where(condition, x, y) returns x if condition is True, otherwise y
        self.mapping_df['ncbi_code'] = np.where(self.mapping_df['RefSeq'].notna(), 
                                                self.mapping_df['RefSeq'], self.mapping_df['EMBL-CDS'])

    def log_failed_finalize(self) -> None:
        """
        Log the failed mappings and finalize the mapping of target sequences to the specified database.
        """        
        df = self.mapping_df.copy()

        # RefSeqs where no UniProtKB-AC was found
        ref_no_uniprot = df.loc[(df['RefSeq'].notna()) & (df['UniProtKB-AC'].isna()) ,'target'].to_list()
        if len(ref_no_uniprot) > 0:
            message = '{} RefSeq ids not mapped to UniProtKB-AC but included in NCBI-Code.'.format(len(ref_no_uniprot))
            self.console.print_warning(message)
            for id in ref_no_uniprot:
                logger.warning(f'Target sequence {id}') 

        # those where no UniProtKB-AC was found
        no_uniprot = df.loc[df['UniProtKB-AC'].isna(), 'target'].to_list()
        no_uniprot = list(set(no_uniprot).difference(set(ref_no_uniprot)))
        if len(no_uniprot) > 0:
            message = '{} ids not mapped to UniProtKB-AC.'.format(len(no_uniprot))
            self.console.print_warning(message)
            for id in no_uniprot:
                logger.warning(f'Target sequence {id}')   

        # those where no RefSeq or EMBL-CDS was found:
        no_ncbi = df.loc[df['ncbi_code'].isna(), 'target'].to_list()
        no_ncbi = list(set(no_ncbi).difference(set(no_uniprot)))
        if len(no_ncbi) > 0:
            message = '{} ids mapped to UniProtKB-AC but not to NCBI-Code.'.format(len(no_ncbi))
            self.console.print_warning(message)
            for id in no_ncbi:
                logger.warning(f'Target sequence {id}')          

    def log_failed(self) -> None:
        """
        Log the failed mappings of target sequences to the specified id standard.
        """        
        df = self.mapping_df.copy()
        no_hits = df.loc[df[self.to_type].isna(), 'target'].to_list()
        if len(no_hits) > 0:
            message = '{} ids not mapped to {}.'.format(len(no_hits),self.to_type)
            self.console.print_warning(message)
            for id in no_hits:
                logger.warning(f'Target sequence {id}')            