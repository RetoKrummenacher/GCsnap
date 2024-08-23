import os
import sys
import time
import datetime
import pandas as pd

# import database handlers
from gcsnap.db_handler_uniprot_mappings import UniprotMappingsDBHandler
from gcsnap.uniprot_dbs_dict import UniprotDict


def read_uniprot_mappings(data_path, columns_to_keep = 'keep_in_df'):

    # default value of blocksize is computed as:
        # https://github.com/dask/dask/issues/1147
        # https://github.com/dask/dask/pull/1328
    # psutil.virtual_memory().total / psutil.cpu_count() / 10 / (2**20)
    # On an ACER Nitro 5 this is the upper limit of 64MB defined in Dask, as with 32GB Ram and 12 cores, block size is over 2 GB.
    # Even when devided by 10 its still 266 MB. So we set it manually to be larger. Then with have fewer parquet files
    # With default value of 64MiB, this results in 703 partitions, with parquet, each arround 90 MiB, we want to be between 100MB and 300 MB:
    # https://docs.dask.org/en/latest/dataframe-parquet.html#dataframe-parquet
    
    mapping_dict = UniprotDict()
    uniprot_dbs = mapping_dict.get_uniprot_dict()
    
    # only keep cols that are desired for the mapping to work with
    cols_to_keep = [key for key,values in uniprot_dbs.items() if values.get(columns_to_keep, False)]
    
    # datatypes
    dtypes = {key: values['dtype'] for key, values in uniprot_dbs.items()}

    mappings_df = pd.read_table(data_path, 
                                 sep='\t', 
                                 header=None, 
                                 names=list(uniprot_dbs), 
                                 dtype=dtypes, 
                                 usecols=cols_to_keep)

    return mappings_df

def create_dbs(path: str) -> None:    

    print('Creating Uniprot mappings database...')
    sys.stdout.flush()  # Flush stdout to ensure output is captured

    st1 = time.time()        
    # open database handler and create tables
    mapping_db_handler = UniprotMappingsDBHandler(path, 'db', 'uniprot_mappings.db')
    mapping_db_handler.create_tables()
    
    # read raw file
    mappings_df = read_uniprot_mappings(os.path.join(path, 'mappings','idmapping_selected.tab'))

    elapsed_time = time.time() - st1
    formatted_time = str(datetime.timedelta(seconds=round(elapsed_time)))    
    print('Table read with pandas {}'.format(formatted_time))
    st1 = time.time()

    # batch size
    batch_size = 10000

    # Compute the number of batches
    num_batches = int(len(mappings_df) / batch_size) 

    # Iteratively load and insert data in batches
    for i in range(num_batches):
        batch_start = i * batch_size
        batch_end = batch_start + batch_size

        # Extract a chunk of data
        batch_df = mappings_df.iloc[batch_start:batch_end]

        # Convert the batch DataFrame to a list of tuples
        batch = [tuple(row) for row in batch_df.values]

        # Insert the batch into the database
        mapping_db_handler.batch_insert_mappings(batch)

    # # Insert remaining rows that might not fit into a complete batch
    remaining_batch_df = mappings_df.loc[num_batches * batch_size:]
    remaining_batch = [tuple(row) for row in remaining_batch_df.values]
    mapping_db_handler.batch_insert_mappings(remaining_batch)

    elapsed_time = time.time() - st1
    formatted_time = str(datetime.timedelta(seconds=round(elapsed_time)))    
    print('All writing done in {}'.format(formatted_time))
    st1 = time.time()

    # Reindex the database
    mapping_db_handler.reindex()
    
    elapsed_time = time.time() - st1
    formatted_time = str(datetime.timedelta(seconds=round(elapsed_time)))    
    print('Index created in {}'.format(formatted_time))        

if __name__ == "__main__":    
    path = '/scicore/home/schwede/GROUP/gcsnap_db'
            
    st = time.time()
    create_dbs(path)
    
    elapsed_time = time.time() - st
    formatted_time = str(datetime.timedelta(seconds=round(elapsed_time)))
    
    print('Uniprot mapping DB createn in {}'.format(formatted_time))
