"""
This script is used to update the sequence databases.
The original db_create_assemblies.py scriptwas using a wrong implementation
of dh_handler_sequences.py. This was correct but the sequences need to be updated.
The Primary Key of the sequences table is the ncbi_code, so if the same ncbi_code.
The indexes are already created, so we just update the seqeucences.
"""

import os
import sys
import glob
import time
import datetime

# import database handlers
from gcsnap.db_handler_sequences import SequenceDBHandler
from gcsnap.utils import processpool_wrapper


def split_into_batches(lst: list, batch_size: int = 1000):
    for i in range(0, len(lst), batch_size):
        # yield produces a generator
        yield lst[i:i + batch_size]
        
        
def split_into_parts(lst: list, n: int) -> list[list]:
    # k, how many times n fits in len(list), m the reminder
    q, r = divmod(len(lst), n)
    return [lst[i * q + min(i, r):(i + 1) * q + min(i + 1, r)] for i in range(n)]
    

def execute_handler(args: tuple):
    handler, batch = args
    # each handler returns a list of tuples
    sequences, _ = handler.parse_sequences_from_faa_files(batch)
    return sequences


def execute_handlers(handler: SequenceDBHandler, batch: list, n_processes: int) -> tuple[list[tuple[str,str]], list[tuple[str,str]]]:
    # split the batch in equal sized parts
    batches = split_into_parts(batch, n_processes) 
    parallel_args = [(handler, subbatch) for subbatch in batches]
    
    # retrurns a list of tuples containing lists of tuples
    result = processpool_wrapper(n_processes, parallel_args, execute_handler)
    
    # flatten the list of tuples containing lists of tuples
    #mappings = [item for sublist in result for item in sublist[1]]
    sequences = [item for sublist in result for item in sublist]
    
    return sequences

def update_dbs(path: str, n_processes: int) -> None:    
    # where to create the databases
    db_dir = os.path.join(path, 'db')           
    # open assembly database handler and create tables
    seq_db_handler = SequenceDBHandler(db_dir)   
    seq_db_handler.disable_indices()
    
    # number of files to write to in parallel
    batch_size = n_processes * 500
    # keep track of sequences and assemblies
    n_assemblies = 0 
    n_sequences = 0
        
    #for loop_var in ['genbank','refseq']:    
                 
        #db_type = loop_var

    # 1. Extract information from .faa files       
    # a) (ncbi_code, sequence) go into sepearte databases (as this will get huge)

    # Directory containing the .ffa files
    #faa_data_dir = os.path.join(path, db_type, 'data')
    # list of all files to parse      
    #file_paths = glob.glob(os.path.join(faa_data_dir,'*_protein.faa.gz'))

    # read needed assembly files for experiments: file was created handisch with assemblies_for_cluster.ipynb
    with open('/scicore/home/schwede/kruret00/MT/assemblies.txt', 'r') as file:
        content = file.read()
    lines = content.splitlines()
    file_names = [line.strip() + '_protein.faa.gz' for line in lines]

    # add path structer
    file_paths =  [os.path.join(path, 'refseq', 'data', file_name) for file_name in file_names if file_name.startswith('GCF')]
    file_paths += [os.path.join(path, 'genbank', 'data', file_name) for file_name in file_names if file_name.startswith('GCA')]
            
    # we loop over all those files in batches, each database takes 
    # indexing is switched off to speed up
    for batch in split_into_batches(file_paths, batch_size):
        
        # start the parsing for each handler, in parallel 
        sequence_list = execute_handlers(seq_db_handler, batch, n_processes)
        
        # keep track of done things
        n_sequences += len(sequence_list)
        n_assemblies += len(batch)
                
        # add sequences
        seq_db_handler.batch_update_sequences(sequence_list)      
                    
        # Format numbers with thousand separators
        formatted_assemblies = "{:,}".format(n_assemblies)
        formatted_sequences = "{:,}".format(n_sequences)
        
        print('{} assemblies and {} sequences done so far'.format(formatted_assemblies, formatted_sequences))           
            
    print('All Updating done')
    
    # get number of unique sequences
    n_unique = seq_db_handler.select_number_of_entries()
    seq_db_handler.enable_indices()

    return n_assemblies, n_sequences, n_unique                   

# Example usage
if __name__ == "__main__":    
   
    n_processes = int(sys.argv[1])
    # n_processes = 10

    path = '/scicore/home/schwede/GROUP/gcsnap_db'
            
    st = time.time()
    n_assemblies, n_sequences, n_unique = update_dbs(path, n_processes)
    
    elapsed_time = time.time() - st
    formatted_time = str(datetime.timedelta(seconds=round(elapsed_time)))

    # Format numbers with thousand separators
    formatted_assemblies = "{:,}".format(n_assemblies)
    formatted_sequences = "{:,}".format(n_sequences)
    formatted_unique = "{:,}".format(n_unique)
    
    print('{} assemblies with {} sequences ({} unique sequences found) done in {}'.
          format(formatted_assemblies, formatted_sequences, formatted_unique, formatted_time))
 