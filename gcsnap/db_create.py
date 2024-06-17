import os
import sys
import glob
import time
import datetime

# import database handlers
from gcsnap.db_handler_assemblies import AssembliesDBHandler
from gcsnap.db_handler_sequences import SequenceDBHandler

from gcsnap.utils import processpool_wrapper

# database limit
SEQ_DB_LIMIT = 15 * (2**30) # 15 GB


def split_into_batches(lst: list, batch_size: int = 1000):
    for i in range(0, len(lst), batch_size):
        # yield produces a generator
        yield lst[i:i + batch_size]
        
        
def split_into_parts(lst: list, n: int) -> list[list]:
    # k, how many times n fits in len(list), m the reminder
    q, r = divmod(len(lst), n)
    return [lst[i * q + min(i, r):(i + 1) * q + min(i + 1, r)] for i in range(n)]
    
    
def create_db_handlers(db_name_prefix: str, db_dir: str, db_counter: list[int], parallel_db: int) -> list:
    handlers = []
    for i in range(parallel_db):
        seq_db_name = db_name_prefix + '{:03d}'.format(db_counter + i) + '.db'
        # open sequence db handler
        seq_db_handler = SequenceDBHandler(db_dir, seq_db_name)
        handlers.append(seq_db_handler)
        
    return handlers


def create_db_handlers_last(db_name_prefix: str, db_dir: str, db_counter: list[int], parallel_db: int) -> list:
    handlers = []    
    for i in range(parallel_db):
        seq_db_name = db_name_prefix + '{:03d}'.format(db_counter + i) + '.db'
        # open sequence db handler
        seq_db_handler = SequenceDBHandler(db_dir, seq_db_name)
        handlers.append(seq_db_handler)
        
    return handlers


def execute_handler(args: tuple):
    handler, batch = args
    # each handler returns a list of tuples
    mapping = handler.insert_sequences_from_faa_files(batch)
    # add database name to result
    return [(item[0], item[1], handler.db_name) for item in mapping]


def execute_handlers(handlers: list[SequenceDBHandler], batch: list) -> list[tuple[str,str,str]]:
    # split the batch in equal sized parts
    batches = split_into_parts(batch, len(handlers)) 
    parallel_args = zip(handlers, batches)
    
    # retrurns a list of lists with tuples
    mapping_list = processpool_wrapper(len(handlers), parallel_args, execute_handler)
    
    # flatten the list of lists
    return [item for sublist in mapping_list for item in sublist]
       

def check_db_size(handlers: list[SequenceDBHandler], db_counter: int) -> list[int]:
    for i in range(len(handlers)):
             if handlers[i].get_db_size() > SEQ_DB_LIMIT:
                 db_counter += 1
        
    return db_counter


def parallel_reindex(args: tuple[str,str]) -> None:
    db_dir, seq_db_name = args
    handler = SequenceDBHandler(db_dir, seq_db_name)
    handler.reindex_sequences()
    return 1


def print_assemblies(assemblies: list[str]) -> None:
    for assembly in assemblies:
        assembly_accession = '_'.join(os.path.basename(assembly).split('_')[:2])
        print('Assembly {} done'.format(assembly_accession))
  

def create_dbs(path: str, parallel_db: int) -> None:    
    # where to create the databases
    db_dir = os.path.join(path, 'ncbi_db')        
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)
        
    # open assembly database handler and create tables
    asse_db_handler = AssembliesDBHandler(db_dir, 'assemblies.db')
    asse_db_handler.create_tables()
    
    # number of databases to write to in parallel
    batch_size = parallel_db * 1000
    # keep track of sequences and assemblies
    n_assemblies = 0 
    n_sequences = 0
        
    for loop_var in [('genbank','GCA_'),('refseq','GCF_')]:    
                 
        db_type, db_name_prefix = loop_var
                
        # 1. fill assemblies table from assembly_summary_{db_type}.txt
        # (assembly_accession, url, taxid, species)
        summary_file = os.path.join(path, db_type, 'assembly_summary_{}.txt'.format(db_type))
        asse_db_handler.insert_assemblies_from_summary(summary_file)
        print('Summary {} done'.format(os.path.basename(summary_file)))

        # 2. Extract information from .faa files       
        # a) (ncbi_code, assembly_accession and seq_db) go into the mapping table of assemblies.db
        # b) (ncbi_code, sequence) go into sepearte databases (as this will get huge)
            # We use new databses for next file when reaching the limit
        # Directory containing the .ffa files
        faa_data_dir = os.path.join(path, db_type, 'data')
        # list of all files to parse
        # file_paths = [os.path.join(faa_data_dir, file_name) for file_name 
        #               in os.listdir(faa_data_dir)]        
        file_paths = glob.glob(os.path.join(faa_data_dir,'*_protein.faa.gz'))
        
        # counter to keep track of databases
        db_counter = 1
        db_counter_old = 1
        
        # we loop over all those files in batches, each database takes 
        # indexing is switched off to speed up
        for batch in split_into_batches(file_paths, batch_size):
                        
            # for the last, just distributed evenly among the existing handlers
            if len(batch) < batch_size:   
                db_counter = db_counter_old                     

            # create 4 db handlers. Name of the database: e.g., GCA_001.db
            handlers = create_db_handlers(db_name_prefix, db_dir, db_counter, parallel_db)
            
            # start the parsing for each handler, in parallel 
            mapping_list = execute_handlers(handlers, batch)
            
            n_sequences += len(mapping_list)
            n_assemblies += len(batch)
            
            # add mappings to assembly db
            asse_db_handler.insert_mappings(mapping_list)           
            
            # increase db counter if size limit is reached
            db_counter_old = db_counter
            db_counter = check_db_size(handlers, db_counter)
            
            # Format numbers with thousand separators
            formatted_assemblies = "{:,}".format(n_assemblies)
            formatted_sequences = "{:,}".format(n_sequences)
            
            print_assemblies(batch)
            print('{} assemblies and {} sequences done so far'.format(formatted_assemblies, formatted_sequences))           
            
    print('All writing done')

    # reindex assembly database
    asse_db_handler.reindex_mappings()
    
    print('Assembly indexings done')

    # reindex all sequence databases   
    parallel_args = [(db_dir, db_name) for db_name in os.listdir(db_dir)
                     if db_name.startswith('G')]
    result = processpool_wrapper(parallel_db, parallel_args, parallel_reindex)
        
    print('{} sequence indexings done'.format(sum(result)))

    return n_assemblies, n_sequences                   

# Example usage
if __name__ == "__main__":
    
   
    n_processes = int(sys.argv[1])
    # n_processes = 4
    
    # set out path
    if os.name == 'nt':  # Windows
        path = r'C:\MT\data'
    else:  # Linux sciCORE
        path = '/scicore/home/schwede/GROUP/gcsnap_db'
            
    st = time.time()
    n_assemblies, n_sequences = create_dbs(path, n_processes)
    
    elapsed_time = time.time() - st
    formatted_time = str(datetime.timedelta(seconds=round(elapsed_time)))

    # Format numbers with thousand separators
    formatted_assemblies = "{:,}".format(n_assemblies)
    formatted_sequences = "{:,}".format(n_sequences)
    
    print('{} assemblies with {} sequences done in {} seconds'.
          format(formatted_assemblies, formatted_sequences, formatted_time))
    
        
    # some selecting
    # ncbi_codes = ['AAK02085.1','AAK02086.1','AAC70070.1','AAC70069.1']
    # res = db_handler.select(ncbi_codes)
    
    # # select all assemblies
    # res = db_handler.select_all('assemblies')
    
    # # select all sequences (DEFAULT)
    # res = db_handler.select_all()
    
    # # Retrieve a sample sequence
    # sequence = db_handler.get_sequence('AAV61620.1')
    # print(sequence)