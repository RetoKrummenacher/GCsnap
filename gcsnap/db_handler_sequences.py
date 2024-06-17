import os
import gzip
import sqlite3

class SequenceDBHandler:
    def __init__(self, db_path: str ,db_name: str):
        self.db = os.path.join(db_path, db_name)
        self.db_name = db_name
        self._check_and_create()
        
    def _check_and_create(self) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='sequences'")
        result = cursor.fetchone()
        conn.close()
        
        if result is None:
            self.create_sequence_table()
            self._disable_indices()
    
    def create_sequence_table(self) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS sequences')
        cursor.execute('''
            CREATE TABLE sequences (
                ncbi_code TEXT PRIMARY KEY,
                sequence TEXT
            )
        ''')
        conn.commit()
        conn.close()    
        
    def _disable_indices(self):
        conn = sqlite3.connect(self.db)
        # Disabling indices and other performance-related settings
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('PRAGMA journal_mode = OFF')
        conn.execute('PRAGMA temp_store = MEMORY')
        conn.execute('PRAGMA cache_size = -1048576')  # 1 GB cache size (1,048,576 KB) 
        conn.commit()
        conn.close()
        
    def _enable_indices(self):
        conn = sqlite3.connect(self.db)
        # Re-enabling indices and other settings
        conn.execute('PRAGMA synchronous = NORMAL')
        conn.execute('PRAGMA journal_mode = WAL')
        conn.execute('PRAGMA temp_store = DEFAULT')        
        conn.execute('PRAGMA cache_size = -2000')  # Reset to default cache size  
        conn.commit()
        conn.close()      
        
    def reindex_sequences(self):
        self._enable_indices()
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('REINDEX sequences')
        conn.commit()
        conn.close()         

    def _batch_insert_sequences(self, sequences: list[tuple[str,str]]) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('PRAGMA journal_mode = WAL')
        cursor.executemany('INSERT OR REPLACE INTO sequences (ncbi_code, sequence) VALUES (?, ?)', sequences)
        conn.commit()
        conn.close()         

    def _read_gzip_file(self, file_path: str) -> str:
        with gzip.open(file_path, 'rt', encoding='utf-8') as file:
            content = file.read()
        return content.splitlines()

    def _read_txt_file(self, file_path: str) -> str:
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        return lines
            
    def _parse_and_insert_sequences(self, file_paths: list[str]) -> list[tuple[str,str]]:
        sequence_list = []  
        mapping_list = []
        for file_path in file_paths:
            # accession taken from file name
            assembly_accession = '_'.join(os.path.basename(file_path).split('_')[:2])
                    
            # read the file in once
            if file_path.endswith('.gz'):
                lines = self._read_gzip_file(file_path)
            else:
                lines = self._read_txt_file(file_path)
            
            # parse the lines, the challange, the sequence can be several lines long
            # make one string of it
            content = ''.join(lines)
            # split that string to extract each sequence id
            # EFB12766.1 hypothetical protein PANDA_022614, partial [Ailuropoda melanoleuca]WSDGHLIYYDDQTRQSVEDKVHMPVDCINIRTGHECRGT
            # the first is an empty result
            entries = content.split('>')[1:]
            
            for entry in entries:
                # split the info from the acutal sequence str
                entry_split = entry.split('\n')
                sequence = ''.join(entry_split[1:])
                # split the organism name in []
                info_split = entry_split[0].split('[') 
                # orgname is same for the assembly, put it into assembly table
                #orgname = info_split[-1].replace(']','')
                ncbi_code = info_split[0].split(' ')[0]
                
                sequence_list.append((ncbi_code, sequence))
                mapping_list.append((ncbi_code, assembly_accession))
        
        # insert all as batch
        self._batch_insert_sequences(sequence_list)     
        
        # return what is needed from the .faa file
        return mapping_list      

    def insert_sequences_from_faa_files(self, file_paths: list[str]) -> list[tuple[str,str]]:
        # new use to do this in parallel, as the database writing is the bottleneck
        # SQLite does not support parallel write, so just batches to reduce write calls          
        return self._parse_and_insert_sequences(file_paths) 
                    
    def select(self, ncbi_codes: list[str], return_fields: list[str] = None) -> list[tuple]:
        # combine query
        
        # which fields to return
        if return_fields is None:
            select_fields = '*'
        else:
            select_fields = ', '.join(return_fields)

        # which records based on ncbi_codes
        records = f"({','.join(['?']*len(ncbi_codes))})" # (?,?,?,?) for each entry in ncbi_code       
        query = 'SELECT {} FROM sequences WHERE ncbi_code IN {}'.format(select_fields,records)
        
        # execute query
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute(query, ncbi_codes)
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()  
        
        return result
    
    def select_all(self) -> list[tuple]:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM sequences')
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()   
        
        return result
    
    def get_db_size(self) -> int:
        return os.path.getsize(self.db)