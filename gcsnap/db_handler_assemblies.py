import os
import sqlite3

class AssembliesDBHandler:
    def __init__(self, db_path: str ,db_name: str ='assemblies.db'):
        self.db = os.path.join(db_path, db_name)
        
    def create_tables(self) -> None:
        self.create_assembly_table()
        self.create_mapping_table()
        self._disable_indices()
    
    def create_mapping_table(self) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS mappings')
        cursor.execute('''
            CREATE TABLE mappings (
                ncbi_code TEXT PRIMARY KEY,
                assembly_accession TEXT,
                sequence_db TEXT
            )
        ''')
        conn.commit()
        conn.close()
        
    def create_assembly_table(self) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS assemblies')
        cursor.execute('''
            CREATE TABLE assemblies (
                assembly_accession TEXT PRIMARY KEY,
                url TEXT,
                taxid INTEGER,                 
                species TEXT
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
        
    def reindex_mappings(self):
        self._enable_indices()        
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('REINDEX mappings')
        conn.commit()
        conn.close()    
        
    def _reindex_assemblies(self):
        self._enable_indices()
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('REINDEX assemblies')
        conn.commit()
        conn.close()         

    def _batch_insert_mappings(self, mappings: list[tuple[str,str,str]]) -> None:
        self._disable_indices()
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.executemany('INSERT OR REPLACE INTO mappings (ncbi_code, assembly_accession, sequence_db) VALUES (?, ?, ?)', mappings)
        conn.commit()
        conn.close()
            
    def _batch_insert_assemblies(self, assemblies: list[tuple[str,str,int,str]]) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.executemany('INSERT OR REPLACE INTO assemblies (assembly_accession, url , taxid, species) VALUES (?, ?, ?, ?)', assemblies)
        conn.commit()
        conn.close()                       
        
    def _parse_and_insert_assemblies(self, file_path: str) -> None:     
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            
        assembly_list = []    
        
        for line in lines:
            if not line.startswith('#'):
                data = line.strip().split('\t')
                # the first column has the Assemlby id, column 19 the url to it
                assembly_accession = data[0]
                taxid = data[5]
                species = data[7]
                url = data[19]
                
                assembly_list.append((assembly_accession, url, taxid, species))        

                if len(assembly_list) >= 100000:  # Adjust batch size as needed
                    self._batch_insert_assemblies(assembly_list)
                    assembly_list = []

        # write the last part not complete 1000 assemblies
        if assembly_list:
            self._batch_insert_assemblies(assembly_list)     
            
        # index assemblies table and set PRAGMA back to fast wrtiting
        self._reindex_assemblies()
        self._disable_indices()
        
    def insert_assemblies_from_summary(self, file_path: str) -> None:
        self._parse_and_insert_assemblies(file_path)
        
    def insert_mappings(self, mappings: list[tuple[str,str,str]]) -> None:
        self._batch_insert_mappings(mappings)
                
    def select(self, ncbi_codes: list[str], return_fields: list[str] = None, table: str = 'mappings') -> list[tuple]:
        # combine query
        
        # which fields to return
        if return_fields is None:
            select_fields = '*'
        else:
            select_fields = ', '.join(return_fields)

        # which records based on ncbi_codes
        records = f"({','.join(['?']*len(ncbi_codes))})" # (?,?,?,?) for each entry in ncbi_code       
        query = 'SELECT {} FROM {} WHERE ncbi_code IN {}'.format(select_fields,table,records)
        
        # execute query
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute(query, ncbi_codes)
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()
        
        return result
    
    def select_all(self, table: str = 'mappings') -> list[tuple]:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute(f'SELECT * FROM {table}')
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()
        
        return result
