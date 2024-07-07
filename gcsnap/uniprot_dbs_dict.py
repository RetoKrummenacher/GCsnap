

class UniprotDict:
    
    def __init__(self) -> None:
        # create empty dictionary with mapping information
        self.get_empty_uniprot_dbs_dict()
        # get list of supported databases
        self.supported = self.get_supported_databases()

    def assign_target_to_type(self, target_list: list) -> None:        
        # add each target to the corresponding dictionary entry
        for target in target_list:       
            target = target.strip()    
            if 'UniRef' in target:
                # UniRef are clusters, just split the ending which is the UniProtKB-AC
                self.uniprot_dbs['UniProtKB-AC']['targets'].append(target.split('_')[1])
            elif target.startswith('UPI'):
                self.uniprot_dbs['UniParc']['targets'].append(target)
            elif target.startswith('ENSG'):
                self.uniprot_dbs['Ensembl']['targets'].append(target)
            elif target.isnumeric():
                self.uniprot_dbs['GeneID']['targets'].append(target)
            elif len(target.split('_')[0]) == 2:
                self.uniprot_dbs['RefSeq']['targets'].append(target)
            elif len(target.split('.')) == 2:
                self.uniprot_dbs['EMBL-CDS']['targets'].append(target)                
            elif len(target.split('_')) == 2:
                self.uniprot_dbs['UniProtKB-ID']['targets'].append(target)
            else:
                self.uniprot_dbs['UniProtKB-AC']['targets'].append(target)  

        # get list of target types
        self.get_target_types()                             

    def get_empty_uniprot_dbs_dict(self) -> None:
        # columns according to: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
        # it is not possible to map everything to everything: Check the online mapping tool https://www.uniprot.org/id-mapping
        # Every UniProtKB_AC-ID can be mapped to all, but for instance not GeneID to EMBL-CDS, only to Uniprot-KB or Swiss-Prot
        # Meaning, there is need to map anything to UniProtKB first and than to whatever you like.
        # However, it seems possible with the REST API: https://www.uniprot.org/help/id_mapping
        # UniRef is a comuted cluster from UniProt genomes, hence there are many hits.
        # The number (50, 90, 100) is the percentag of matching within the cluster.
        # To decrease size of the file uniprot mapping files to store, we included
        # only supported columns (supported == True).
        # UniRef is supported, but those information is not needed as this are
        # clusters already and we use the split part which is the UniProtKB-AC
        self.uniprot_dbs = {'UniProtKB-AC': {'targets' : [], 
                                        'from_dbs' : 'UniProtKB_AC-ID',
                                        'to_dbs' : 'UniProtKB',
                                        'dtype' : str,
                                        'supported': True},
                    'UniProtKB-ID': {'targets' : [], 
                                        'from_dbs' : 'UniProtKB_AC-ID',
                                        'to_dbs' : 'UniProtKB-Swiss-Prot',
                                        'dtype' : str,
                                        'supported': True},
                    'GeneID': {'targets' : [], #EntrezID
                                'from_dbs' : 'GeneID',
                                'to_dbs' : 'GeneID',
                                'dtype' : str,
                                'supported': True},
                    'RefSeq': {'targets' : [], 
                                'from_dbs' : 'RefSeq_Protein',
                                'to_dbs' : 'RefSeq_Protein', 
                                'dtype' : str,
                                'supported': True},
                    'GI': {'targets' : [],
                            'from_dbs' : 'GI_number',
                            'to_dbs' : 'GI_number', 
                            'dtype' : str,
                            'supported': False},
                    'PDB': {'targets' : [],
                            'from_dbs' : 'PDB', 
                            'to_dbs' : 'PDB', 
                            'dtype' : str,
                            'supported': False},
                    'GO': {'targets' : [],
                            'from_dbs' : '', 
                            'to_dbs' : '', 
                            'dtype' : str,
                            'supported': False},
                    'UniRef100': {'targets' : [],
                                    'from_dbs' : 'UniRef100', 
                                    'to_dbs' : 'UniRef100', 
                                    'dtype' : str,
                                    'supported': False},
                    'UniRef90': {'targets' : [],
                                    'from_dbs' : 'UniRef90', 
                                    'to_dbs' : 'UniRef90', 
                                    'dtype' : str,
                                    'supported': False},
                    'UniRef50': {'targets' : [],
                                    'from_dbs' : 'UniRef50',
                                    'to_dbs' : 'UniRef50',
                                    'dtype' : str,
                                    'supported': False},
                    'UniParc': {'targets' : [],
                                'from_dbs' : 'UniParc', 
                                'to_dbs' : 'UniParc', 
                                'dtype' : str,
                                'supported': True},
                    'PIR': {'targets' : [],
                            'from_dbs' : 'PIR', 
                            'to_dbs' : 'PIR', 
                            'dtype' : str,
                            'supported': False},
                    'NCBI-taxon': {'targets' : [],
                                    'from_dbs' : None, 
                                    'to_dbs' : None, 
                                    'dtype' : str,
                                    'supported': False},
                    'MIM': {'targets' : [],
                            'from_dbs' : 'MIM', 
                            'to_dbs' : 'MIM', 
                            'dtype' : str,
                            'supported': False},
                    'UniGene': {'targets' : [],
                                'from_dbs' : '', 
                                'to_dbs' : '', 
                                'dtype' : str,
                                'supported': False},
                    'PubMed': {'targets' : [], # https://pubmed.ncbi.nlm.nih.gov/10840038/
                                'from_dbs' : None, 
                                'to_dbs' : None, 
                                'dtype' : str,
                                'supported': False},
                    'EMBL': {'targets' : [],
                                'from_dbs' : 'EMBL-GenBank-DDBJ', 
                                'to_dbs' : 'EMBL-GenBank-DDBJ', 
                                'dtype' : str,
                                'supported': False},
                    'EMBL-CDS': {'targets' : [],
                                    'from_dbs' : 'EMBL-GenBank-DDBJ_CDS', 
                                    'to_dbs' : 'EMBL-GenBank-DDBJ_CDS', 
                                    'dtype' : str,
                                    'supported': True},
                    'Ensembl': {'targets' : [],
                                'from_dbs' : 'Ensembl',
                                'to_dbs' : 'Ensembl',
                                'dtype' : str,
                                'supported': True},
                    'Ensembl_TRS': {'targets' : [],
                                    'from_dbs' : 'Ensembl_Transcript', 
                                    'to_dbs' : 'Ensembl_Transcript', 
                                    'dtype' : str,
                                    'supported': False},
                    'Ensembl_PRO': {'targets' : [],
                                    'from_dbs' : 'Ensembl_Protein', 
                                    'to_dbs' : 'Ensembl_Protein',
                                    'dtype' : str,
                                    'supported': False},
                    'Additional PubMed': {'targets' : [], # https://pubmed.ncbi.nlm.nih.gov/10840038/
                                            'from_dbs' : None,
                                            'to_dbs' : None, 
                                            'dtype' : str,
                                            'supported': False}} 
        
    def get_supported_databases(self) -> list:
        return [key for key in self.uniprot_dbs if self.uniprot_dbs[key]['supported'] == True]
    
    def get_target_types(self) -> list:
        self.target_types = [key for key in self.uniprot_dbs if len(self.uniprot_dbs[key]['targets']) > 0]