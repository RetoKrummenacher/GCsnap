## GCsnap.py - devoloped by Joana Pereira, Structural Computational Biology, Biozentrum University of Basel, Basel Switzerland
## Last changed: 18.08.2023

import subprocess as sp
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

import matplotlib
import json
import os
import time
import random
import urllib.parse
import urllib.request
import requests
import sys
import statistics
import argparse
import requests_cache
import pacmap

# for UniProt id mapping
import re
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import Phylo
from scipy.cluster import hierarchy
from scipy.spatial import distance
from scipy import stats
from collections import Counter
from multiprocessing.pool import ThreadPool

from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture

from bokeh.plotting import figure, output_file, gridplot, save
from bokeh.colors import RGB
from bokeh.models import HoverTool, TapTool, LassoSelectTool, Range1d, LinearAxis, WheelZoomTool, Circle, MultiLine, Panel, Tabs
from bokeh.models import ColumnDataSource, DataTable, DateFormatter, TableColumn, Legend, HTMLTemplateFormatter
from bokeh.models.callbacks import OpenURL
from bokeh.models.graphs import from_networkx, NodesAndLinkedEdges
from bokeh.models.widgets import Div
from bokeh.layouts import row, column

# for timings
from codetiming import Timer
from humanfriendly import format_timespan
from datetime import datetime

now = datetime.now()
current_time = now.strftime("%d-%m-%Y_%H-%M-%S")
timing_log_file = f'times_{current_time}.tsv'

# import bokeh.layouts 

plt.switch_backend('agg')

# set parameters for uniprot requests

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504, 400])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# HELPING ROUTINES

# 0. General routines

def chunk_list(l, n):


	chunks = np.array_split(np.array(l), n)

	chunks = [list(chunk) for chunk in chunks]
	return chunks

@Timer("find_clusters_in_distance_matrix", text=lambda secs: f"Time to find clusters in distance matrix: {format_timespan(secs)}")
def find_clusters_in_distance_matrix(distance_matrix, t = 0):

	distance_matrix = distance.squareform(distance_matrix)
	linkage = hierarchy.linkage(distance_matrix, method = 'single')
	clusters = hierarchy.fcluster(linkage, t, criterion = 'distance')
	clusters = [int(i) for i in clusters]

	return clusters

#@Timer("mask_singleton_clusters", text=lambda secs: f"Time to mask singleton clusters: {format_timespan(secs)}")
def mask_singleton_clusters(clusters_list, mask = 0):
	
	new_clusters_list = []

	for value in clusters_list:
		if list(clusters_list).count(value) == 1:
			new_clusters_list.append(mask)
		else:
			new_clusters_list.append(value)

	return new_clusters_list

def merge_intervals(intervals):

	intervals = np.array(intervals)
	intervals = intervals[intervals[:, 0].argsort()]	

	new_intervals = []
	for interval in intervals:
		if len(new_intervals) == 0:
			new_intervals.append(interval)
		else:
			previous_interval = list(new_intervals[-1])
			if interval[0] <= previous_interval[1]:
				overlap_range = interval+previous_interval
				new_intervals[-1] = [min(overlap_range), max(overlap_range)]
			else:
				new_intervals.append(interval)

	return new_intervals

# 1. Routines to get the assembly for a given ncbi_code (or entrezID), download it and parse it

#@Timer("map_uniprot_to_ncbi", text=lambda secs: f"Time to map uniprot to ncbi: {format_timespan(secs)}")
def map_uniprot_to_ncbi(uniprot_code, search_database = 'EMBL-GenBank-DDBJ_CDS'):

	if 'UniRef' in uniprot_code:
		uniprot_label = uniprot_code.split('_')[-1]
	elif uniprot_code.startswith('UPI'):  # if it is a UniParc ID
		uniprot_label = map_codes_to_uniprot([uniprot_code], from_database = 'UniParc')[uniprot_code]
	elif uniprot_code.startswith('ENSG'): # if it is an ensemble gene ID
		uniprot_label = map_codes_to_uniprot([uniprot_code], from_database = 'Ensembl')[uniprot_code]
	elif uniprot_code.isnumeric():		# it is a gene id
		uniprot_label = map_codes_to_uniprot([uniprot_code], from_database = 'GeneID')[uniprot_code]
	else:
		uniprot_label = uniprot_code
	
	results = get_mappings_through_uniprot([uniprot_code], from_database = 'UniProtKB_AC-ID', to_database = search_database)
	
	ncbi_code = 'nan'
	if results != []:
		ncbi_code = results[0]['to']

	if ncbi_code == 'nan':
		if search_database != 'RefSeq_Protein':
			ncbi_code, search_database =  map_uniprot_to_ncbi(uniprot_code, search_database = 'RefSeq_Protein')
	else:
		print(" ... {} corresponds to {} in {} database".format(uniprot_code, ncbi_code, search_database))

	return ncbi_code, search_database

#@Timer("find_ncbi_code_assembly", text=lambda secs: f"Time to find ncbi code assembly: {format_timespan(secs)}")
def find_ncbi_code_assembly(ncbi_code, database_assembly_mapping):

	assembly_id = 'nan'
	assembly_source = 'nan'
	assembly_link = 'nan'
	uniprot_code = 'nan'
	search_database = 'nan'

	if '.' not in ncbi_code:
		uniprot_code = ncbi_code
		t1 = Timer('map_uniprot_to_ncbi',
				   text=lambda secs: f"Time to map uniprot id to ncbi code for {uniprot_code}: {format_timespan(secs)}")
		t1.start()
		ncbi_code, search_database = map_uniprot_to_ncbi(uniprot_code)
		t1.stop()

	try:
		t2 = Timer('get_assembly_id_from_ncbi_code',
				   text=lambda secs: f"Time to get assembly id from ncbi code for {ncbi_code}: {format_timespan(secs)}")
		t2.start()
		handle = Entrez.efetch(db="protein", id=ncbi_code, rettype="ipg", retmode="xml")
		record = Entrez.read(handle)

		print(' ... ... Searching for assembly for {}'.format(ncbi_code))

		assemblies_found = {}
		for report in record:
			if 'ProteinList' in record[report]:
				products = record[report]['ProteinList']
				for product in products:
					if 'CDSList' in product and "assembly" in str(product):
						cds = product['CDSList']
						cds = str(cds)
						cds = cds.replace('[','').replace(']','').replace('{', '').replace('}', '').replace('(','').replace(')','')

						assemblyId = cds.split("assembly")
						assemblyId = assemblyId[-1].split()[1].replace(' ','').replace("'",'')
						assemblyId = assemblyId.replace(',','')

						protein_acc = str(product).split("}, attributes={'accver': ")[1].split(',')[0].replace("'","")
						source = str(product).split("}, attributes={'accver': ")[1].split(", 'source': ")[1].split(',')[0].replace("'","")
						
						if protein_acc == ncbi_code:
							assembly_id = assemblyId
							assembly_source = source

							assembly_link = database_assembly_mapping[assembly_id]
		t2.stop()
	except:
		assembly_id = 'nan'
		assembly_source = 'nan'
		assembly_link = 'nan'

	if uniprot_code != 'nan' and assembly_id == 'nan' and search_database == 'EMBL':
		ncbi_code, search_database = map_uniprot_to_ncbi(uniprot_code, search_database = 'P_REFSEQ_AC')
		ncbi_code, assembly_id, assembly_link = find_ncbi_code_assembly(ncbi_code, database_assembly_mapping)

	return ncbi_code, assembly_id, assembly_link

#@Timer("download_and_extract_assembly", text=lambda secs: f"Time to download and extract assembly: {format_timespan(secs)}")
def download_and_extract_assembly(assembly_id, assembly_link, tmp_folder = None, label = '', get_chunk = False, chunk_size = 0, target = None):

	print(' ... ... Downloading and extracting assembly {} annotated gff file'.format(assembly_id))

	assembly_label = assembly_link.split('/')[-1]
	
	assembly_genomic_file_link = '{}/{}_genomic.gff.gz'.format(assembly_link, assembly_label)
	out_assembly_file_gz = '{}/{}'.format(tmp_folder, assembly_genomic_file_link.split('/')[-1])
	out_assembly_file = out_assembly_file_gz[:-3]

	if not os.path.isfile(out_assembly_file) and not os.path.isfile(out_assembly_file_gz):
		download_assembly = sp.Popen(['wget', assembly_genomic_file_link, '-O', out_assembly_file_gz], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		t1 = Timer('download_assembly',
				   text=lambda secs: f"Time to download assembly {assembly_id}: {format_timespan(secs)}")
		t1.start()
		stdout, stderr = download_assembly.communicate()
		t1.stop()

	if os.path.isfile(out_assembly_file_gz) or os.path.isfile(out_assembly_file):

		print(' ... ... ... Downloaded assembly {} ....'.format(assembly_id))

		extract_assembly_file = sp.Popen(['gunzip', out_assembly_file_gz], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		t2 = Timer('extract_assembly_file',
				   text=lambda secs: f"Time to extract assembly file {assembly_id}: {format_timespan(secs)}")
		t2.start()
		stdout, stderr = extract_assembly_file.communicate()
		t2.stop()

		if os.path.isfile(out_assembly_file):

			t3 = Timer('parse_assembly',
					   text=lambda secs: f"Time to parse assembly {assembly_id}: {format_timespan(secs)}")
			t3.start()
			assembly = parse_assembly(out_assembly_file, get_chunk = get_chunk, chunk_size = chunk_size, target = target)
			t3.stop()

			if not get_chunk:
				print(' ... ... ... Finished parsing assembly {} and found {} CDS entries'.format(assembly_id, len(assembly['ncbi_codes'])))
			else:
				print(' ... ... ... Finished parsing assembly {} and collected {} CDS entries around the target'.format(assembly_id, len(assembly['ncbi_codes'])))
			return assembly

		else:
			print(' ... ... ... It was not possible to extract assembly {}'.format(assembly_id))
			return 'nan'
	else:
		print(' ... ... ... It was not possible to save assembly {}'.format(assembly_id))
		return 'nan'

#@Timer("parse_assembly", text=lambda secs: f"Time to parse assembly: {format_timespan(secs)}")
def parse_assembly(assembly_file, get_chunk = False, chunk_size = 0, target = None):

	assembly = {'ncbi_codes': [], 'starts': [], 'ends': [], 'directions': [], 'names': [], 'scaffolds': []}
	curr_scaffold = 0
	found_chunk = False
	chunk_timer = 0
	
	with open(assembly_file, 'r') as in_assembly:
		for line in in_assembly:
			if not line.startswith('#'):
				line_data = line.split('\t')

				if line_data[2] == 'CDS':

					start = int(line_data[3])
					end = int(line_data[4])
					direction = line_data[6]

					if 'cds-' in line:
						ncbi_code = line_data[8].split('ID=cds-')[1].split(';')[0]
					else:
						if 'Name=' in line:
							ncbi_code = line_data[8].split('Name=')[1].split(';')[0]
						else:
							ncbi_code = 'unk'

					if 'pseudo=' not in line and 'product=' in line and 'fragment' not in line:
						prot_name = line_data[8].split('product=')[1].split(';')[0]
					else:
						prot_name = 'pseudogene'

					if ncbi_code in assembly['ncbi_codes'] and assembly['ncbi_codes'][-1] == ncbi_code: # it means that this is some kind of fragmented gene (has introns?...) and so we have to collect the largest interval encompassing it
						if start < assembly['starts'][-1]:
							assembly['starts'][-1] = start
						if end > assembly['ends'][-1]:
							assembly['ends'][-1] = end
					else:
						if '|' in ncbi_code:
							ncbi_code = ncbi_code.replace('|','_')

						assembly['ncbi_codes'].append(ncbi_code)
						assembly['names'].append(prot_name)
						assembly['scaffolds'].append(curr_scaffold)
						assembly['starts'].append(start)
						assembly['ends'].append(end)
						assembly['directions'].append(line_data[6])

					if get_chunk:
						if ncbi_code == target:
							chunk_timer = 1
						elif chunk_timer > 0:
							chunk_timer += 1

						if chunk_timer == round(chunk_size/2):
							break
					
			elif line.startswith('##sequence-region'):
				curr_scaffold += 1

	if get_chunk:
		chunked_assembly = {}
		for key in assembly:
			chunked_assembly[key] = assembly[key][-chunk_size:]
		assembly = chunked_assembly

	return assembly				
				
# 2. Routines to get the N flanking genes for a given ncbi_code 

#@Timer("get_n_flanking_genes", text=lambda secs: f"Time to get {n_flanking5 + n_flanking3} flanking_genes: {format_timespan(secs)}")
def get_n_flanking_genes(target_ncbi_code, assembly,  n_5 = None, n_3 = None, exclude_partial = None):

	if n_5 == n_3:
		print(' ... ... Extracting {} flanking genes ({} to each side) of {}'.format(n_5+n_3, n_5, target_ncbi_code))
	else:
		print(" ... ... Extracting {} flanking genes ({} to 5' and {} to 3') of {}".format(n_5+n_3, n_5, n_3, target_ncbi_code))
	
	flanking_genes = {'relative_starts': [], 'relative_ends': []}

	if target_ncbi_code in assembly['ncbi_codes']:
		index_of_target = assembly['ncbi_codes'].index(target_ncbi_code)
		direction_of_target = assembly['directions'][index_of_target]

		if direction_of_target == '+':
			genomic_context_block = [index_of_target-n_5, index_of_target+n_3]
		else:
			genomic_context_block = [index_of_target-n_3, index_of_target+n_5]

		for i in range(genomic_context_block[0], genomic_context_block[1]+1):
			if i>= 0 and i < len(assembly['scaffolds']):
				if assembly['scaffolds'][i] == assembly['scaffolds'][index_of_target]:
					for key in assembly.keys():
						if key != 'scaffolds':
							if key not in flanking_genes:
								flanking_genes[key] = []

							flanking_genes[key].append(assembly[key][i])

							if key == 'starts' or key == 'ends':
								flanking_genes['relative_{}'.format(key)].append(assembly[key][i] - assembly['starts'][index_of_target] + 1)

		print(' ... ... ... Found {} flanking genes for {}'.format(len(flanking_genes['ncbi_codes'])-1, target_ncbi_code))

		if direction_of_target == '-':
			index_of_target_in_flanking = flanking_genes['ncbi_codes'].index(target_ncbi_code)
			current_relative_starts = flanking_genes['relative_starts']
			current_relative_ends = flanking_genes['relative_ends']

			flanking_genes['relative_starts'] = [value*(-1)+current_relative_ends[index_of_target_in_flanking]+1 for value in current_relative_ends]
			flanking_genes['relative_ends'] = [value*(-1)+current_relative_ends[index_of_target_in_flanking]+1 for value in current_relative_starts]

			for key in flanking_genes:
				flanking_genes[key] = flanking_genes[key][::-1]
				if key == 'directions':
					flanking_genes[key] = list(''.join(flanking_genes[key]).replace('+','p'))
					flanking_genes[key] = list(''.join(flanking_genes[key]).replace('-','+'))
					flanking_genes[key] = list(''.join(flanking_genes[key]).replace('p','-'))

		if exclude_partial and len(flanking_genes['ncbi_codes']) < n_5+n_3+1:
			print(' ... ... ... This represents a partial syntenic block and as "exclude_partial" is on, it will not be considered.')
			return 'nan'
		else:   
			return flanking_genes
		
	else:
		print(' ... ... ... {} was not found in assembly'.format(target_ncbi_code))
		return 'nan'

#@Timer("add_sequences_to_flanking_genes", text=lambda secs: f"Time to add sequences to flanking_genes: {format_timespan(secs)}")
def add_sequences_to_flanking_genes(flanking_genes, target_ncbi_code):

	print(' ... ... Collecting sequences for flanking proteins')

	flanking_genes['sequences'] = []
	flanking_genes['species'] = []
	flanking_genes['taxID'] = []

	for ncbi_code in flanking_genes['ncbi_codes']:
		seq = ''
		species = ''
		uniprot_code = ''
		swiss_model_link = ''
		taxid = ''
		try:
			handle = Entrez.efetch(db="protein", id=ncbi_code, retmode="xml", rettype='fasta')
			record = Entrez.read(handle)
			
			for report in record:
				if seq == '':
					seq = report['TSeq_sequence']
				if species == '' and ncbi_code == target_ncbi_code:
					taxid = report['TSeq_taxid']
					species = report['TSeq_orgname']
		except:
			seq = 'FAKESEQUENCEFAKESEQUENCEFAKESEQUENCEFAKESEQUENCE'

		flanking_genes['sequences'].append(seq)

		if species != '':
			flanking_genes['species'] = species
		if taxid != '':
			flanking_genes['taxID'] = taxid

	return flanking_genes

# 3. Routines to compute all-against-all distance matrix and find protein families 

#@Timer("write_flanking_sequences_to_fasta", text=lambda secs: f"Time to write flanking sequences to fasta: {format_timespan(secs)}")
def write_flanking_sequences_to_fasta(all_syntenies, out_dir, out_label, exclude_pseudogenes = False, mode = 'flanking_sequences'):

	out_fasta = '{}/{}_{}.fasta'.format(out_dir, out_label, mode)
	seqs_lens = {}

	with open(out_fasta, 'w') as outfst:
		for target in all_syntenies.keys():

			if mode == 'flanking_sequences':
				flanking_genes = all_syntenies[target]['flanking_genes']
				for i, ncbi_code in enumerate(flanking_genes['ncbi_codes']):

					if flanking_genes['names'][i] != 'pseudogene' or not exclude_pseudogenes:
						outfst.write('>{}|{}\n{}\n'.format(ncbi_code, flanking_genes['names'][i], flanking_genes['sequences'][i]))
						seqs_lens[ncbi_code] = len(flanking_genes['sequences'][i])
			
			if mode == 'operon':
				outfst.write('>{}\n'.format(target))
				for sequence in all_syntenies[target]['flanking_genes']['sequences']:
					outfst.write(sequence)
				outfst.write('\n')

	return out_fasta, seqs_lens

# 3.1. Routines for when the method is based on BLASTp

@Timer("make_blast_database_from_fasta", text=lambda secs: f"Time to make blast database from fasta: {format_timespan(secs)}")
def make_blast_database_from_fasta(infasta, blast = None):

	print(' ... ... Making BLAST database')
	blastDB = '{}_blastDB'.format(infasta[:-6])

	make_blastDB = sp.Popen(['makeblastdb','-in', infasta, '-dbtype', 'prot', '-out', blastDB], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
	stdout, stderr = make_blastDB.communicate()
	
	if len(stderr) > 0:
		print(stderr)
		if 'BLAST engine error' not in stderr:
			make_blastDB = sp.Popen(['{}makeblastdb'.format(blast),'-in', infasta, '-dbtype', 'prot', '-out', blastDB], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
			stdout, stderr = make_blastDB.communicate()

	return blastDB

@Timer("run_blast_for_flanking_sequences", text=lambda secs: f"Time to run blast for flanking sequences: {format_timespan(secs)}")
def run_blast_for_flanking_sequences(seq_fasta, database, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, blast = None, tmp_folder = None):

	print(' ... ... Running BLAST')
	blast_outfile = '{}/{}_{}.xml'.format(tmp_folder, seq_fasta.split('/')[-1][:-6], max_evalue)

	if not os.path.isfile(blast_outfile):
		run_blast = sp.Popen(['psiblast', '-query', seq_fasta, '-db', database, '-out', blast_outfile, '-num_threads', str(num_threads), '-evalue', str(max_evalue), '-inclusion_ethresh', str(max_evalue), '-num_iterations', str(num_iterations),'-num_alignments', str(num_alignments), '-outfmt', '5'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		stdout, stderr = run_blast.communicate()

		if len(stderr) > 0:
			print(stderr)
			if 'BLAST engine error' not in stderr:
				run_blast = sp.Popen(['{}psiblast'.format(blast), '-query', seq_fasta, '-db', database, '-out', blast_outfile, '-num_threads', str(num_threads), '-evalue', str(max_evalue), '-inclusion_ethresh', str(max_evalue), '-num_iterations', str(num_iterations),'-num_alignments', str(num_alignments), '-outfmt', '5'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
				stdout, stderr = run_blast.communicate()

	return blast_outfile

@Timer("extract_distance_matrix_from_blast_output", text=lambda secs: f"Time to extract distance matrix from blast output: {format_timespan(secs)}")
def extract_distance_matrix_from_blast_output(blast_results, default_base = None, min_coverage = 70, sequences_lengths = {}):

	print(' ... ... Computing sequences similarity matrix')
	result_handle = open(blast_results)
	blast_records = NCBIXML.parse(result_handle)
	
	all_queries = sorted(list(set([record.query.split('|')[0] for record in blast_records])))
	distance_matrix = [[default_base for query in all_queries] for query in all_queries]

	result_handle = open(blast_results)
	blast_records = NCBIXML.parse(result_handle)
	
	for record in blast_records:
		query_name = record.query.split('|')[0]
		query_index = all_queries.index(query_name)

		for alignment in record.alignments:
			target_name = alignment.title.split('|')[2].split()[1]
			target_index = all_queries.index(target_name)

			if query_name != target_name:
				query_intervals = []
				sbjct_intervals = []

				for hsp in alignment.hsps:
					query_intervals.append(np.array([hsp.query_start, hsp.query_end]))
					sbjct_intervals.append(np.array([hsp.sbjct_start, hsp.sbjct_end]))

				query_intervals  = merge_intervals(query_intervals)
				target_intervals = merge_intervals(sbjct_intervals)

				if query_name in sequences_lengths and target_name in sequences_lengths:
					query_length  = sequences_lengths[query_name]
					target_lenght = sequences_lengths[target_name]

					query_coverage  = sum([i[-1]-i[0] for i in query_intervals])*100.0/float(query_length)
					target_coverage = sum([i[-1]-i[0] for i in target_intervals])*100.0/float(target_lenght)

					if query_coverage >= min_coverage and target_coverage >= min_coverage:
						distance_matrix[query_index][target_index] = 0
						distance_matrix[target_index][query_index] = 0

		distance_matrix[query_index][query_index] = 0

	return np.array(distance_matrix), all_queries

# 3.2. Routines for when the method is based on MMseqs

#@Timer("get_queries_labels", text=lambda secs: f"Time to get queries labels: {format_timespan(secs)}")
def get_queries_labels(seq_fasta):

	all_queries = []
	with open(seq_fasta, 'r') as infasta:
		for line in infasta:
			if line.startswith('>'):
				query = line.split('|')[0].strip('>')
				all_queries.append(query)
	return all_queries

#@Timer("run_mmseqs_for_flanking_sequences", text=lambda secs: f"Time to run mmseqs for flanking sequences: {format_timespan(secs)}")
def run_mmseqs_for_flanking_sequences(seq_fasta, num_threads = None, max_evalue = None, num_iterations = None, mmseqs = None, tmp_folder = None, sensitivity=7.5, min_coverage=None):

	print(' ... ... Running MMseqs')
	mmseqs_outfile = '{}/{}_{}.mmseqs'.format(tmp_folder, seq_fasta.split('/')[-1][:-6], max_evalue)

	if not os.path.isfile(mmseqs_outfile):
		run_mmseqs = sp.Popen(['mmseqs', 'easy-search', seq_fasta, seq_fasta, mmseqs_outfile, tmp_folder, '-e', str(max_evalue), '-s', str(sensitivity), '-c', str(min_coverage), '--num-iterations', str(num_iterations), '--threads', str(num_threads), '--format-output', 'query,target,evalue'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		stdout, stderr = run_mmseqs.communicate()

		if len(stderr) > 0:
			try:
				run_mmseqs = sp.Popen(['{}mmseqs'.format(mmseqs), 'easy-search', seq_fasta, seq_fasta, mmseqs_outfile, tmp_folder, '-e', str(max_evalue), '-s', str(sensitivity), '-c', str(min_coverage), '--num-iterations', str(num_iterations), '--threads', str(num_threads), '--format-output', 'query,target,evalue'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
				stdout, stderr = run_mmseqs.communicate()
			except:
				print("\n--> ERROR:  There's no MMseqs installation")
				exit()

	return mmseqs_outfile

#@Timer("extract_distance_matrix_from_mmseqs_output", text=lambda secs: f"Time to extract distance matrix from mmseqs output: {format_timespan(secs)}")
def extract_distance_matrix_from_mmseqs_output(mmseqs_results, all_queries, default_base = None):

	print(' ... ... Computing sequences similarity matrix')
	distance_matrix = [[default_base if i!=j else 0 for i in all_queries] for j in all_queries]

	all_queries = {query: i for i, query in enumerate(all_queries)}

	with open(mmseqs_results, 'r') as mmseqs_records:
		for hsp in mmseqs_records:
			hsp = hsp.split()
			if len(hsp) > 0:
				query = hsp[0].split('|')[0]
				query_index = all_queries[query]

				target = hsp[1].split('|')[0]
				if target != query:
					target_index = all_queries[target]
					distance_matrix[query_index][target_index] = 0
					distance_matrix[target_index][query_index] = 0

	return np.array(distance_matrix)

# 3.3. Wrapping routines

#@Timer("compute_all_agains_all_distance_matrix", text=lambda secs: f"Time to compute all agains all distance matrix: {format_timespan(secs)}")
def compute_all_agains_all_distance_matrix(in_syntenies, out_label = None, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, blast = None, mmseqs = None, default_base = None, tmp_folder = None, mode = 'flanking_sequences', method = 'blast', min_coverage=None):
	
	out_dir = '{}/{}_all_against_all_searches'.format(os.getcwd(), out_label)
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	flanking_fasta, sequences_len = write_flanking_sequences_to_fasta(in_syntenies, out_dir, out_label, mode = mode)

	if method == 'psiblast':

		sequences_database = make_blast_database_from_fasta(flanking_fasta, blast = blast)
		blast_results = run_blast_for_flanking_sequences(flanking_fasta, sequences_database, num_threads = num_threads, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, blast = blast, tmp_folder = tmp_folder)
		distance_matrix, queries_labels = extract_distance_matrix_from_blast_output(blast_results, default_base = default_base, sequences_lengths = sequences_len, min_coverage = min_coverage)
	
	elif method == 'mmseqs':
		queries_labels = get_queries_labels(flanking_fasta)
		mmseqs_results = run_mmseqs_for_flanking_sequences(flanking_fasta, num_threads = num_threads, max_evalue = max_evalue, num_iterations = num_iterations, min_coverage = min_coverage/100, mmseqs = mmseqs, tmp_folder = tmp_folder)
		distance_matrix = extract_distance_matrix_from_mmseqs_output(mmseqs_results, queries_labels, default_base = default_base)

	return distance_matrix, queries_labels

#@Timer("get_protein_families_summary", text=lambda secs: f"Time to get protein families summary: {format_timespan(secs)}")
def get_protein_families_summary(in_syntenies, write_to_file = True, out_label = None):

	families = {}

	for target in in_syntenies:
		for i, family in enumerate(in_syntenies[target]['flanking_genes']['families']):
			curr_ncbi_code = in_syntenies[target]['flanking_genes']['ncbi_codes'][i]

			if family not in families:
				families[family] = {'name': [], 'members': [], 'all_names': []}

			if curr_ncbi_code not in families[family]['members']:

				if family != 0:
					families[family]['all_names'].append(in_syntenies[target]['flanking_genes']['names'][i])
				else:
					families[family]['all_names'] = ['Non-conserved']

				families[family]['members'].append(curr_ncbi_code)

	for family in families:
		if len(set(families[family]['all_names'])) > 1 and 'hypothetical protein' in set(families[family]['all_names']):
			families[family]['name'] = [name for name in families[family]['all_names'] if name != 'hypothetical protein']
		else:
			families[family]['name'] = families[family]['all_names']
		
		try:
			families[family]['name'] = statistics.mode(families[family]['name'])
		except:
			families[family]['name'] = families[family]['name'][0]

	if 10000 in families:
		n_pseudogenes = len(families[10000]['members'])
	else:
		n_pseudogenes = 0

	if 0 in families:
		n_nonconserved = len(families[0]['members'])
	else:
		n_nonconserved = 0
		
	print(' ... Found {} conserved protein families, {} pseudogenes and {} non-conserved protein coding regions'.format(len([i for i in families if i != 0 and i != 10000]), n_pseudogenes, n_nonconserved))

	if write_to_file:
		out_file = '{}_protein_families_summary.txt'.format(out_label)
		with open(out_file, 'w') as outf:
			for family in sorted(list(families.keys())):
				if families[family]['name'] != 'Non-conserved':
					outf.write('\n ### Family: {} -> {}\n\n'.format(family, families[family]['name']))
					for i, member in enumerate(families[family]['members']):
						outf.write('	 {}\t{}\n'.format(member, families[family]['all_names'][i]))
			
	return families

# 4. Routines to annotate functions and find structures for protein families

#@Timer("get_mappings_through_uniprot", text=lambda secs: f"Time to get mappings through uniprot: {format_timespan(secs)}")
def get_mappings_through_uniprot(codes, from_database = 'RefSeq_Protein', to_database = 'UniProtKB', POLLING_INTERVAL = 10, API_URL = "https://rest.uniprot.org"):

	# The code below is copy-pasted from the help page of UniProt API (https://www.uniprot.org/help/id_mapping)


	def check_response(response):
		try:
			response.raise_for_status()
		except requests.HTTPError:
			print(response.json())
			raise


	def submit_id_mapping(from_db, to_db, ids):
		request = requests.post(
			f"{API_URL}/idmapping/run",
			data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
		)
		check_response(request)
		return request.json()["jobId"]


	def get_next_link(headers):
		re_next_link = re.compile(r'<(.+)>; rel="next"')
		if "Link" in headers:
			match = re_next_link.match(headers["Link"])
			if match:
				return match.group(1)


	def check_id_mapping_results_ready(job_id):
		while True:
			request = session.get(f"{API_URL}/idmapping/status/{job_id}")
			check_response(request)
			j = request.json()
			if "jobStatus" in j:
				if j["jobStatus"] == "RUNNING":
					print(f"Retrying in {POLLING_INTERVAL}s")
					time.sleep(POLLING_INTERVAL)
				else:
					raise Exception(j["jobStatus"])
			else:
				try:
					return bool(j["results"] or j["failedIds"])
				except:
					return True


	def get_batch(batch_response, file_format, compressed):
		batch_url = get_next_link(batch_response.headers)
		while batch_url:
			batch_response = session.get(batch_url)
			batch_response.raise_for_status()
			yield decode_results(batch_response, file_format, compressed)
			batch_url = get_next_link(batch_response.headers)


	def combine_batches(all_results, batch_results, file_format):
		if file_format == "json":
			for key in ("results", "failedIds"):
				if key in batch_results and batch_results[key]:
					all_results[key] += batch_results[key]
		elif file_format == "tsv":
			return all_results + batch_results[1:]
		else:
			return all_results + batch_results
		return all_results


	def get_id_mapping_results_link(job_id):
		url = f"{API_URL}/idmapping/details/{job_id}"
		request = session.get(url)
		check_response(request)
		return request.json()["redirectURL"]


	def decode_results(response, file_format, compressed):
		if compressed:
			decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
			if file_format == "json":
				j = json.loads(decompressed.decode("utf-8"))
				return j
			elif file_format == "tsv":
				return [line for line in decompressed.decode("utf-8").split("\n") if line]
			elif file_format == "xlsx":
				return [decompressed]
			elif file_format == "xml":
				return [decompressed.decode("utf-8")]
			else:
				return decompressed.decode("utf-8")
		elif file_format == "json":
			return response.json()
		elif file_format == "tsv":
			return [line for line in response.text.split("\n") if line]
		elif file_format == "xlsx":
			return [response.content]
		elif file_format == "xml":
			return [response.text]
		return response.text


	def get_xml_namespace(element):
		m = re.match(r"\{(.*)\}", element.tag)
		return m.groups()[0] if m else ""


	def merge_xml_results(xml_results):
		merged_root = ElementTree.fromstring(xml_results[0])
		for result in xml_results[1:]:
			root = ElementTree.fromstring(result)
			for child in root.findall("{http://uniprot.org/uniprot}entry"):
				merged_root.insert(-1, child)
		ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
		return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


	def print_progress_batches(batch_index, size, total):
		n_fetched = min((batch_index + 1) * size, total)
		print(f"Fetched: {n_fetched} / {total}")


	def get_id_mapping_results_search(url):
		parsed = urlparse(url)
		query = parse_qs(parsed.query)
		file_format = query["format"][0] if "format" in query else "json"
		if "size" in query:
			size = int(query["size"][0])
		else:
			size = 500
			query["size"] = size
		compressed = (
			query["compressed"][0].lower() == "true" if "compressed" in query else False
		)
		parsed = parsed._replace(query=urlencode(query, doseq=True))
		url = parsed.geturl()
		request = session.get(url)
		check_response(request)
		results = decode_results(request, file_format, compressed)
		total = int(request.headers["x-total-results"])
		print_progress_batches(0, size, total)
		for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
			results = combine_batches(results, batch, file_format)
			print_progress_batches(i, size, total)
		if file_format == "xml":
			return merge_xml_results(results)
		return results


	def get_id_mapping_results_stream(url):
		if "/stream/" not in url:
			url = url.replace("/results/", "/results/stream/")
		request = session.get(url)
		check_response(request)
		parsed = urlparse(url)
		query = parse_qs(parsed.query)
		file_format = query["format"][0] if "format" in query else "json"
		compressed = (
			query["compressed"][0].lower() == "true" if "compressed" in query else False
		)
		return decode_results(request, file_format, compressed)


	job_id = submit_id_mapping(
		from_db=from_database, to_db=to_database, ids=codes
	)
	if check_id_mapping_results_ready(job_id):
		link = get_id_mapping_results_link(job_id)
		results = get_id_mapping_results_search(link)["results"]

	return results


#@Timer("map_codes_to_uniprot", text=lambda secs: f"Time to map codes to uniprot: {format_timespan(secs)}")
def map_codes_to_uniprot(codes, uniprot_codes = {}, from_database = 'RefSeq_Protein'):
	
	if uniprot_codes == {}:
		uniprot_codes = {code: 'nan' for code in codes}
		
	results = get_mappings_through_uniprot(codes, from_database = from_database, to_database = 'UniProtKB')
	
	if results is not None:
		for result in results:
			in_code = result['from']
			uniprot_code = result['to']['primaryAccession']
			uniprot_codes[in_code] = uniprot_code 

	unmapped_codes = [code for code in uniprot_codes if uniprot_codes[code] == 'nan']
	if len(unmapped_codes) > 0 and from_database == 'RefSeq_Protein':
		uniprot_codes =  map_codes_to_uniprot(unmapped_codes, uniprot_codes = uniprot_codes, from_database = 'EMBL-GenBank-DDBJ_CDS')  

	upi_codes = [uniprot_codes[code] for code in uniprot_codes if 'UPI' in uniprot_codes[code]]
	if len(upi_codes) > 0:
		upi_codes =  map_codes_to_uniprot(upi_codes, from_database = 'UniParc')
		for code in uniprot_codes:
			if uniprot_codes[code] in upi_codes:
				uniprot_codes[code] = upi_codes[uniprot_codes[code]]

	return uniprot_codes

#@Timer("find_uniprot_in_swiss_model_repository", text=lambda secs: f"Time to find uniprot in swiss model repository: {format_timespan(secs)}")
def find_uniprot_in_swiss_model_repository(uniprot_code):

	link = 'nan'

	if uniprot_code != 'nan':

		link = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
		json_link = json_link = '{}.json'.format(link)

		try:
			swissmodel_req = requests.get(json_link)
			if swissmodel_req.ok:
				swiss_repository_data = swissmodel_req.text
				swiss_repository_data = json.loads(swiss_repository_data)

				if len(swiss_repository_data['result']) == 0 or len(swiss_repository_data['result']['uniprot_entries']) == 0:
					link = 'nan*'
				elif len(swiss_repository_data['result']['structures']) == 0:
					link = 'nan'
			else:
				link = 'nan*'
		except:
			link = 'nan*'

	return link

#@Timer("find_uniprot_in_alphafold_database", text=lambda secs: f"Time to find uniprot in alphafold database: {format_timespan(secs)}")
def find_uniprot_in_alphafold_database(uniprot_code):

	link = 'nan'

	if uniprot_code != 'nan':

		link = 'https://alphafold.ebi.ac.uk/entry/{}'.format(uniprot_code)
		protein_link = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v3.pdb'.format(uniprot_code)
		try:
			afdb_req = requests.get(protein_link)
			if not afdb_req.ok:
				link = 'nan*'
			else:
				link = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
		except:
			link = 'nan*'

	return link

#@Timer("get_uniprot_annotations", text=lambda secs: f"Time to get uniprot annotations: {format_timespan(secs)}")
def get_uniprot_annotations(uniprot_code, previous_annotations = ''):

	uniprot_annotations = 'nan'

	if uniprot_code != 'nan':
		uniprot_accession = uniprot_code.split('_')[0]
		uniprot_link = 'https://www.ebi.ac.uk/proteins/api/proteins/{}'.format(uniprot_accession)

		try:
			uniprot_req = requests.get(uniprot_link, headers={ "Accept" : "application/json"})

			if uniprot_req.ok:
				uniprot_data = uniprot_req.text
				uniprot_data = json.loads(uniprot_data)
				uniprot_annotations = parse_uniprot_data(uniprot_data, previous_annotations = previous_annotations)
		except:
			uniprot_annotations = 'nan'
	
	return uniprot_annotations

#@Timer("parse_uniprot_data", text=lambda secs: f"Time to parse uniprot data: {format_timespan(secs)}")
def parse_uniprot_data(uniprot_data, previous_annotations = ''):

	if previous_annotations == '':
		uniprot_annotations = {'TM_topology': '', 'GO_terms': [], 'Keywords': [], 'Function_description': ''}
	else:
		uniprot_annotations = previous_annotations

	if 'features' in uniprot_data:
		for feature in uniprot_data['features']:
			if feature['type'] == 'TRANSMEM' and uniprot_annotations['TM_topology'] == '':
				tm_topology = feature['description']

				if 'Helical' in tm_topology:
					tm_topology = 'alpha-helical'
				elif 'Beta' in tm_topology:
					tm_topology = 'beta-stranded'
				else:
					tm_topology = 'transmembrane'

				uniprot_annotations['TM_topology'] = tm_topology

	if 'dbReferences' in uniprot_data:
		for dbReference in uniprot_data['dbReferences']:
			if dbReference['type'] == 'GO':
				go_term = dbReference['properties']['term']
				if go_term not in uniprot_annotations['GO_terms']:
					uniprot_annotations['GO_terms'].append(go_term)

	if 'keywords' in uniprot_data:
		for keyword in uniprot_data['keywords']:
			keyword_term = keyword['value']
			if keyword_term not in uniprot_annotations['Keywords'] and keyword_term != 'Reference proteome':
				uniprot_annotations['Keywords'].append(keyword_term)

	if 'comments' in uniprot_data:
		for comment in uniprot_data['comments']:
			if comment['type'] == 'FUNCTION' and uniprot_annotations['Function_description'] == '':
				uniprot_annotations['Function_description'] = comment['text'][0]['value']

	return uniprot_annotations

# 5. Routines to find operon types

#@Timer("compute_operons_distance_matrix", text=lambda secs: f"Time to compute operons distance matrix: {format_timespan(secs)}")
def compute_operons_distance_matrix(in_syntenies, label = None):

	distance_matrix = [[0 for target in in_syntenies] for target in in_syntenies]
	sorted_ncbi_codes = sorted(list(in_syntenies.keys()))

	for i, operon_i in enumerate(sorted_ncbi_codes):
		for j, operon_j in enumerate(sorted_ncbi_codes):
			if i < j:
				vector_i = in_syntenies[operon_i]['flanking_genes']['families']
				reference_family_i = in_syntenies[operon_i]['target_family']
				
				vector_j = in_syntenies[operon_j]['flanking_genes']['families']
				reference_family_j = in_syntenies[operon_j]['target_family']

				if len(vector_i) >= len(vector_j):
					reference_vector = vector_i
					curr_vector = vector_j
					index_reference_family_reference_vector = list(reference_vector).index(reference_family_i)
					index_reference_family_curr_vector = list(curr_vector).index(reference_family_j)
				else:
					reference_vector = vector_j
					curr_vector = vector_i
					index_reference_family_reference_vector = list(reference_vector).index(reference_family_j)
					index_reference_family_curr_vector = list(curr_vector).index(reference_family_i)

				number_fillings = 0
				if len(curr_vector) < len(reference_vector):
					if index_reference_family_curr_vector < index_reference_family_reference_vector:
						for a in range(index_reference_family_reference_vector - index_reference_family_curr_vector):
							curr_vector = np.insert(curr_vector, 0, -1)
							number_fillings += 1
					
					if len(curr_vector) < len(reference_vector):
						for b in range(len(reference_vector) - len(curr_vector)):
							curr_vector = np.append(curr_vector, -1)
							number_fillings += 1


				families_present = set(np.concatenate([reference_vector,curr_vector]))
				overlapping_families = set(reference_vector).intersection(set(curr_vector))
				dist = 1-len(overlapping_families)/len(families_present)

				distance_matrix[i][j] = dist
				distance_matrix[j][i] = dist
	
	return np.array(distance_matrix), sorted_ncbi_codes
	
@Timer("get_operon_types_summary", text=lambda secs: f"Time to get operon types summary: {format_timespan(secs)}")
def get_operon_types_summary(in_syntenies, label = None, write_to_file = True):

	operon_types = {}
	advanced = False
	if 'operon_PaCMAP' in in_syntenies[list(in_syntenies.keys())[0]]:
		advanced = True

	for target in in_syntenies:
		curr_operon_type = in_syntenies[target]['operon_type']
		if curr_operon_type not in operon_types:
			operon_types[curr_operon_type] = {'target_members': [], 'operon_protein_families_structure': []}
			if advanced:
				operon_types[curr_operon_type]['operon_PaCMAP'] = []
				operon_types[curr_operon_type]['operon_filtered_PaCMAP'] = []


		operon_types[curr_operon_type]['target_members'].append(target)
		operon_types[curr_operon_type]['operon_protein_families_structure'].append(in_syntenies[target]['flanking_genes']['families'])
		if advanced:
			operon_types[curr_operon_type]['operon_PaCMAP'].append(in_syntenies[target]['operon_PaCMAP'])
			operon_types[curr_operon_type]['operon_filtered_PaCMAP'].append(in_syntenies[target]['operon_filtered_PaCMAP'])

	if advanced:
		for curr_operon_type in operon_types:
			centroid_coords = np.mean(operon_types[curr_operon_type]['operon_filtered_PaCMAP'], axis=0)
			operon_types[curr_operon_type]['operon_centroid_PaCMAP'] = list(centroid_coords)

	print(' ... Found {} operon types (out of a total of {} input targets)'.format(len(operon_types), len(in_syntenies)))

	if write_to_file:
		out_file = '{}_operon_types_summary.txt'.format(label)
		with open(out_file, 'w') as outf:
			for operon in sorted(list(operon_types.keys())):
				outf.write('\n ### Operon: {}\n\n'.format(operon))
				for i, member in enumerate(operon_types[operon]['target_members']):
					outf.write('	 {}\t{}\n'.format(member, operon_types[operon]['operon_protein_families_structure'][i]))
						
	return operon_types

#@Timer("find_most_populated_operon_types", text=lambda secs: f"Time to find most populated operon types: {format_timespan(secs)}")
def find_most_populated_operon_types(operon_types_summary, nmax = None):

	operons_count_matrix = []

	for operon in operon_types_summary:
		operons_count_matrix.append([operon, len(operon_types_summary[operon]['target_members'])])

	operons_count_matrix = pd.DataFrame(operons_count_matrix)
	operons_count_matrix = operons_count_matrix.sort_values(by = [1, 0], ascending = [False, True])	
	operons_count_matrix = np.array(operons_count_matrix)
	
	if len(operons_count_matrix) > nmax:
		operons_count_matrix = operons_count_matrix[:nmax+1]

	selected_operons = {}
	most_populated_operon = ''
	for i, line in enumerate(operons_count_matrix):
		label = 'GC Type {:05d}'.format(line[0])
		if i == 0:
			most_populated_operon = label
		
		selected_operons[label] = operon_types_summary[line[0]]

	print(' ... Selected {} operon/genomic_context types, with most populated corresponding to {}'.format(len(selected_operons), most_populated_operon))

	return selected_operons, most_populated_operon


#@Timer("get_family_presence_matrix", text=lambda secs: f"Time to get family presence matrix: {format_timespan(secs)}")
def get_family_presence_matrix(in_syntenies, protein_families_summary, clean = True, min_freq = 2, max_freq = 20):

	sorted_ncbi_codes = sorted(list(in_syntenies.keys()))
	sorted_families   = [i for i in sorted(list(protein_families_summary.keys())) if (i>=0 and i<10000)]

	# select only the protein families that are not very frequenct but also not very rare
	if clean and len(sorted_families) > 10:
		families_frequency = [len(protein_families_summary[family]['members']) for family in sorted_families]
		families_frequency = [i*100/len(in_syntenies) for i in families_frequency]
		sorted_families	= [family for i, family in enumerate(sorted_families) if families_frequency[i] <= max_freq and families_frequency[i] >= min_freq]

	presence_matrix = [[0 for i in sorted_families] for i in sorted_ncbi_codes]
	for i, target_i in enumerate(sorted_ncbi_codes):
		operon_i = in_syntenies[target_i]['flanking_genes']['families']
		for family in operon_i:
			if family in sorted_families:
				presence_matrix[i][sorted_families.index(family)] += 1

	return np.array(presence_matrix), sorted_ncbi_codes, sorted_families

#@Timer("calculate_start_eps", text=lambda secs: f"Time to calculate start eps: {format_timespan(secs)}")
def calculate_start_eps(coordinates):

	distances = []
	for i, vector_i in enumerate(coordinates):
		for j, vector_j in enumerate(coordinates):
			if j>i:
				dist = np.linalg.norm(np.array(vector_i) - np.array(vector_j))
				distances.append(dist)

	distances = np.array(distances)
	mixture = GaussianMixture(n_components=2).fit(distances.reshape(-1,1))
	means = mixture.means_.flatten()
	sds = np.sqrt(mixture.covariances_).flatten()

	mean = min(means)
	sd = sds[list(means).index(mean)]

	eps = mean - sd

	return round(eps, 2)

#@Timer("find_operon_clusters_with_PaCMAP", text=lambda secs: f"Time to find operon clusters with PaCMAP: {format_timespan(secs)}")
def find_operon_clusters_with_PaCMAP(in_syntenies, protein_families_summary, clean = True, coordinates_only = False, min_freq = 2, max_freq = 20, iteration = 0, max_eps_cost = 1.3, eps_step = 0.02):

	presence_matrix, sorted_ncbi_codes, selected_families = get_family_presence_matrix(in_syntenies, protein_families_summary, clean = clean, min_freq = min_freq, max_freq = max_freq)

	paCMAP_embedding = pacmap.PaCMAP(n_components = 2)
	paCMAP_coordinat = paCMAP_embedding.fit_transform(presence_matrix)

	if coordinates_only:
		return paCMAP_coordinat, sorted_ncbi_codes

	else:
		# embed into n-D paCMAP space
		n_dims = len(selected_families)
		if n_dims < 2:
			n_dims = 2

		paCMAP_embedding = pacmap.PaCMAP(n_components = n_dims)
		paCMAP_N_coordinat = paCMAP_embedding.fit_transform(presence_matrix)

		# find clusters in the paCMAP space
		# do this by selecting the best eps based on the number of clusters it creates compared to the number of operons 
		# that are not assigned a clusters (i.e., a given maximum cost)
		
		print(' ... Fine tuning EPS')

		eps = calculate_start_eps(paCMAP_N_coordinat)

		n_clusters = [0]
		n_singletons = [0]
		cost = 0
		while cost <= max_eps_cost and eps > 0:
			eps = eps - eps_step

			model = DBSCAN(eps = eps)
			model.fit(paCMAP_coordinat)
			clusters = model.fit_predict(paCMAP_coordinat)

			n = len(set(clusters))
			delta_n_clusters = n - n_clusters[-1]
			delta_singletons = list(clusters).count(-1) - n_singletons[-1]
			if delta_n_clusters > 0:
				cost = delta_singletons/delta_n_clusters
			else:
				cost = 0

		print(' ... ... EPS:  {}'.format(eps))
		print(' ... ... Cost: {}'.format(cost))
		print(' ... ... N:	{}'.format(n))


		return paCMAP_coordinat, clusters, sorted_ncbi_codes

# 6. Routines to make the genomic_context/operon block figures

#@Timer("define_family_colors", text=lambda secs: f"Time to define family colors: {format_timespan(secs)}")
def define_family_colors(families, reference_family, mode = 'matplotlib', cmap = 'rainbow', print_summary = False):

	colors = {}

	cmap = matplotlib.cm.get_cmap(cmap)
	norm = matplotlib.colors.Normalize(vmin=0, vmax=len(families))

	colours = [cmap(norm(i)) for i in range(len(families))]
	#random.shuffle(colours)

	if print_summary:
		print(" ... Setting the colors of the families identified")
		print(" ... ... Colormap: {}".format(cmap))
		print('\nLabel\tColor (RGBA)')

	for i, label in enumerate(sorted(families)):
		if label not in colors:
			colors[label] = {}

		if label == reference_family:			   # the reference gene
			colors[label]['Color (RGBA)'] = 'grey'
			colors[label]['Color (tuplet)'] = 'grey'
			colors[label]['Line color'] = 'black'
			colors[label]['Line style'] = '-'
		elif label == 0:							# a gene without any other homolog in the figure
			colors[label]['Color (RGBA)'] = 'lightgrey'
			colors[label]['Color (tuplet)'] = 'lightgrey'
			colors[label]['Line color'] = 'lightgrey'
			colors[label]['Line style'] = '-'
		elif type(label) == int and label == 10000:  # a pseudogene
			colors[label]['Color (RGBA)'] = 'white'
			colors[label]['Color (tuplet)'] = 'white'
			colors[label]['Line color'] = 'grey'
			colors[label]['Line style'] = ':'
		else:									   # real genes that occur many times in the figure
			if mode == 'matplotlib':
				colors[label]['Color (RGBA)'] = [int(255*j) for j in colours[i]]
				colors[label]['Color (tuplet)'] = colours[i]
				colors[label]['Line color'] = 'black'
				colors[label]['Line style'] = '-'
			elif mode == 'bokeh':
				colors[label]['Color (RGBA)'] = RGB(int(255*list(colours[i])[0]), int(255*list(colours[i])[1]), int(255*list(colours[i])[2]))
				colors[label]['Color (tuplet)'] = colours[i]
				colors[label]['Line color'] = 'black'
				colors[label]['Line style'] = '-'

		if print_summary:
			print('{}\t{}'.format(label, colors[label]['Color (RGBA)']))

	return colors


#@Timer("draw_genomic_context", text=lambda secs: f"Time to draw genomic context: {format_timespan(secs)}")
def draw_genomic_context(operons, all_syntenies, family_colors, reference_family, label = None, out_format = None):

	curr_y_level = len(operons.keys())
	all_xs = []
	all_populations = []
	yticklabels = []

	plt.clf()
	if len(operons) == 1:
		fig, ax = plt.subplots(1, 2, figsize=(20, 2), gridspec_kw={'width_ratios': [4, 1]})
	elif len(operons) < 5:
		fig, ax = plt.subplots(1, 2, figsize=(20, len(operons)), gridspec_kw={'width_ratios': [4, 1]})
	else:
		fig, ax = plt.subplots(1, 2, figsize=(20, int(len(operons)/1.5)), gridspec_kw={'width_ratios': [4, 1]})

	for operon in sorted(list(operons.keys())):
		current_target = operons[operon]['target_members'][0]
		current_genomic_context_block = all_syntenies[current_target]['flanking_genes']
		current_species = all_syntenies[current_target]['species']
		current_reference_family = all_syntenies[current_target]['target_family']
		operon_population = len(operons[operon]['target_members'])*100/len(all_syntenies)

		for i, flanking_gene in enumerate(current_genomic_context_block['ncbi_codes']):
			family = current_genomic_context_block['families'][i]
			gene_dx = current_genomic_context_block['relative_ends'][i] - current_genomic_context_block['relative_starts'][i]+1
			gene_direction = current_genomic_context_block['directions'][i]

			if gene_direction == '-':
				gene_x_tail = current_genomic_context_block['relative_ends'][i]
				gene_dx = gene_dx*(-1)
			else:
				gene_x_tail = current_genomic_context_block['relative_starts'][i]

			# make genomic_context side			
			if family == 0:
				zorder = family
				facecolor = family_colors[family]['Color (tuplet)']
				edgecolor = family_colors[family]['Line color']
				linestyle = family_colors[family]['Line style']
			else:
				zorder = len(current_genomic_context_block['ncbi_codes']) - i + 1
				facecolor = family_colors[family]['Color (tuplet)']
				edgecolor = family_colors[family]['Line color']
				linestyle = family_colors[family]['Line style']

			ax[0].arrow(gene_x_tail, curr_y_level, gene_dx, 0, width=0.5, head_width=0.5, length_includes_head = True, head_length = 100, zorder = zorder, facecolor = facecolor, edgecolor = edgecolor, linestyle = linestyle)
				
			text_x = gene_x_tail + (gene_dx/2)
			if family != 0 and family != reference_family and family < 10000:
				ax[0].text(text_x, curr_y_level+0.3, str(family), horizontalalignment='center')

			# make histogram side
			ax[1].arrow(0, curr_y_level, operon_population, 0, width=0.5, head_width=0.5, length_includes_head = True, head_length = 0, facecolor = 'black', edgecolor = 'black')
			ax[1].text(operon_population+2.5, curr_y_level, '{}'.format(round(operon_population, 1)), verticalalignment = 'center')
	
			all_xs.append(current_genomic_context_block['relative_ends'][i])
			all_xs.append(current_genomic_context_block['relative_starts'][i])
		
		all_populations.append(operon_population)
		curr_y_level -= 1
		yticklabels.append('{}: {} ({})'.format(operon, current_target, current_species))
		
	yticklabels.append('')
	yticklabels.reverse()

	# format genomic_context side
	ax[0].set_xlim(min(all_xs)-100, max(all_xs)+100)
	ax[0].set_ylim(0, len(operons.keys())+1)

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
	ax[1].set_ylim(0, len(operons.keys())+1)
	ax[1].spines['right'].set_visible(False)
	ax[1].tick_params(axis='y', which='both', left = False, right = False, labelleft=False, labelright=False)
	ax[1].set_xlabel('Operon frequency (%)')
	ax[1].set_title('{}% total operons represented'.format(round(sum(all_populations), 1)))

	try:
		plt.tight_layout()
	except:
		pass
	plt.savefig('{}_genomic_context.{}'.format(label, out_format), format = out_format)
	plt.close('all')

#@Timer("draw_genomic_context_legend", text=lambda secs: f"Time to draw genomic context legend: {format_timespan(secs)}")
def draw_genomic_context_legend(families_summary, family_colors, reference_family, label = None, out_format = None):

	curr_y_level = len(family_colors.keys())
	x_tail = 0
	dx = 5

	plt.clf()
	fig, ax = plt.subplots(figsize=(10, int(len(family_colors)/1.5)))
	for family in sorted(list(family_colors.keys())):
		plt.arrow(x_tail, curr_y_level, dx, 0, width=0.5, head_width=0.5, length_includes_head = True, head_length = 0.5, facecolor = family_colors[family]['Color (tuplet)'], edgecolor = family_colors[family]['Line color'], linestyle = family_colors[family]['Line style'])
		if family == 0:
			plt.text(dx + 2, curr_y_level, 'Non-conserved gene')
		elif family == reference_family:
			plt.text(dx + 2, curr_y_level, 'Target protein: {}'.format(families_summary[family]['name']) )
		elif family == 10000:
			plt.text(dx + 2, curr_y_level, 'Pseudogene')
		elif family in families_summary:
			plt.text(2.25, curr_y_level+0.3, str(family), horizontalalignment='center')
			plt.text(dx + 2, curr_y_level, families_summary[family]['name'])

		curr_y_level -= 1

	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	plt.yticks(fontsize = 10)
	plt.tick_params(axis='y', which='both', left = False, right = False, labelleft=False, labelright=False)
	plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
	
	plt.xlim(0, 50)
	plt.ylim(0, len(family_colors.keys())+1)
	
	plt.savefig('{}_genomic_context_legend.{}'.format(label, out_format), format = out_format)
	plt.close('all')

# 7. Routines to make taxonomy distribution figures

#@Timer("merge_taxonomy_dictionaries", text=lambda secs: f"Time to merge taxonomy dictionaries: {format_timespan(secs)}")
def merge_taxonomy_dictionaries(taxonomy, dic):

	for superkingdom in dic.keys():
		if superkingdom not in taxonomy.keys():
			taxonomy[superkingdom] = dic[superkingdom]
		else:
			for phylum in dic[superkingdom].keys():
				if phylum not in taxonomy[superkingdom].keys():
					taxonomy[superkingdom][phylum] = dic[superkingdom][phylum]
				else:
					for taxclass in dic[superkingdom][phylum].keys():
						if taxclass not in taxonomy[superkingdom][phylum].keys():
							taxonomy[superkingdom][phylum][taxclass] = dic[superkingdom][phylum][taxclass]
						else:
							for order in dic[superkingdom][phylum][taxclass].keys():
								if order not in taxonomy[superkingdom][phylum][taxclass].keys():
									taxonomy[superkingdom][phylum][taxclass][order] = dic[superkingdom][phylum][taxclass][order]
								else:
									for genus in dic[superkingdom][phylum][taxclass][order].keys():
										if genus not in taxonomy[superkingdom][phylum][taxclass][order].keys():
											taxonomy[superkingdom][phylum][taxclass][order][genus] = dic[superkingdom][phylum][taxclass][order][genus]
										else:
											for species in dic[superkingdom][phylum][taxclass][order][genus].keys():
												if species not in taxonomy[superkingdom][phylum][taxclass][order][genus].keys():
													taxonomy[superkingdom][phylum][taxclass][order][genus][species] = dic[superkingdom][phylum][taxclass][order][genus][species]
												else:
													for key in taxonomy[superkingdom][phylum][taxclass][order][genus][species].keys():
														taxonomy[superkingdom][phylum][taxclass][order][genus][species][key] += dic[superkingdom][phylum][taxclass][order][genus][species][key]
														taxonomy[superkingdom][phylum][taxclass][order][genus][species][key] = list(set(taxonomy[superkingdom][phylum][taxclass][order][genus][species][key]))

	return taxonomy

#@Timer("map_taxonomy_to_targets", text=lambda secs: f"Time to map taxonomy to targets: {format_timespan(secs)}")
def map_taxonomy_to_targets(in_syntenies, mode = 'taxonomy', threads = 1):

	# Prepare all parallel jobs
	separated_jobs = chunk_list(list(in_syntenies.keys()), threads)

	list_arguments = [i for i in zip(separated_jobs, [in_syntenies for job in separated_jobs], [mode for job in separated_jobs], range(threads))]

	pool = ThreadPool(threads)
	results = pool.imap_unordered(map_taxonomy, list_arguments)

	taxonomy = {}
	for dic in results:
		taxonomy = merge_taxonomy_dictionaries(taxonomy, dic)

	return taxonomy

#@Timer("map_taxonomy", text=lambda secs: f"Time to map taxonomy: {format_timespan(secs)}")
def map_taxonomy(arguments):

	targets = arguments[0]
	in_syntenies = arguments[1]
	mode = arguments[2]
	thread_id = arguments[3]

	out_json = '{}_taxonomy.json'.format(thread_id)
	if os.path.isfile(out_json) and mode == 'taxonomy':
		taxonomy = json.load(open(out_json, 'r'))
	else:
		taxonomy = {}

		for curr_target in targets:
			ncbi_code = in_syntenies[curr_target]['assembly_id'][0]
			taxID = in_syntenies[curr_target]['flanking_genes']['taxID']
			species = in_syntenies[curr_target]['flanking_genes']['species']

			try:
				t1 = Timer('map_taxonomy',
						   text=lambda secs: f"Time to map taxonomy for {ncbi_code}: {format_timespan(secs)}")
				t1.start()
				superkingdom = 'na'
				phylum = 'na'
				taxclass = 'na'
				order = 'na'
				genus = 'na'

				if mode == 'taxonomy':
					if species.split()[0] == 'uncultured':
						genus = species.split()[1]
					else:
						genus = species.split()[0]
			   
					taxsearch = Entrez.efetch(id = taxID, db = "taxonomy", retmode = "xml")
					taxrecords = Entrez.parse(taxsearch)
					for taxrecord in taxrecords:
						taxrecord = taxrecord

					for level in taxrecord[u"LineageEx"]:
						if level[u"Rank"] == "superkingdom":
							superkingdom = level[u"ScientificName"]
						elif level[u"Rank"] == "phylum":
							phylum = level[u"ScientificName"]
						elif level[u"Rank"] == "class":
							taxclass = level[u"ScientificName"]
						elif level[u"Rank"] == "order":
							order = level[u"ScientificName"]

					if phylum == 'na':
						phylum = '{}_na'.format(superkingdom)
					if taxclass == 'na':
						taxclass = '{}_na'.format(phylum)
					if order == 'na':
						order = '{}_na'.format(taxclass)

				if superkingdom not in taxonomy.keys():
					taxonomy[superkingdom] = {}
				   
				if phylum not in taxonomy[superkingdom].keys():
					taxonomy[superkingdom][phylum] = {}

				if taxclass not in taxonomy[superkingdom][phylum].keys():
					taxonomy[superkingdom][phylum][taxclass] = {}

				if order not in taxonomy[superkingdom][phylum][taxclass].keys():
					taxonomy[superkingdom][phylum][taxclass][order] = {}

				if genus not in taxonomy[superkingdom][phylum][taxclass][order].keys():
					taxonomy[superkingdom][phylum][taxclass][order][genus] = {}
				   
				if species not in taxonomy[superkingdom][phylum][taxclass][order][genus].keys():
					taxonomy[superkingdom][phylum][taxclass][order][genus][species] = {'ncbi_codes': [], 'target_members': []}

				taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members'].append(curr_target)
				taxonomy[superkingdom][phylum][taxclass][order][genus][species]['ncbi_codes'].append(ncbi_code)

				if mode == 'taxonomy':
					json.dump(taxonomy, open(out_json, 'w'), indent=4)

				t1.stop()

			except:
				print(' ... ... Not possible to find taxonomy for {} ({})'.format(curr_target, ncbi_code))

	return taxonomy

# 8. Routines to annotate transmembrane segments and signal peptides

#@Timer("run_TM_signal_peptide_annotation", text=lambda secs: f"Time to run TM signal peptide annotation in {annotation_TM_mode} mode: {format_timespan(secs)}")
def run_TM_signal_peptide_annotation(in_fasta, annotation_TM_mode = None):

	out_file = '{}_{}.out'.format(in_fasta[:-6], annotation_TM_mode)

	if not os.path.isfile(out_file):
		if annotation_TM_mode == 'phobius':
			try:
				run_phobius = sp.Popen(['phobius.pl', in_fasta, '-short'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
				stdout, stderr = run_phobius.communicate()
			except:
				stderr = " --> ERROR:  There's no phobius installation"
				stdout = ""

			if len(stderr) > 0 or len(stdout) == 0:
				print(stderr)
				print('	 Run phobius online and come back with the input file')
				pass
			else:
				with open(out_file, 'w') as outf:
					out_file.write(stdout)


		elif annotation_TM_mode == 'tmhmm':

			try:
				run_tmhmm = sp.Popen(['tmhmm', in_fasta, '-short'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
				stdout, stderr = run_tmhmm.communicate()
			except:
				stderr = " --> ERROR:  There's no tmhmm installation"
				stdout = ""

			if len(stderr) > 0 or len(stdout) == 0:
				print(stderr)
				print('	 Run tmhmm online and come back with the input file')
				pass
			else:
				with open(out_file, 'w') as outf:
					out_file.write(stdout)

		elif annotation_TM_mode == 'uniprot':
			
			ncbi_codes = list(parse_fasta_file(in_fasta).keys())
			uniprot_codes = map_codes_to_uniprot(ncbi_codes)

			with open(out_file, 'w') as outf:
				for ncbi_code in ncbi_codes:
					uniprot_code = uniprot_codes[ncbi_code]
					t1 = Timer('get_uniprot_annotation_for_each_peptide',
						   text=lambda secs: f"Time to get uniprot tm annotation for {ncbi_code}: {format_timespan(secs)}")
					t1.start()
					curr_uniprot_annotations = get_uniprot_annotations(uniprot_code)
					t1.stop()

					if curr_uniprot_annotations != 'nan':
						tm_annotation = curr_uniprot_annotations['TM_topology']
						if tm_annotation == '':
							tm_annotation = 'nan'
					else:
						tm_annotation = 'nan'
					outf.write('{}\t{}\n'.format(ncbi_code, tm_annotation))

		else:
			print(' ... ... ... {} mode not allowed')

	return out_file

#@Timer("parse_fasta_file", text=lambda secs: f"Time to parse fasta file: {format_timespan(secs)}")
def parse_fasta_file(infasta, sep = '|'):

	sequences = {}

	with open(infasta, 'r') as infst:
		for line in infst:
			if line.startswith('>'):
				ncbi_code = line.split('|')[0].replace('>','')
			else:
				sequences[ncbi_code] = line.strip()

	return sequences

#@Timer("parse_annotation_TM_file", text=lambda secs: f"Time to parse annotation TM file: {format_timespan(secs)}")
def parse_annotation_TM_file(annotation_TM_file, annotation_TM_mode):

	if annotation_TM_mode == 'phobius':
		annotations = parse_phobius_output(annotation_TM_file)

	elif annotation_TM_mode == 'tmhmm':
		annotations = parse_tmhmm_output(annotation_TM_file)

	elif annotation_TM_mode == 'uniprot':
		annotations = parse_tm_uniprot_output(annotation_TM_file)

	else:
		print(' ... ... ... {} mode not allowed')
		annotations = {}

	return annotations

#@Timer("parse_phobius_output", text=lambda secs: f"Time to parse phobius output: {format_timespan(secs)}")
def parse_phobius_output(annotation_TM_file):

	annotations = {}

	with open(annotation_TM_file, 'r') as outphobius:
		for line in outphobius:
			if 'PREDICTION' not in line:
				ncbi_code = line.split('|')[0]
				tm_pred = int(line.split()[1])
				sp_pred = line.split()[2]

				if tm_pred > 0:
					annotations[ncbi_code] = 'TM' # it has a transmembrane helix -> it is transmembrane
				elif sp_pred == 'Y':
					annotations[ncbi_code] = 'SP' # it has a signal peptide
				else:
					annotations[ncbi_code] = '' # does not have any special feature

	return annotations

#@Timer("parse_tmhmm_output", text=lambda secs: f"Time to parse tmhmm output: {format_timespan(secs)}")
def parse_tmhmm_output(annotation_TM_file):

	annotations = {}

	with open(annotation_TM_file, 'r') as outphobius:
		for line in outphobius:
			if 'PREDICTION' not in line:
				ncbi_code = line.split('|')[0]
				tm_pred = line.split('Topology=')[-1].strip()

				if tm_pred != 'o' and tm_pred != 'i':
					annotations[ncbi_code] = 'TM' # it has a transmembrane helix -> it is transmembrane
				else:
					annotations[ncbi_code] = '' # does not have any special feature (signal peptides are not identified)

	return annotations

#@Timer("parse_tm_uniprot_output", text=lambda secs: f"Time to parse_tm_uniprot_output: {format_timespan(secs)}")
def parse_tm_uniprot_output(annotation_TM_file):

	annotations = {}

	with open(annotation_TM_file, 'r') as outuniprot:
		for line in outuniprot:
			ncbi_code = line.split()[0]
			tm_pred = line.split()[-1].strip()

			if tm_pred != 'nan':
				annotations[ncbi_code] = 'TM' 
			else:
				annotations[ncbi_code] = '' 

	return annotations

#@Timer("add_TM_annotations_to_flanking_genes", text=lambda secs: f"Time to add_TM_annotations_to_flanking_genes: {format_timespan(secs)}")
def add_TM_annotations_to_flanking_genes(in_syntenies, protein_annotations):

	for curr_target in in_syntenies:
		flanking_genes = in_syntenies[curr_target]['flanking_genes']
		flanking_genes['TM_annotations'] = []

		for ncbi_code in flanking_genes['ncbi_codes']:
			if ncbi_code in protein_annotations:
				flanking_genes['TM_annotations'].append(protein_annotations[ncbi_code])
			else:
				flanking_genes['TM_annotations'].append('')

		in_syntenies[curr_target]['flanking_genes'] = flanking_genes

		json.dump(in_syntenies, open('all_syntenies.json', 'w'), indent = 4)
		
	return in_syntenies

# 9. Routines to make interactive output html file

# 9.1. Routines to draw most conserved features

#@Timer("find_most_common_genomic_context", text=lambda secs: f"Time to find most common genomic context: {format_timespan(secs)}")
def find_most_common_genomic_context(operons, all_syntenies, n_flanking5=None, n_flanking3=None):
	
	# will use only the complete genomic contexts and ignore the partial ones	
	operon_matrix = []
	for operon in operons:
		for curr_context in operons[operon]['operon_protein_families_structure']:
			if len(curr_context) == n_flanking5+n_flanking3+1:
				operon_matrix.append(curr_context)
	
	operon_matrix = np.array(operon_matrix).T
	
	most_common_context = {'selected_context': [],
						   'families_frequency': [],
						   'average_starts': [],
						   'average_ends': [],
						   'average_size': [],
						   'stdev_size': [],
						   'directions': [],
						   'tm_annotations': []}

	for i, column in enumerate(operon_matrix):
		occurence_count = Counter(column) 
		most_common_family = occurence_count.most_common(1)[0][0]

		most_common_context['selected_context'].append(most_common_family)
		most_common_context['families_frequency'].append(round(occurence_count.most_common(1)[0][1]*100/len(column), 1))

		all_starts_of_most_common = []
		all_ends_of_most_common = []
		all_orientations = []
		all_sizes = []
		all_tm_annotations = []

		for operon in operons:
			for j, curr_context in enumerate(operons[operon]['operon_protein_families_structure']):
				curr_target = operons[operon]['target_members'][j]
			
				if len(curr_context) == n_flanking5+n_flanking3+1:			
					if operons[operon]['operon_protein_families_structure'][j][i] == most_common_family:
						all_starts_of_most_common.append(all_syntenies[curr_target]['flanking_genes']['relative_starts'][i])
						all_ends_of_most_common.append(all_syntenies[curr_target]['flanking_genes']['relative_ends'][i])
						all_sizes.append((all_syntenies[curr_target]['flanking_genes']['relative_ends'][i] - all_syntenies[curr_target]['flanking_genes']['relative_starts'][i])/3)
						all_orientations.append(all_syntenies[curr_target]['flanking_genes']['directions'][i])
						
						if 'TM_annotations' in all_syntenies[curr_target]['flanking_genes']:
							all_tm_annotations.append(all_syntenies[curr_target]['flanking_genes']['TM_annotations'][i])
		
		most_common_context['average_starts'].append(int(statistics.median(all_starts_of_most_common)))
		most_common_context['average_ends'].append(int(statistics.median(all_ends_of_most_common)))
		most_common_context['average_size'].append(int(statistics.median(all_sizes)))
		most_common_context['stdev_size'].append(int(stats.median_abs_deviation(all_sizes)))
		
		try:
			most_common_context['directions'].append(statistics.mode(all_orientations))
		except:
			most_common_context['directions'].append('+')
		
		try:
			most_common_context['tm_annotations'].append(statistics.mode(all_tm_annotations))
		except:
			most_common_context['tm_annotations'].append('')
		
	return most_common_context

#@Timer("create_most_common_genomic_context_features", text=lambda secs: f"Time to create most common genomic context features: {format_timespan(secs)}")
def create_most_common_genomic_context_features(most_common_context, families_summary, reference_family, family_colors):

	data = {'xs': [],
			'ys': [],
			'edgecolor': [],
			'facecolor': [],
			'text_x': [],
			'text_y': [],
			'family': [],
			'tm_text_x': [],
			'tm_text_y': [],
			'tm_text': [],
			'tm_pred_text': [],
			'protein_name': [],
			'protein_size': [],
			'family_frequency': [],
			'transparency': [],
			'relative_start': [],
			'relative_end': [],
			'found_models': [],
			'model_links': []}
	
	for i, family in enumerate(most_common_context['selected_context']):
		gene_dx = most_common_context['average_ends'][i] - most_common_context['average_starts'][i]+1
		gene_direction = most_common_context['directions'][i]

		if gene_direction == '-':
			gene_x_tail = most_common_context['average_ends'][i]
			gene_dx = gene_dx*(-1)
			gene_x_head = gene_x_tail + gene_dx
			gene_x_head_start = gene_x_head+100
			text_x = gene_x_tail - (gene_x_tail-gene_x_head_start)/2
		else:
			gene_x_tail = most_common_context['average_starts'][i]
			gene_x_head = gene_x_tail + gene_dx
			gene_x_head_start = gene_x_head-100
			text_x = gene_x_tail + (gene_x_head_start-gene_x_tail)/2

		if family == 0:
			facecolor = family_colors[family]['Color (RGBA)']
			edgecolor = family_colors[family]['Line color']
			linestyle = family_colors[family]['Line style']
		else:
			facecolor = family_colors[family]['Color (RGBA)']
			edgecolor = family_colors[family]['Line color']
			linestyle = family_colors[family]['Line style'] 
		
		if family == 0:
			relative_start = 'n.a.'
			relative_end = 'n.a.'
			family_frequency = 'n.a.'
			protein_size = 'n.a.' 
			protein_name = 'n.a.'
			transparency = 0.2
			tm_annotation = ''
			tm_pred_text = 'n.a.'
		else:
			relative_start = format(gene_x_tail, ',d')
			relative_end = format(gene_x_head, ',d')
			family_frequency = '{}%'.format(most_common_context['families_frequency'][i])
			protein_size = r'{} ({})'.format(most_common_context['average_size'][i], most_common_context['stdev_size'][i])
			protein_name = families_summary[family]['name']
			transparency = most_common_context['families_frequency'][i]/100  
			
			if 'tm_annotations' in most_common_context:
				tm_annotation = most_common_context['tm_annotations'][i]
				if tm_annotation == 'TM':
					tm_pred_text = 'Yes'
				elif tm_annotation == 'SP':
					tm_pred_text = 'Contains signal peptide'
				else:
					tm_pred_text = 'No'
			else:
				tm_annotation = ''
				tm_pred_text = 'n.a.' 
				
		data['relative_start'].append(relative_start)
		data['relative_end'].append(relative_end)
		data['facecolor'].append(facecolor)
		data['edgecolor'].append(edgecolor)
		data['family_frequency'].append(family_frequency)
		data['protein_size'].append(protein_size)
		data['protein_name'].append(protein_name)
		data['xs'].append([gene_x_tail, gene_x_tail, gene_x_head_start, gene_x_head, gene_x_head_start])
		data['ys'].append([1-0.25, 1+0.25, 1+0.25, 1, 1-0.25])
		data['text_x'].append(text_x)
		data['text_y'].append(1+0.25)
		data['transparency'].append(transparency)
		
		data['tm_text'].append(tm_annotation)
		data['tm_text_x'].append(text_x)
		data['tm_text_y'].append(1)
		data['tm_pred_text'].append(tm_pred_text)

		if family != 0 and family != reference_family and family < 10000:
			data['family'].append(family)
		else:
			data['family'].append(str(''))
		
		if 'model_state' in families_summary[family]:
			model_state = families_summary[family]['model_state']

			if model_state == 'Model exists':
				model_state = 'Yes (click to view in Swiss-Model repository)'
			elif model_state == 'Model does not exist':
				model_state = 'No (click to model with Swiss-Model)'
			else:
				if family > 0 and family < 10000:
					model_state = 'Not possible to find'
				else:
					model_state = ''
			
			structure = families_summary[family]['structure']
			if structure == '':
				uniprot_code = families_summary[family]['uniprot_code']
				structure = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
			
		else:
			model_state = 'n.a.'
			structure = 'n.a.'
		
		data['found_models'].append(model_state)
		data['model_links'].append(structure)
		
	tooltips = [('Protein name', "@protein_name"),
				("Predicted membrane protein", "@tm_pred_text"),
				('Structural model found', '@found_models'),
				('Frequency in position', '@family_frequency'),
				('Median protein size', '@protein_size'),
				('Median starting position', '@relative_start'),
				('Median end position', '@relative_end')] 
	
	return tooltips, data

#@Timer("create_most_common_genomic_features_figure", text=lambda secs: f"Time to create most common genomic features figure: {format_timespan(secs)}")
def create_most_common_genomic_features_figure(operons, all_syntenies, families_summary, reference_family, family_colors, n_flanking5=None, n_flanking3=None):
	
	most_common_context = find_most_common_genomic_context(operons, all_syntenies, n_flanking5=n_flanking5, n_flanking3=n_flanking3)
			
	gc_tooltips, gc_data = create_most_common_genomic_context_features(most_common_context, families_summary, reference_family = reference_family, family_colors = family_colors)
		
	gc = figure(plot_width=2000, plot_height=200, y_range = [0, 4], title = 'Most conserved gene per position', toolbar_location="left")

	for i, xs in enumerate(gc_data['xs']):
		gc.patch(xs, gc_data['ys'][i], fill_color = gc_data['facecolor'][i], line_color = gc_data['edgecolor'][i], fill_alpha = gc_data['transparency'][i], line_alpha = gc_data['transparency'][i], line_width = 1)	
	
	gc.patches('xs', 'ys', fill_color = None, line_color = None, line_width = 0, source = gc_data, 
			  hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
			  selection_fill_color='facecolor', selection_line_color='edgecolor',
			  nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', nonselection_fill_alpha=0.2)

	gc.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", text_font_size = {'value': '6pt'}, source = gc_data)
	gc.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, source = gc_data)
	
	gc.yaxis.ticker = [1]
	gc.yaxis.major_label_overrides = {1: max([all_syntenies[target]['species'] for target in all_syntenies], key=len)}
#	 gc.yaxis.major_label_text_font_size = {'value': '8pt'}
	
	gc.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	gc.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	gc.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
	gc.yaxis.axis_line_width = 0

	# define general features
	gc.grid.visible = False
	gc.outline_line_width = 0
	
	# define xticks
	gc.xaxis.axis_label = "Position relative to target (bp)"
	
	gc.add_tools(HoverTool(tooltips=gc_tooltips))
	gc.add_tools(TapTool(callback = OpenURL(url='@model_links')))
	
	return gc

#@Timer("create_most_common_genomic_features_figure_2", text=lambda secs: f"Time to create most common genomic features figure 2: {format_timespan(secs)}")
def create_most_common_genomic_features_figure(operons, all_syntenies, families_summary, reference_family, family_colors, n_flanking5=None, n_flanking3=None):
	
	most_common_context = find_most_common_genomic_context(operons, all_syntenies, n_flanking5=n_flanking5, n_flanking3=n_flanking3)
			
	gc_tooltips, gc_data = create_most_common_genomic_context_features(most_common_context, families_summary, reference_family = reference_family, family_colors = family_colors)
		
	gc = figure(plot_width=2000, plot_height=200, y_range = [0, 4], title = 'Most conserved gene per position', toolbar_location="left")

	for i, xs in enumerate(gc_data['xs']):
		gc.patch(xs, gc_data['ys'][i], fill_color = gc_data['facecolor'][i], line_color = gc_data['edgecolor'][i], fill_alpha = gc_data['transparency'][i], line_alpha = gc_data['transparency'][i], line_width = 1)	
	
	gc.patches('xs', 'ys', fill_color = None, line_color = None, line_width = 0, source = gc_data, 
			  hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
			  selection_fill_color='facecolor', selection_line_color='edgecolor',
			  nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', nonselection_fill_alpha=0.2)

	gc.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", text_font_size = {'value': '6pt'}, source = gc_data)
	gc.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, source = gc_data)
	
	gc.yaxis.ticker = [1]
	gc.yaxis.major_label_overrides = {1: max([all_syntenies[target]['species'] for target in all_syntenies], key=len)}
	gc.yaxis.major_label_text_font_size = {'value': '8pt'}
	
	gc.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	gc.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	gc.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
	gc.yaxis.axis_line_width = 0

	# define general features
	gc.grid.visible = False
	gc.outline_line_width = 0
	
	# define xticks
	gc.xaxis.axis_label = "Position relative to target (bp)"
	
	gc.add_tools(HoverTool(tooltips=gc_tooltips))
	gc.add_tools(TapTool(callback = OpenURL(url='@model_links')))
	
	return gc

# 9.2. Routines to make the dendogram

#@Timer("get_taxonomy_distance_matrix", text=lambda secs: f"Time to get taxonomy distance matrix: {format_timespan(secs)}")
def get_taxonomy_distance_matrix(taxonomy, operons, input_targets = None, mode = 'taxonomy'):
	
	targets_taxonomy_vector = {}

	for superkingdom in taxonomy:
		for phylum in taxonomy[superkingdom]:
			for taxclass in taxonomy[superkingdom][phylum]:
				for order in taxonomy[superkingdom][phylum][taxclass]:
					for genus in taxonomy[superkingdom][phylum][taxclass][order]:
						for species in taxonomy[superkingdom][phylum][taxclass][order][genus]:
							for target in taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
								# find if this target is in the operons
								operon = [operon for operon in operons if target in operons[operon]['target_members']]
								if len(operon) > 0:
									targets_taxonomy_vector[target] = [superkingdom, phylum, taxclass, order, genus, species]

	if mode == 'taxonomy':

		labels = sorted(list(targets_taxonomy_vector.keys()))
		taxonomy_distance_matrix = [[0 for label in targets_taxonomy_vector] for label in targets_taxonomy_vector]

		for i, target_i in enumerate(labels):
			vector_i = targets_taxonomy_vector[target_i]
			for j, target_j in enumerate(labels):
				vector_j = targets_taxonomy_vector[target_j]
				if i >= j:
					dist = 0
					for k, level in enumerate(vector_i):
						if level != vector_j[k]:
							dist = 6-k
							break
					
					if dist == 1: # it means they correspond to different species. check if they are just not different strains and fuse them if so
						if len(list(set(vector_i[-1].split()) & set(vector_j[-1].split()))) >= 2: # it means they share the genus and the species at least
							dist = 0						
						
					taxonomy_distance_matrix[i][j] = dist
					taxonomy_distance_matrix[j][i] = dist

	elif mode == 'as_input':
		labels = list(targets_taxonomy_vector.keys())
		taxonomy_distance_matrix = [[0 for label in targets_taxonomy_vector] for label in targets_taxonomy_vector]

		if input_targets != None:
			labels = [in_target for in_target in input_targets if in_target in labels]
	
	taxonomy_distance_matrix = np.array(taxonomy_distance_matrix)
	
	return taxonomy_distance_matrix, labels 

#@Timer("get_phylogeny_distance_matrix", text=lambda secs: f"Time to get phylogeny distance matrix: {format_timespan(secs)}")
def get_phylogeny_distance_matrix(in_tree, tree_format = None, input_targets = None, print_tree = True):

	if starting_directory not in in_tree:
		in_tree = '{}/{}'.format(starting_directory, in_tree)

	tree = Phylo.read(in_tree, tree_format)  

	if print_tree:
		print(' ... Input tree: \n')
		Phylo.draw_ascii(tree)

	tree.ladderize()
	tree_depths = tree.depths()	
	if not max(tree_depths.values()): 
		tree_depths = tree.depths(unit_branch_lengths=True) 

	max_depth = max([tree_depths[clade] for clade in tree_depths if clade.name != None])

	# align everything to the right by extending the branches
	all_leafs = []
	for clade in tree.get_terminals():
		curr_depth = tree_depths[clade]
		clade.branch_length += max_depth - curr_depth
		all_leafs.append(clade.name)

	# remove branches that are not in our dataset
	for leaf in all_leafs:
		if leaf not in input_targets:
			tree.prune(leaf)

	labels = [leaf.name for leaf in tree.get_terminals()]

	distance_matrix = [[0 for leaf in labels] for leaf in labels]
	for i, leaf_i in enumerate(labels):
		for j, leaf_j in enumerate(labels):
			distance_matrix[i][j] = tree.distance(leaf_i, leaf_j)

	distance_matrix = np.array(distance_matrix)

	Phylo.write(tree, 'out_tree.nwk', 'newick')

	return distance_matrix, labels

#@Timer("get_operons_distance_matrix", text=lambda secs: f"Time to get operons distance matrix: {format_timespan(secs)}")
def get_operons_distance_matrix(operons, input_targets = None, advanced = False):

	labels = []

	if 'operon_filtered_PaCMAP' in operons[list(operons.keys())[0]]:
		paCMAP_coordinat = []
		for operon_type in operons:
			for i, target in enumerate(operons[operon_type]['target_members']):
				if input_targets is None or target in input_targets:
					paCMAP_coordinat.append(operons[operon_type]['operon_filtered_PaCMAP'][i])
					labels.append(target)

		distance_matrix = [[0 for i in paCMAP_coordinat] for j in paCMAP_coordinat]
		for i, vector_i in enumerate(paCMAP_coordinat):
			for j, vector_j in enumerate(paCMAP_coordinat):
				if i < j:
					dist = np.linalg.norm(np.array(vector_i) - np.array(vector_j))
					distance_matrix[i][j] = dist
					distance_matrix[j][i] = dist

	else:
		genomic_contexts = []
		for operon_type in operons:
			for i, target in enumerate(operons[operon_type]['target_members']):
				if input_targets is None or target in input_targets:
					genomic_contexts.append(operons[operon_type]['operon_protein_families_structure'][i])
					labels.append(target)

		distance_matrix = [[0 for i in genomic_contexts] for j in genomic_contexts]
		for i, vector_i in enumerate(genomic_contexts):
			for j, vector_j in enumerate(genomic_contexts):
				if i < j:
					families_present = set(vector_i+vector_j)
					overlapping_families = set(vector_i).intersection(set(vector_j))
					dist = 1-len(overlapping_families)/len(families_present)
					distance_matrix[i][j] = dist
					distance_matrix[j][i] = dist

	return distance_matrix, labels

#@Timer("get_avgoperons_distance_matrix", text=lambda secs: f"Time to get avgoperons distance matrix: {format_timespan(secs)}")
def get_avgoperons_distance_matrix(operons):
		
	matrix = [[0 for i in operons if '-' not in i] for i in operons if '-' not in i]
	
	selected_operons = [i for i in operons if '-' not in i]
	selected_operons = sorted(selected_operons)
	for i, i_operon_type in enumerate(selected_operons):
		for j, j_operon_type in enumerate(selected_operons):
			if i > j:
				dists = []
				for i_member_pacmap in operons[i_operon_type]['operon_filtered_PaCMAP']:
					for j_member_pacmap in operons[j_operon_type]['operon_filtered_PaCMAP']:
						dist = np.linalg.norm(np.array(i_member_pacmap)-np.array(j_member_pacmap))
						dists.append(dist)
#				 dist = np.linalg.norm(np.array(operons[i_operon_type]['operon_centroid_PaCMAP'])-np.array(operons[j_operon_type]['operon_centroid_PaCMAP']))
				dist = min(dists)
				matrix[i][j] = dist
				matrix[j][i] = dist
	
	return np.array(matrix), selected_operons

#@Timer("normalize_matrix", text=lambda secs: f"Time to normalize matrix: {format_timespan(secs)}")
def normalize_matrix(similarity_matrix, power=30):

	min_dist = pd.DataFrame(similarity_matrix).min().min()
	max_dist = pd.DataFrame(similarity_matrix).max().max()
	
	normalised_matrix = 1-(np.array(similarity_matrix)-min_dist)/max_dist
	normalised_matrix = np.power(normalised_matrix, power)
	
	return normalised_matrix

#@Timer("compute_dendogram", text=lambda secs: f"Time to compute dendogram: {format_timespan(secs)}")
def compute_dendogram(taxonomy, operons, input_targets = None, mode = 'taxonomy', tree = None, tree_format = None, distance_matrix = None, labels = None):
	
	if distance_matrix is None:
		if tree != None:
			distance_matrix, labels = get_phylogeny_distance_matrix(tree, input_targets = input_targets, tree_format = tree_format)
		elif mode == 'operon':
			distance_matrix, labels = get_operons_distance_matrix(operons, input_targets = input_targets)
		elif mode == 'operon clusters':
			distance_matrix, labels = get_avgoperons_distance_matrix(operons)		  
		else:
			distance_matrix, labels = get_taxonomy_distance_matrix(taxonomy, operons, input_targets = input_targets, mode = mode)

	distance_matrix = distance.squareform(distance_matrix)

	Z = hierarchy.linkage(distance_matrix, method = 'single')
	results = hierarchy.dendrogram(Z, no_plot=True, count_sort = 'descending')

	if mode == 'taxonomy' or tree != None or 'operon' in mode:
		labels_indeces = list(map(int, results['ivl']))
		labels = [labels[i] for i in labels_indeces]

	return results, labels

#@Timer("create_dendogram_features", text=lambda secs: f"Time to create dendogram features: {format_timespan(secs)}")
def create_dendogram_features(dendogram, leaf_labels, taxonomy, operons = None, mode = 'taxonomy', colors = None):

	if mode == 'operon clusters':
		data = {'cluster size': [],
				'x': [],
				'y': [],
				'color': [],
				'leaf_label': leaf_labels}
		
		tooltips = [('Cluster type', '@leaf_label'),
					('Cluster size', '@cluster_size')]
		
	else:
		data = {'superkingdom': [],
				'phylum': [],
				'class': [],
				'order': [],
				'genus': [],
				'species': [],
				'x': [],
				'y': [],
				'color': [],
				'leaf_label': leaf_labels}
		
		tooltips = [('Superkingdom', '@superkingdom'),
					('Phylum', '@phylum'),
					('Class', '@class'),
					('Order', '@order'),
					('Genus', '@genus'),
					('Species', '@species')]

	icoord, dcoord = dendogram['icoord'], dendogram['dcoord']

	data['y'] = list(np.linspace(min([num for sublist in icoord for num in sublist]), max([num for sublist in icoord for num in sublist]), len(leaf_labels)))
	data['x'] = [1 for y in data['y']] 

	for label in leaf_labels:
		if colors is not None and label in colors:
			color = colors[label]['Color (RGBA)']
		else:
			color = 'darkgrey'
		
		data['color'].append(color)
		
		if mode == 'operon clusters':
			data['cluster size'].append(str(len(operons[label]['target_members'])))
		else:
			found_taxonomy = False
			for superkingdom in taxonomy:
				for phylum in taxonomy[superkingdom]:
					for taxclass in taxonomy[superkingdom][phylum]:
						for order in taxonomy[superkingdom][phylum][taxclass]:
							for genus in taxonomy[superkingdom][phylum][taxclass][order]:
								for species in taxonomy[superkingdom][phylum][taxclass][order][genus]:
									if species == label:
										data['superkingdom'].append(superkingdom)
										data['phylum'].append(phylum)
										data['class'].append(taxclass)
										data['order'].append(order)
										data['genus'].append(genus)
										found_taxonomy = True
									else:
										for target in taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
											if target == label:
												data['superkingdom'].append(superkingdom)
												data['phylum'].append(phylum)
												data['class'].append(taxclass)
												data['order'].append(order)
												data['genus'].append(genus)
												data['species'].append(species)
												found_taxonomy = True
			if not found_taxonomy:
				data['superkingdom'].append('na')
				data['phylum'].append('na')
				data['class'].append('na')
				data['order'].append('na')
				data['genus'].append('na')
				data['species'].append('na')

	return tooltips, data

#@Timer("make_dendogram_figure", text=lambda secs: f"Time to make dendogram figure: {format_timespan(secs)}")
def make_dendogram_figure(taxonomy, operons, input_targets = None, tree = None, mode = 'taxonomy', show_leafs = True, height_factor = 20, tree_format = None, distance_matrix = None, labels = None, colors = None):

	dendogram, leaf_labels = compute_dendogram(taxonomy, operons, input_targets = input_targets, mode = mode, tree = tree, tree_format = tree_format, distance_matrix = distance_matrix, labels = labels)
	den_tooltips, den_data = create_dendogram_features(dendogram, leaf_labels, taxonomy, operons = operons, mode = mode, colors = colors)

	y_range = [0, max(den_data['y'])+min(den_data['y'])+(den_data['y'][1]-den_data['y'][0])]

	if len(leaf_labels) < 5:
		height = int(len(leaf_labels)*height_factor)*3
	elif len(leaf_labels) < 10:
		height = int(len(leaf_labels)*height_factor)*2
	else:
		height = int(len(leaf_labels)*height_factor)

	if tree != None:
		title = 'Input phylogeny/hierarchy'
	elif mode == 'taxonomy':
		title = 'Taxonomy hierarchy'
	elif mode == 'operons':
		title = 'Genomic contexts similarity hierarchy'
	else:
		title = 'Genomic contexts clusters hierarchy'

	if show_leafs:
		den = figure(title = title, height=height, width=500, x_range=[-max(max(dendogram['dcoord']))-(max(max(dendogram['dcoord']))-min(min(dendogram['dcoord'])))*0.2, 40], y_range = y_range, toolbar_location="left")
	else:
		den = figure(title = title, height=height, width=250, x_range=[-max(max(dendogram['dcoord']))-(max(max(dendogram['dcoord']))-min(min(dendogram['dcoord'])))*0.2, 0.2], y_range = y_range, toolbar_location="left")

	for i, d in zip(dendogram['icoord'], dendogram['dcoord']):
		d = list(map(lambda x: -x, d))
		den.line(x=d, y=i, line_color='black')

	if show_leafs:
		leafs = den.text(x='x', y='y', text = 'leaf_label', text_baseline='middle', text_font_size='8pt', source = den_data)
	else:
		leafs = den.circle(x=0.1, y='y', line_color = 'black', fill_color = 'color', source = den_data, size = (den_data['y'][1] - den_data['y'][0])*0.7)
	den.add_tools(HoverTool(tooltips=den_tooltips, renderers=[leafs]))

	den.title.align = "right"
	den.axis.major_tick_line_color = None
	den.axis.minor_tick_line_color = None
	den.axis.major_label_text_color = None
	den.axis.major_label_text_font_size = '0pt'
	den.axis.axis_line_color = None
	den.grid.grid_line_color = None
	den.outline_line_color = None

	return den, den_data

# 9.3. Routines to draw all genomic contexts

#@Timer("create_genomic_context_features", text=lambda secs: f"Time to create genomic context features: {format_timespan(secs)}")
def create_genomic_context_features(operons, all_syntenies, family_colors, syn_den_data, reference_family, legend_mode = 'ncbi_code'):
	
	data = {'operon':[],
			'target_id':[],
			'ncbi_code':[],
			'ncbi_link':[],
			'assembly': [],
			'name':[],
			'family': [],
			'superkingdom': [],
			'phylum': [],
			'class': [],
			'order': [],
			'genus': [],
			'species': [],
			'relative_start': [],
			'relative_end': [],
			'facecolor': [],
			'edgecolor': [],
			'linestyle': [],
			'xs':[],
			'ys':[],
			'half_heights': [],
			'text_x': [],
			'text_y': [],
			'tm_text_x': [],
			'tm_text_y': [],
			'tm_text': [],
			'tm_pred_text': []}
	
	yyticklabels = []
	yys = []
	y_step = syn_den_data['y'][1] - syn_den_data['y'][0]
	y_half_height = y_step/4
	
	for i, current_target in enumerate(syn_den_data['leaf_label']):
		operon = [operon for operon in operons if current_target in operons[operon]['target_members']][0]
		current_assembly = all_syntenies[current_target]['assembly_id'][1]
		current_species = all_syntenies[current_target]['species']
		current_genomic_context_block = all_syntenies[current_target]['flanking_genes']
		current_species = all_syntenies[current_target]['species']
		current_reference_family = all_syntenies[current_target]['target_family']
		
		curr_y = syn_den_data['y'][i]

		for j, flanking_gene in enumerate(current_genomic_context_block['ncbi_codes']):
			family = current_genomic_context_block['families'][j]
			name = current_genomic_context_block['names'][j]
			gene_dx = current_genomic_context_block['relative_ends'][j] - current_genomic_context_block['relative_starts'][j]+1
			gene_direction = current_genomic_context_block['directions'][j]

			if gene_direction == '-':
				gene_x_tail = current_genomic_context_block['relative_ends'][j]
				gene_dx = gene_dx*(-1)
				gene_x_head = gene_x_tail + gene_dx
				gene_x_head_start = gene_x_head+100
				text_x = gene_x_tail - (gene_x_tail-gene_x_head_start)/2
			else:
				gene_x_tail = current_genomic_context_block['relative_starts'][j]
				gene_x_head = gene_x_tail + gene_dx
				gene_x_head_start = gene_x_head-100
				text_x = gene_x_tail + (gene_x_head_start-gene_x_tail)/2

			if family == 0:
				facecolor = family_colors[family]['Color (RGBA)']
				edgecolor = family_colors[family]['Line color']
				linestyle = family_colors[family]['Line style']
			else:
				facecolor = family_colors[family]['Color (RGBA)']
				edgecolor = family_colors[family]['Line color']
				linestyle = family_colors[family]['Line style'] 
				
			data['operon'].append(operon)
			data['target_id'].append(current_target)
			data['ncbi_code'].append(flanking_gene)
			data['assembly'].append(current_assembly)
			data['name'].append(name)
			data['relative_start'].append(format(gene_x_tail, ',d'))
			data['relative_end'].append(format(gene_x_head, ',d'))
			data['facecolor'].append(facecolor)
			data['edgecolor'].append(edgecolor)
			data['linestyle'].append(linestyle)
			data['xs'].append([gene_x_tail, gene_x_tail, gene_x_head_start, gene_x_head, gene_x_head_start])
			data['ys'].append([curr_y-y_half_height, curr_y+y_half_height, curr_y+y_half_height, curr_y, curr_y-y_half_height])
			data['half_heights'].append(y_half_height)
			data['text_x'].append(text_x)
			data['text_y'].append(curr_y+y_half_height)
			data['ncbi_link'].append('https://www.ncbi.nlm.nih.gov/protein/{}'.format(flanking_gene))
			data['tm_text_x'].append(text_x)
			data['tm_text_y'].append(curr_y)
			
			if 'TM_annotations' in current_genomic_context_block:
				tm_annotation = current_genomic_context_block['TM_annotations'][j]
				if tm_annotation == 'TM':
					data['tm_pred_text'].append('Yes')
				elif tm_annotation == 'SP':
					data['tm_pred_text'].append('Contains signal peptide')
				else:
					data['tm_pred_text'].append('No')
			else:
				tm_annotation = ''
				data['tm_pred_text'].append('n.a.')
				
			data['tm_text'].append(tm_annotation)
								  
			if family != 0 and family != reference_family and family < 10000:
				data['family'].append(family)
			else:
				data['family'].append(str(''))
				
			data['species'].append(current_species)
			for level in ['superkingdom', 'phylum', 'class', 'order', 'genus']:
				data[level].append(syn_den_data[level][i])
		
		if legend_mode in ['superkingdom', 'phylum', 'class', 'order', 'genus', 'species']:
			label = data[legend_mode][-1]
			if label == []:
				label = 'n.a'
			yyticklabels.append(label)
		else:
			yyticklabels.append('{} | {}'.format(current_target, operon))
			
		yys.append(curr_y)
		
	yyticklabels = {int(yys[i]): yyticklabels[i] for i in range(len(yyticklabels))}

	tooltips = [('GC type', "@operon"),
				('InputID', "@target_id"),
				('EntrezID', '@ncbi_code'),
				('Genome assembly', '@assembly'),
				('Gene relative start', '@relative_start'),
				('Gene relative end', '@relative_end'),
				("Protein name", "@name"),
				("Protein family code", "@family"),
				("Predicted membrane protein", "@tm_pred_text"),
				('Superkingdom', '@superkingdom'),
				('Phylum', '@phylum'),
				('Class', '@class'),
				('Order', '@order'),
				('Genus', '@genus'),
				("Species", "@species")] 
		
	return tooltips, data, yyticklabels


#@Timer("create_genomic_context_figure", text=lambda secs: f"Time to create genomic context figure: {format_timespan(secs)}")
def create_genomic_context_figure(operons, all_syntenies, family_colors, syn_den_data, syn_dendogram, most_common_gc_figure, reference_family, legend_mode = 'ncbi_code', height_factor = 25*1.2):

	p_tooltips, p_data, p_yyticklabels = create_genomic_context_features(operons, all_syntenies, family_colors, syn_den_data, reference_family = reference_family, legend_mode = legend_mode)
	
	p = figure(plot_width=most_common_gc_figure.plot_width, plot_height=syn_dendogram.height, x_range = most_common_gc_figure.x_range, y_range = syn_dendogram.y_range, toolbar_location="left", title = 'Representative genomic contexts (hover to get more information)')  # the genomic_context figure

	for i, xs in enumerate(p_data['xs']):
		p.patch(xs, p_data['ys'][i], fill_color = p_data['facecolor'][i], line_color = p_data['edgecolor'][i], line_width = 1)

	p.patches('xs', 'ys', fill_color = None, line_color = None, line_width = 0, source = p_data, 
			  hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
			  selection_fill_color='facecolor', selection_line_color='edgecolor',
			  nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', nonselection_fill_alpha=0.2)

	p.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", text_font_size = {'value': '6pt'}, source = p_data)
	p.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, source = p_data)
	
	# define yticks on the left
	p.yaxis.ticker = [int(n) for n in list(p_yyticklabels.keys())]
	p.yaxis.major_tick_line_color = None
	p.yaxis.major_label_overrides = {int(i): p_yyticklabels[i] for i in [int(n) for n in p_yyticklabels.keys()]}
#	 p.yaxis.major_label_text_font_size = {'value': '8pt'}
#	 p.yaxis.major_label_text_font_style = 'italic'
	p.yaxis.axis_line_width = 0
	
	# define xticks
	p.xaxis.axis_label = "Position relative to target (bp)"

	# define general features
	p.grid.visible = False
	p.outline_line_width = 0

	p.add_tools(HoverTool(tooltips=p_tooltips))
	p.add_tools(TapTool(callback = OpenURL(url='@ncbi_link')))
		
	return p

# 9.4. Routines to draw legend figure

#@Timer("create_legend_features", text=lambda secs: f"Time to create legend features: {format_timespan(secs)}")
def create_legend_features(operons, families_summary, reference_family, family_colors, dx = 50):
	
	data = {'xs': [],
			'ys': [],
			'edgecolor': [],
			'facecolor': [],
			'text': [],
			'label_x': [],
			'label_y': [],
			'text_x': [],
			'text_y': [],
			'tm_text_x': [],
			'tm_text_y': [],
			'tm_text': [],
			'tm_type': [],
			'tm_mode': [],
			'family': [],
			'found_models': [],
			'model_links': [],
			'keywords': [],
			'go_terms': [],
			'function': []}
	
	
	curr_y = len(family_colors.keys())
	
	for family in sorted(list(families_summary.keys())):
		if family == 0:
			data['text'].append('Non-conserved gene')
		elif family == reference_family:
			data['text'].append('Target protein: {}'.format(families_summary[family]['name']))
		elif family == 10000:
			data['text'].append('Pseudogene')
		elif family in families_summary:
			data['text'].append(families_summary[family]['name'])
		
		if family != 0 and family != reference_family and family < 10000:
			data['family'].append(family)
		else:
			data['family'].append(str(''))
		
		if 'model_state' in families_summary[family]:
			model_state = families_summary[family]['model_state']

			if model_state == 'Model exists':
				model_state = families_summary[family]['structure']
			elif model_state == 'Model does not exist':
				model_state = 'click to model/view with Swiss-Model'
			else:
				if family > 0 and family < 10000:
					model_state = 'Not possible to find'
				else:
					model_state = 'n.a.'
		else:
			model_state = ''
			
		data['found_models'].append(model_state)
		
		if family > 0 and family < 10000 and 'function' in families_summary[family]:
			if 'TM_topology' in families_summary[family]['function']:
				tm_type = families_summary[family]['function']["TM_topology"]
				keywords = ', '.join(sorted(families_summary[family]['function']['Keywords']))
				go_terms = '; '.join(sorted(families_summary[family]['function']['GO_terms']))
				function = families_summary[family]['function']['Function_description']

				if len(tm_type) > 0:
					tm_text = 'TM'
					tm_mode = 'Yes -> type:'
				else:
					tm_text = ''
					tm_mode = 'No'
			else:
				tm_type = ''
				tm_text = ''
				tm_mode = ''
				keywords = ''
				go_terms = ''   
				function = ''
		else:
			tm_type = 'n.a.'
			tm_text = ''
			tm_mode = 'n.a.'
			keywords = 'n.a.'
			go_terms = 'n.a.'   
			function = 'n.a.'
		
		if 'structure' in families_summary[family]:
			structure = families_summary[family]['structure']
			if structure == '':
				uniprot_code = families_summary[family]['uniprot_code']
				structure = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
		else:
			structure = 'n.a.'
			
		data['model_links'].append(structure)
		
		data['facecolor'].append(family_colors[family]['Color (RGBA)'])
		data['edgecolor'].append(family_colors[family]['Line color'])
		data['xs'].append([5, 5, dx-5, dx, dx-5])
		data['ys'].append([curr_y-0.25, curr_y+0.25, curr_y+0.25, curr_y, curr_y-0.25])
		data['text_x'].append(((dx-10)/2)+5)
		data['text_y'].append(curr_y+0.25)
		
		data['label_x'].append(dx + 10) 
		data['label_y'].append(curr_y) 
		
		data['tm_text_x'].append(((dx-10)/2)+5)
		data['tm_text_y'].append(curr_y)
		data['tm_text'].append(tm_text)
		data['tm_type'].append(tm_type)
		data['tm_mode'].append(tm_mode)
		
		data['go_terms'].append(go_terms)
		data['keywords'].append(keywords)
		data['function'].append(function)
		
		curr_y -= 1
		
	tooltips = [('Protein family', '@text'),
				('Protein family code', '@family'),
				('Structure', '@found_models'),
				('Predicted membrane protein', '@tm_mode @tm_type'),
				('Keywords', '@keywords'),
				('GO terms', '@go_terms'),
				('Function', '@function')]
	
	return tooltips, data

#@Timer("create_legend_figure", text=lambda secs: f"Time to create legend figure: {format_timespan(secs)}")
def create_legend_figure(operons, families_summary, reference_family, family_colors, grid, rescale_height = False):
	
	l_tooltips, l_data = create_legend_features(operons, families_summary, reference_family, family_colors)
	
	height = int(len(family_colors)*25*1.2)
	
	if len(l_data['ys']) > len(operons) or rescale_height:
		height = grid.height
		
	l = figure(plot_width=500, plot_height=height, x_range = [0,500], title = 'Protein families (click to view in Swiss-Model Repository)', toolbar_location="right")  # the legend figure
	
	for i, xs in enumerate(l_data['xs']):
		l.patch(xs, l_data['ys'][i], fill_color = l_data['facecolor'][i], line_color = l_data['edgecolor'][i], line_width = 1)	
	
	l.patches('xs', 'ys', fill_color = None, line_color = None, line_width = 0, source = l_data, 
			  hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
			  selection_fill_color='facecolor', selection_line_color='edgecolor',
			  nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', nonselection_fill_alpha=0.2)

	l.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", text_font_size = {'value': '6pt'}, source = l_data)
	l.text('label_x', 'label_y', text = 'text', text_baseline="middle", text_align="left", text_font_size = {'value': '8pt'}, source = l_data)
	l.text('tm_text_x', 'tm_text_y', text = 'tm_text',  text_color = "white", text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, source = l_data)
	
	l.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
	l.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
	l.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	l.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	l.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
	l.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
	l.yaxis.axis_line_width = 0
	l.xaxis.axis_line_width = 0
	# define general features
	l.grid.visible = False
	l.outline_line_width = 0
	
	l.add_tools(HoverTool(tooltips=l_tooltips))
	l.add_tools(TapTool(callback = OpenURL(url='@model_links')))
			
	return l

# 9.5. Routines to make the gene correlation figure

#@Timer("get_coocurrence_matrix", text=lambda secs: f"Time to get coocurrence matrix: {format_timespan(secs)}")
def get_coocurrence_matrix(operons, families_summary, min_coocc = 0.40):
	
	family_labels = sorted([family for family in families_summary.keys() if family > 0 and family < 10000])
	matrix = [[0 for family in family_labels] for family in family_labels]

	context_count = 0
	for operon in operons:
		for genomic_context in operons[operon]['operon_protein_families_structure']:
			for i in range(len(set(genomic_context))):
				for j in range(len(set(genomic_context))):
					if i > j:
						out_family = list(set(genomic_context))[i]
						in_family = list(set(genomic_context))[j]

						if all(family > 0 and family < 10000 for family in [out_family, in_family]):
							matrix[family_labels.index(out_family)][family_labels.index(in_family)] += 1
							matrix[family_labels.index(in_family)][family_labels.index(out_family)] += 1


	matrix = np.array(matrix)
	matrix = matrix/np.amax(matrix)
	matrix = np.where(matrix < min_coocc, 0, matrix)
	matrix = np.where(matrix!=0, (np.exp(matrix*2)-1)*1, matrix)

	# remove all columns and lines with all zero
	indices_with_all_zeros = np.all(matrix == 0, axis=1)

	matrix = matrix[~indices_with_all_zeros]
	matrix = matrix.T[~indices_with_all_zeros]

	selected_families = np.array(family_labels)[~indices_with_all_zeros]
	selected_families_summary = {family: families_summary[family] for family in selected_families}
	
	return matrix, selected_families_summary

#@Timer("get_adjcency_matrix", text=lambda secs: f"Time to get adjcency matrix: {format_timespan(secs)}")
def get_adjcency_matrix(operons, families_summary):
	
	family_labels = sorted([family for family in families_summary.keys() if family > 0 and family < 10000])
	matrix = [[0 for family in family_labels] for family in family_labels]
	
	context_count = 0
	for operon in operons:
		for genomic_context in operons[operon]['operon_protein_families_structure']:
			for i in range(len(genomic_context)-1):
				out_family = genomic_context[i]
				in_family = genomic_context[i+1]

				if all(family > 0 and family < 10000 for family in [out_family, in_family]):
					matrix[family_labels.index(out_family)][family_labels.index(in_family)] += 1
					matrix[family_labels.index(in_family)][family_labels.index(out_family)] += 1

	matrix = np.array(matrix)
	
	return matrix, family_labels

#@Timer("get_graph_from_matrix", text=lambda secs: f"Time to get graph from matrix: {format_timespan(secs)}")
def get_graph_from_matrix(matrix, selected_families_summary, family_colors):
	
	G=nx.from_numpy_array(matrix)

	# take care of the edges
	edge_params = {'color': {}, 'weight': {}, 'line_width': {}}
	edge_cmap = matplotlib.cm.get_cmap('Greys')
	edge_norm = matplotlib.colors.Normalize(vmin = 0, vmax = 4)

	for start_node, end_node, params in G.edges(data=True):
		edge_color = edge_cmap(edge_norm(round(params['weight'])))
		edge_params['line_width'][(start_node, end_node)] = params['weight']
		edge_params['color'][(start_node, end_node)] = RGB(int(255*list(edge_color)[0]), int(255*list(edge_color)[1]), int(255*list(edge_color)[2]))
		edge_params['weight'][(start_node, end_node)] = 1
		
	# take care of the nodes
	node_params = {'color': {}, 'label': {}}
	for node, _ in G.nodes(data=True):
		node_label = sorted(list(selected_families_summary.keys()))[node]
		node_params['label'][node] = node_label

		node_color = family_colors[node_label]['Color (RGBA)']
		node_params['color'][node] = node_color

	nx.set_node_attributes(G, node_params['color'], "node_color")
	nx.set_node_attributes(G, node_params['label'], "node_label")
	nx.set_edge_attributes(G, edge_params['color'], "edge_color")
	nx.set_edge_attributes(G, edge_params['line_width'], "line_width")
	nx.set_edge_attributes(G, edge_params['weight'], "weight")
	
	return G

#@Timer("remove_non_adjacent_edges", text=lambda secs: f"Time to remove non adjacent edges: {format_timespan(secs)}")
def remove_non_adjacent_edges(operons, families_summary, G, families_present):
	
	adjacency_matrix, family_labels = get_adjcency_matrix(operons, families_summary)
	
	edges_to_remove = []
	for start_node, end_node, params in G.edges(data=True):
		start_node_index = family_labels.index(families_present[start_node])
		end_node_index = family_labels.index(families_present[end_node])
		
		if adjacency_matrix[start_node_index][end_node_index] == 0:
			edges_to_remove.append((start_node, end_node))
	
	for edge in edges_to_remove:
		G.remove_edge(*edge)
	
	return G

#@Timer("create_node_features", text=lambda secs: f"Time to create node features: {format_timespan(secs)}")
def create_node_features(families_summary, reference_family, node_graph_coords, G):
	
	data = {'text_x': [],
			'text_y': [],
			'family': [],
			'tm_text_x': [],
			'tm_text_y': [],
			'tm_text': [],
			'tm_type': [],
			'tm_pred_text': [],
			'protein_name': [],
			'found_models': [],
			'model_links': []}
	
	for node in node_graph_coords:
		
		family = G.nodes[node]['node_label']
		coords = node_graph_coords[node]
		protein_name = families_summary[family]['name']
		
		if family > 0 and family < 10000 and 'function' in families_summary[family]:
			if 'TM_topology' in families_summary[family]['function']:
				tm_type = families_summary[family]['function']["TM_topology"]
				
				if len(tm_type) > 0:
					tm_text = 'TM'
					tm_mode = 'Yes -> type:'
				else:
					tm_text = ''
					tm_mode = 'No'
			else:
				tm_type = ''
				tm_text = ''
				tm_mode = ''
		else:
			tm_type = 'n.a.'
			tm_text = ''
			tm_mode = 'n.a.'
		
		if 'model_state' in families_summary[family]:
			model_state = families_summary[family]['model_state']

			if model_state == 'Model exists':
				model_state = 'Yes (click to view in Swiss-Model repository)'
			elif model_state == 'Model does not exist':
				model_state = 'No (click to model with Swiss-Model)'
			else:
				if family > 0 and family < 10000:
					model_state = 'Not possible to find'
				else:
					model_state = ''
			
			structure = families_summary[family]['structure']
			if structure == '':
				uniprot_code = families_summary[family]['uniprot_code']
				structure = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
			
		else:
			model_state = 'n.a.'
			structure = 'n.a.'
		
		if family != 0 and family != reference_family and family < 10000:
			data['family'].append(family)
		else:
			data['family'].append(str(''))
			
		data['found_models'].append(model_state)
		data['model_links'].append(structure)
		
		y_range = (min([node_graph_coords[node][1] for node in node_graph_coords])-0.5, max([node_graph_coords[node][1] for node in node_graph_coords])+0.5)
		y_step = (y_range[1]-y_range[0])*0.09

		data['text_x'].append(coords[0])
		data['text_y'].append(coords[1]+y_step)
		
		data['tm_text_x'].append(coords[0])
		data['tm_text_y'].append(coords[1])
		data['tm_text'].append(tm_text)
		data['tm_pred_text'].append(tm_mode)
		data['tm_type'].append(tm_type)
		
		data['protein_name'].append(protein_name)
		
	
	tooltips = [('Protein name', "@protein_name"),
				('Protein family code', '@family'),
				("Predicted membrane protein", "@tm_pred_text @tm_type"),
				('Structural model found', '@found_models')] 
	
	return data, tooltips

#@Timer("create_graph_figure", text=lambda secs: f"Time to create graph figure: {format_timespan(secs)}")
def create_graph_figure(operons, reference_family, families_summary, family_colors, gc, min_coocc = 0.40, mode = 'coocurrence', graph_coord = {}, previous_net = ''):
	
	matrix, selected_families_summary = get_coocurrence_matrix(operons, families_summary, min_coocc = min_coocc)
	graph = get_graph_from_matrix(matrix, selected_families_summary, family_colors)
	
	if mode == 'coocurrence':
		title = 'Gene co-occurrence network'
	elif mode == 'adjacency':
		graph = remove_non_adjacent_edges(operons, families_summary, graph, families_present = sorted(list(selected_families_summary.keys())))
		title = 'Gene adjcency network'
	
	if len(graph_coord) == 0:
		graph_coord = nx.spring_layout(graph)
	
	node_data, node_tooltips = create_node_features(selected_families_summary, reference_family, graph_coord, graph)
		
	if previous_net != '':
		x_range = previous_net.x_range
		y_range = previous_net.y_range
		
	else:
		x_range = (min([graph_coord[node][0] for node in graph_coord])-0.5, max([graph_coord[node][0] for node in graph_coord])+0.5)
		y_range = (min([graph_coord[node][1] for node in graph_coord])-0.5, max([graph_coord[node][1] for node in graph_coord])+0.5)

	g = figure(width = gc.plot_width, height = gc.plot_height, x_range=x_range, y_range=y_range, title = title)

	graph_renderer = from_networkx(graph, graph_coord, scale=1, center=(0, 0))
	graph_renderer.edge_renderer.glyph = MultiLine(line_width="line_width", line_color = "edge_color")
	graph_renderer.node_renderer.glyph = Circle(size=22, fill_color = "node_color")
	
	g.renderers.append(graph_renderer)
	
	g.text('text_x', 'text_y', text = 'family', text_baseline="bottom", text_align="center", text_font_size = {'value': '6pt'}, source = node_data)
	g.text('tm_text_x', 'tm_text_y', text = 'tm_text', text_color = "white", text_baseline="middle", text_align="center", text_font_size = {'value': '6pt'}, source = node_data)
	g.circle('tm_text_x', 'tm_text_y', color = None, size = 22, source = node_data)
	
	g.add_tools(HoverTool(tooltips=node_tooltips))
	g.add_tools(TapTool(callback = OpenURL(url='@model_links')))

	g.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
	g.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
	g.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	g.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	g.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
	g.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
	g.yaxis.axis_line_width = 0
	g.xaxis.axis_line_width = 0
	# define general features
	g.grid.visible = False
	g.outline_line_width = 0
	
	return g, graph_coord

# 9.6. Routines for advenced interactive plotting

#@Timer("define_operons_colors", text=lambda secs: f"Time to define operons colors: {format_timespan(secs)}")
def define_operons_colors(clusters, mode = 'matplotlib', cmap = 'rainbow', print_summary = False):

	colors = {}

	cmap = matplotlib.cm.get_cmap(cmap)
	norm = matplotlib.colors.Normalize(vmin=0, vmax=len(clusters))

	colours = [cmap(norm(i)) for i in range(len(clusters))]
	#random.shuffle(colours)

	if print_summary:
		print(" ... Setting the colors of the operon clusters identified")
		print(" ... ... Colormap: {}".format(cmap))
		print('\nLabel\tColor (RGBA)')

	for i, label in enumerate(sorted(clusters)):
		if label not in colors:
			colors[label] = {}

		if '-' in label:   # singleton clusters
			colors[label]['Color (RGBA)'] = 'grey'
			colors[label]['Color (tuplet)'] = 'grey'
			colors[label]['Line color'] = 'black'
			colors[label]['Line style'] = '-'
			colors[label]['Size'] = '2'

		else:  
			if mode == 'matplotlib':
				colors[label]['Color (RGBA)'] = [int(255*j) for j in colours[i]]
				colors[label]['Color (tuplet)'] = colours[i]
				colors[label]['Line color'] = 'black'
				colors[label]['Line style'] = '-'
			elif mode == 'bokeh':
				colors[label]['Color (RGBA)'] = RGB(int(255*list(colours[i])[0]), int(255*list(colours[i])[1]), int(255*list(colours[i])[2]))
				colors[label]['Color (tuplet)'] = colours[i]
				colors[label]['Line color'] = 'black'
				colors[label]['Line style'] = '-'
			colors[label]['Size'] = '5'

		if print_summary:
			print('{}\t{}'.format(label, colors[label]['Color (RGBA)']))

	return colors 

#@Timer("parse_coordinates_from_clans", text=lambda secs: f"Time to parse coordinates from clans: {format_timespan(secs)}")
def parse_coordinates_from_clans(clans_file):
	
	clans_coords = {}
	seq_map = {}
	
	found_seq_block = False
	found_coords = False
	seq_count = 0
	with open(clans_file, 'r') as inclans:
		for line in inclans:
			if '<seq>' in line:
				found_seq_block = True
			elif '</seq>' in line:
				found_seq_block = False
			elif found_seq_block and line.startswith('>'):
				line = line[1:]
				ncbi_code = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()
				seq_map[seq_count] = ncbi_code
				
				if '[' in line and ']' in line:
					species = line.split('[')[1].split(']')[0]
				else:
					species = 'n.a.'
				
				clans_coords[ncbi_code]={'species': species, 'xy': None}
				
				seq_count += 1
			
			elif '<pos>' in line:
				found_coords = True
			elif '</pos>' in line:
				found_coords = False
			elif not found_seq_block and found_coords:
				coords = [float(i) for i in line.strip().split(' ')]
				ncbi_code = seq_map[int(coords[0])]
				clans_coords[ncbi_code]['xy'] = coords[1:3]
			
	return clans_coords

#@Timer("generate_coordinates_for_clans", text=lambda secs: f"Time to generate coordinates for clans: {format_timespan(secs)}")
def generate_coordinates_for_clans(all_syntenies, out_label = None, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, min_coverage = None, method = None, mmseqs = None, blast = None, default_base = None, tmp_folder = None):

	in_syntenies = {}
	for target in all_syntenies:
		in_syntenies[target] = {'flanking_genes': {}}

		assembly_targetid = all_syntenies[target]['assembly_id'][0]
		context_idx = all_syntenies[target]['flanking_genes']['ncbi_codes'].index(assembly_targetid)

		for key in all_syntenies[target]['flanking_genes']:
			if type(all_syntenies[target]['flanking_genes'][key]) == list:
				if key == 'ncbi_codes':
					in_syntenies[target]['flanking_genes'][key] = [target]
				else:
					in_syntenies[target]['flanking_genes'][key] = [all_syntenies[target]['flanking_genes'][key][context_idx]]

	distance_matrix, ordered_ncbi_codes = compute_all_agains_all_distance_matrix(in_syntenies, out_label = '{}_targets'.format(out_label), num_threads = num_threads, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, min_coverage = min_coverage, method = method, mmseqs = mmseqs, blast = blast, default_base = default_base, tmp_folder = tmp_folder)	

	paCMAP_embedding = pacmap.PaCMAP(n_components = 2)
	paCMAP_coordinat = paCMAP_embedding.fit_transform(distance_matrix)

	clans_coords = {ordered_ncbi_codes[i]: {'xy': paCMAP_coordinat[i]} for i in range(len(paCMAP_coordinat))}

	return clans_coords

#@Timer("create_data_structure", text=lambda secs: f"Time to create data structure: {format_timespan(secs)}")
def create_data_structure(operons, clans_file, clusters_colors, all_syntenies, out_label = None, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, min_coverage = None, method = None, mmseqs = None, blast = None, default_base = None, tmp_folder = None):

	if clans_file is not None:
		clans_coords = parse_coordinates_from_clans(clans_file)
		seq_in_clans_without_operon = list(set(clans_coords.keys()).difference(list(all_syntenies.keys())))
	else:
		clans_coords = generate_coordinates_for_clans(all_syntenies, out_label = out_label, num_threads = num_threads, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, min_coverage = min_coverage, method = method, mmseqs = mmseqs, blast = blast, default_base = default_base, tmp_folder = tmp_folder)
		seq_in_clans_without_operon = []
	
	data = {'x': [],		  # pacmap x coord for individual operon
			'y': [],		  # pacmap y coord for individual operon
			'avg_x': [],	  # average x coord of the corresponding operon type
			'avg_y': [],	  # average y coord of the corresponding operon type
			'clans_x': [],	# clans x coord for the input target
			'clans_y': [],	# clans y coord for the input target
			'edgecolor': [],
			'facecolor': [],
			'size': [],
			'node_size': [],
			'text': [],
			'type': [],
			'species': [],
			'target': []}
	
	# first go through the sequences in the clans map that were not assigned any gene cluster type 
	# (eg, the clans map is for a much larger set of sequences, our set of sequences is a subset of the clans map
	# or sequences were filteresd out because they were in a partial genome context)
	for member in seq_in_clans_without_operon:
		clans_x, clans_y = clans_coords[member]['xy']
		text = 'Not analysed'
		operon_type = 'n.a.'
		species	 = clans_coords[member]['species']
		facecolor   = 'grey'
		edgecolor   = 'grey'
		size		= 1
	
		data['clans_x'].append(clans_x)
		data['clans_y'].append(-clans_y)
		data['x'].append(np.nan)
		data['y'].append(np.nan)
		data['avg_x'].append(np.nan)
		data['avg_y'].append(np.nan)
		data['text'].append(text)
		data['edgecolor'].append(edgecolor)
		data['facecolor'].append(facecolor)
		data['size'].append(size)
		data['node_size'].append(np.nan)
		data['type'].append(operon_type)
		data['species'].append(species)
		data['target'].append(member)
	
	# now go through the sequences in the gene clusters and collect their xy clans coordinates if they
	# are in the clans map
	for operon_type in sorted(list(operons.keys())):
		facecolor = clusters_colors[operon_type]['Color (RGBA)']
		edgecolor = clusters_colors[operon_type]['Line color']
		size = clusters_colors[operon_type]['Size']
		
		for i, member in enumerate(operons[operon_type]['target_members']):
			if member in clans_coords:
				clans_x, clans_y = clans_coords[member]['xy']
				clans_y = -clans_y
			else:
				clans_x, clans_y = np.nan, np.nan
			
			x, y = operons[operon_type]['operon_filtered_PaCMAP'][i]
			text = operons[operon_type]['operon_protein_families_structure'][i]

			data['clans_x'].append(clans_x)
			data['clans_y'].append(clans_y)
			data['x'].append(x)
			data['y'].append(y)
			data['avg_x'].append(np.nan)
			data['avg_y'].append(np.nan)
			data['text'].append(text)
			data['edgecolor'].append(edgecolor)
			data['facecolor'].append(facecolor)
			data['size'].append(size)
			data['node_size'].append(np.nan)
			data['type'].append(operon_type.split()[2])
			data['species'].append(all_syntenies[member]['species'])
			data['target'].append(member)
		
		# add the data for the average position of the GC in the pacmap space if it is not a -00001 type
		if '-' not in operon_type:
			avg_x, avg_y = operons[operon_type]['operon_centroid_PaCMAP']
			
			data['clans_x'].append(np.nan)
			data['clans_y'].append(np.nan)
			data['x'].append(avg_x)
			data['y'].append(avg_y)
			data['avg_x'].append(avg_x)
			data['avg_y'].append(avg_y)
			data['text'].append(np.nan)
			data['edgecolor'].append(edgecolor)
			data['facecolor'].append(facecolor)
			data['size'].append('1')
			data['node_size'].append('10')
			data['type'].append(operon_type.split()[2])
			data['species'].append(np.nan)
			data['target'].append(np.nan)
				
	tooltips = [('Operon/GC type', '@type'),
				('EntrezID', '@target'),
				('Species', '@species'),]
	
	return tooltips, ColumnDataSource(data)

#@Timer("create_operons_clusters_scatter", text=lambda secs: f"Time to create operons clusters scatter: {format_timespan(secs)}")
def create_operons_clusters_scatter(scatter_data):

	p_tooltips, p_data = scatter_data

	p = figure(title = 'Genomic context types/clusters', plot_width=500, plot_height=500)
	p.add_layout(Legend(orientation="horizontal"), 'above')

#	 p.circle('x', 'y', size='size', line_color='edgecolor', fill_color='facecolor', legend_field='type', alpha=1, source = p_data)
	p.circle('x', 'y', size='size', line_color='edgecolor', fill_color='facecolor', alpha=1, source = p_data)

	p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
	p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
	p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
	p.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
	p.yaxis.axis_line_width = 0
	p.xaxis.axis_line_width = 0
	# define general features
	p.grid.visible = False
#	 p.outline_line_width = 0

	p.xaxis.axis_label = ""
	p.yaxis.axis_label = ""

	p.add_tools(HoverTool(tooltips=p_tooltips))
	p.add_tools(LassoSelectTool())

	p.background_fill_color = "lightgrey"
	p.background_fill_alpha = 0.2

	p.legend.click_policy="hide"

	return p

#@Timer("create_edges_data", text=lambda secs: f"Time to create edges data: {format_timespan(secs)}")
def create_edges_data(scatter_data, operons, alpha_matrix, labels):

	data = {'x': [],
			'y': [],
			'color': [],
			'alpha':[]}

	for i, i_operon_type in enumerate(labels):
		for j, j_operon_type in enumerate(labels):
			if i > j:
				x_start, y_start = operons[i_operon_type]['operon_centroid_PaCMAP']
				x_end, y_end = operons[j_operon_type]['operon_centroid_PaCMAP']

				data['x'].append([x_start, x_end])
				data['y'].append([y_start, y_end])
				data['color'].append('black')
				data['alpha'].append(round(alpha_matrix[i][j], 1))

	tooltips = [('Relative distance/alpha', '@alpha')]

	return tooltips, data

#@Timer("create_avg_operons_clusters_network", text=lambda secs: f"Time to create avg operons clusters network: {format_timespan(secs)}")
def create_avg_operons_clusters_network(operons, operons_scatter, scatter_data, max_family_freq, min_family_freq):

	p_tooltips, p_data = scatter_data

	similarity_matrix, operons_labels = get_avgoperons_distance_matrix(operons)
	similarity_matrix = normalize_matrix(similarity_matrix)

	edge_tooltips, edge_data = create_edges_data(p_data, operons, similarity_matrix, operons_labels)

	p = figure(title = 'Genomic context types/clusters similarity network', plot_width=operons_scatter.plot_width, 
			   plot_height=operons_scatter.height, x_range = operons_scatter.x_range, y_range = operons_scatter.y_range)
	p.add_layout(Legend(orientation="horizontal"), 'above')

#	 p.circle('x', 'y', size='size', line_color='edgecolor', fill_color='facecolor', legend_field='type', alpha=1, source = p_data)
	p.multi_line('x', 'y', color='color', alpha='alpha', source=edge_data, name='edges')
	p.circle('avg_x', 'avg_y', size='node_size', line_color='edgecolor', fill_color='facecolor', alpha=1, source = p_data, name='nodes')

	p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
	p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
	p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
	p.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
	p.yaxis.axis_line_width = 0
	p.xaxis.axis_line_width = 0
	# define general features
	p.grid.visible = False
#	 p.outline_line_width = 0

	p.xaxis.axis_label = ""
	p.yaxis.axis_label = ""

	p.add_tools(HoverTool(tooltips=p_tooltips, names=['nodes']))
	p.add_tools(HoverTool(tooltips=edge_tooltips, names=['edges']))
	p.add_tools(LassoSelectTool())

	p.background_fill_color = "lightgrey"
	p.background_fill_alpha = 0.2

	p.legend.click_policy="hide"

	return p, similarity_matrix
	
#@Timer("create_clans_map_scatter", text=lambda secs: f"Time to create clans map scatter: {format_timespan(secs)}")
def create_clans_map_scatter(scatter_data):

	p_tooltips, p_data = scatter_data

	p = figure(title = 'Sequence similarity cluster (CLANS) map', plot_width=500, plot_height=500)
	p.add_layout(Legend(orientation="horizontal"), 'above')

	p.circle('clans_x', 'clans_y', size='size', line_color='edgecolor', fill_color='facecolor', alpha=1, source = p_data)

	p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
	p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
	p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
	p.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
	p.yaxis.axis_line_width = 0
	p.xaxis.axis_line_width = 0
	# define general features
	p.grid.visible = False
#	 p.outline_line_width = 0

	p.xaxis.axis_label = ""
	p.yaxis.axis_label = ""

	p.add_tools(HoverTool(tooltips=p_tooltips))
	p.add_tools(LassoSelectTool())

	p.background_fill_color = "lightgrey"
	p.background_fill_alpha = 0.2

#	 p.legend.click_policy="hide"

	return p

#@Timer("create_targets_table_data", text=lambda secs: f"Time to create targets table data: {format_timespan(secs)}")
def create_targets_table_data(all_syntenies, taxonomy, clusters_colors):

	data = {'Target EntrezID': [],
			'Superkingdom': [],
			'Phylum': [],
			'Class': [],
			'Order': [],
			'Genus': [],
			'Species': [],
			'Genomic context type': [],
			'color': []}

	for superkingdom in taxonomy.keys():
		for phylum in taxonomy[superkingdom].keys():
			for taxclass in taxonomy[superkingdom][phylum].keys():
				for order in taxonomy[superkingdom][phylum][taxclass].keys():
					for genus in taxonomy[superkingdom][phylum][taxclass][order].keys():
						for species in taxonomy[superkingdom][phylum][taxclass][order][genus].keys():
							for target in taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
								operon_type = all_syntenies[target]['operon_type']
								operon_type = 'GC Type {:05d}'.format(operon_type)

								if '-' not in operon_type:
									operon_color = clusters_colors[operon_type]['Color (RGBA)']

									data['Target EntrezID'].append(target)
									data['Superkingdom'].append(superkingdom)
									data['Phylum'].append(phylum)
									data['Class'].append(taxclass)
									data['Order'].append(order)
									data['Genus'].append(genus)
									data['Species'].append(species)
									data['Genomic context type'].append(operon_type.split()[-1])
									data['color'].append(operon_color)

	columns = [TableColumn(field=i, title=i) for i in data.keys()]

	columns = [TableColumn(field=i, title=i) if i not in ['color']
		  else TableColumn(field=i, title='Genomic context color', formatter=HTMLTemplateFormatter(template='<span style="color:<%= value %>;font-size:18pt;text-shadow: 1px 1px 2px #000000;">&#9632;</span>'))
		  for i in data.keys()]																			 

	data = pd.DataFrame(data)

	return data, columns

#@Timer("create_targets_summary_table", text=lambda secs: f"Time to create targets summary table: {format_timespan(secs)}")
def create_targets_summary_table(all_syntenies, taxonomy, clusters_colors):

	t_data, t_columns = create_targets_table_data(all_syntenies, taxonomy, clusters_colors)	
	t = DataTable(source=ColumnDataSource(t_data), columns=t_columns, width=1500, height=500)
	return t

#@Timer("compute_families_frequencies", text=lambda secs: f"Time to compute families frequencies: {format_timespan(secs)}")
def compute_families_frequencies(operons):

	families_frequencies = {}
	number_of_operons = 0
	for operon_type in operons:
		for i, target in enumerate(operons[operon_type]['target_members']):
			number_of_operons += 1
			for family in operons[operon_type]['operon_protein_families_structure'][i]:
				if family not in families_frequencies:
					families_frequencies[family] = 1
				else:
					families_frequencies[family] += 1

	return families_frequencies, number_of_operons

#@Timer("compute_families_frequencies_per_operon_cluster", text=lambda secs: f"Time to compute families frequencies per operon cluster: {format_timespan(secs)}")
def compute_families_frequencies_per_operon_cluster(operons):

	families_frequencies = {}
	for operon_type in operons:
		curr_frequencies, curr_number_of_operons = compute_families_frequencies({operon_type: operons[operon_type]})
		families_frequencies[operon_type] = curr_frequencies

		for family in families_frequencies[operon_type]:
			families_frequencies[operon_type][family] = families_frequencies[operon_type][family]/curr_number_of_operons

	return families_frequencies

#@Timer("clean_uncommon_families", text=lambda secs: f"Time to clean uncommon families: {format_timespan(secs)}")
def clean_uncommon_families(families_frequencies, families_summary, operons_labels, min_freq = 0.05):

	families_labels = [i for i in sorted(families_summary.keys()) if i < 10000 and i > 0]
	matrix = [[0 for family in families_labels] for operon_type in operons_labels]

	for i, operon_type in enumerate(operons_labels):
		for j, family in enumerate(families_labels):
			if family in families_frequencies[operon_type]:
				matrix[i][j] = families_frequencies[operon_type][family]

	matrix = pd.DataFrame(matrix, index=operons_labels, columns=families_labels).T
	matrix = matrix.loc[matrix.max(axis=1) > min_freq]
	matrix['sum'] = matrix.sum(axis=1)
	matrix = matrix.sort_values(by='sum', ascending=False)
	matrix = matrix.drop('sum', axis=1).T

	return matrix  

#@Timer("create_family_spectrum_data", text=lambda secs: f"Time to create family spectrum data: {format_timespan(secs)}")
def create_family_spectrum_data(families_frequencies_matrix, families_summary, family_colors, oprn_den_data, reference_family, dx = 1):

	data = {'xs': [],
			'ys': [],
			'edgecolor': [],
			'facecolor': [],
			'text': [],
			'text_x': [],
			'text_y': [],
			'tm_text_x': [],
			'tm_text_y': [],
			'tm_text': [],
			'tm_type': [],
			'tm_mode': [],
			'family': [],
			'found_models': [],
			'model_links': [],
			'keywords': [],
			'go_terms': [],
			'function': []}

	yyticklabels = oprn_den_data['leaf_label']
	yys = []
	y_step = oprn_den_data['y'][1] - oprn_den_data['y'][0]

	for i, operon_type in enumerate(families_frequencies_matrix.index):
		curr_y = oprn_den_data['y'][i]-(y_step/2)
		
		for j, family in enumerate(families_frequencies_matrix.columns):
			dy = families_frequencies_matrix.at[operon_type,family]
			if dy > 1:
				dy = 1
				
			dy = dy*y_step
			curr_x = j
			
			if family == reference_family:
				data['text'].append('Target protein: {}'.format(families_summary[family]['name']))
			elif family in families_summary:
				data['text'].append(families_summary[family]['name'])

			data['family'].append(family)

			if 'model_state' in families_summary[family]:
				model_state = families_summary[family]['model_state']

				if model_state == 'Model exists':
					model_state = families_summary[family]['structure']
				elif model_state == 'Model does not exist':
					model_state = 'click to model/view with Swiss-Model'
				else:
					model_state = 'Not possible to find'
			else:
				model_state = ''

			data['found_models'].append(model_state)

			if 'function' in families_summary[family]:
				if 'TM_topology' in families_summary[family]['function']:
					tm_type = families_summary[family]['function']["TM_topology"]
					keywords = ', '.join(sorted(families_summary[family]['function']['Keywords']))
					go_terms = '; '.join(sorted(families_summary[family]['function']['GO_terms']))
					function = families_summary[family]['function']['Function_description']

					if len(tm_type) > 0:
						tm_text = 'TM'
						tm_mode = 'Yes -> type:'
					else:
						tm_text = ''
						tm_mode = 'No'
				else:
					tm_type = ''
					tm_text = ''
					tm_mode = ''
					keywords = ''
					go_terms = ''   
					function = ''
			else:
				tm_type = 'n.a.'
				tm_text = ''
				tm_mode = 'n.a.'
				keywords = 'n.a.'
				go_terms = 'n.a.'   
				function = 'n.a.'

			if 'structure' in families_summary[family]:
				structure = families_summary[family]['structure']
				if structure == '':
					uniprot_code = families_summary[family]['uniprot_code']
					structure = 'https://swissmodel.expasy.org/repository/uniprot/{}'.format(uniprot_code)
			else:
				structure = 'n.a.'

			data['model_links'].append(structure)

			data['facecolor'].append(family_colors[family]['Color (RGBA)'])
			data['edgecolor'].append(family_colors[family]['Line color'])
			data['xs'].append([curr_x, curr_x, curr_x+dx,curr_x+dx])
			data['ys'].append([curr_y, curr_y+dy, curr_y+dy, curr_y])
			data['text_x'].append(((curr_x+dx)/2)+curr_x)
			data['text_y'].append(curr_y+0.5)

			data['tm_text_x'].append(((curr_x+dx)/2)+curr_x)
			data['tm_text_y'].append(curr_y+0.25)
			data['tm_text'].append(tm_text)
			data['tm_type'].append(tm_type)
			data['tm_mode'].append(tm_mode)

			data['go_terms'].append(go_terms)
			data['keywords'].append(keywords)
			data['function'].append(function)

		curr_y -= 1
		
		yys.append(curr_y+(y_step/2))

	yyticklabels = {yys[i]: yyticklabels[i] for i in range(len(yyticklabels))}

	tooltips = [('Protein family', '@text'),
				('Protein family code', '@family'),
				('Structure', '@found_models'),
				('Predicted membrane protein', '@tm_mode @tm_type'),
				('Keywords', '@keywords'),
				('GO terms', '@go_terms'),
				('Function', '@function')]

	return tooltips, data, yyticklabels

#@Timer("create_family_frequency_per_operon_figure", text=lambda secs: f"Time to create family frequency per operon figure: {format_timespan(secs)}")
def create_family_frequency_per_operon_figure(operons, families_summary, family_colors, oprn_dendogram, oprn_den_data, reference_family, height_factor = 25*1.2, min_freq = 0.05):

	families_frequencies = compute_families_frequencies_per_operon_cluster(operons)
	families_frequencies_matrix = clean_uncommon_families(families_frequencies, families_summary, operons_labels = oprn_den_data['leaf_label'], min_freq = min_freq)

	p_tooltips, p_data, p_yyticklabels = create_family_spectrum_data(families_frequencies_matrix, families_summary, family_colors, oprn_den_data, reference_family = reference_family)
		
	p = figure(plot_width=1250, plot_height=oprn_dendogram.height, x_range = [0, len(families_frequencies_matrix.columns)], y_range = oprn_dendogram.y_range, 
			   toolbar_location="left", title = 'Genomic context family spectrum (hover to get more information and click to model with SWISS-MODEL)') 
	
	p.patches('xs', 'ys', fill_color = 'facecolor', line_color = 'edgecolor', line_width = 1, source = p_data, 
			  hover_fill_color = 'white', hover_line_color = 'edgecolor', hover_fill_alpha = 0.5, 
			  selection_fill_color='facecolor', selection_line_color='edgecolor',
			  nonselection_fill_color='facecolor', nonselection_line_color='edgecolor', nonselection_fill_alpha=0.2)
		
	p.yaxis.ticker = list(p_yyticklabels.keys())
	p.yaxis.major_label_overrides = {int(i): p_yyticklabels[i] for i in p_yyticklabels.keys()}
	
	p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
	p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
	p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
	p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
	p.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
	p.yaxis.axis_line_width = 0
	p.xaxis.axis_line_width = 0
	# define general features
	p.grid.visible = False
	p.outline_line_width = 0

	p.add_tools(HoverTool(tooltips=p_tooltips))
	p.add_tools(TapTool(callback = OpenURL(url='@model_links')))

	return p

#@Timer("create_family_table_data", text=lambda secs: f"Time to create family table data: {format_timespan(secs)}")
def create_family_table_data(operons, families_summary):
		
	data = {'Family': [],
			'Name': [],
			'UniprotKB': [],
			'Frequency': [],
			'Keywords': [],
			'GO terms: Molecular function (F)': [],
			'GO terms: Cellular component (C)': [],
			'GO terms: Biological process (P)': [],
			'Function': [],
			'Predicted membrane protein': [],
			'Structure': []}
	
	families_frequencies, number_of_operons = compute_families_frequencies(operons)
	
	for family in families_frequencies:  
		if family > 0 and family < 10000:
			data['Family'].append(family)
			data['Frequency'].append(round(families_frequencies[family]*100/number_of_operons, 1))
			data['Name'].append(families_summary[family]['name'])
			uniprot_code = 'nan'
			
			if 'function' in families_summary[family]:
				uniprot_code = families_summary[family]['uniprot_code']
				if 'TM_topology' in families_summary[family]['function']:
					tm_type = families_summary[family]['function']["TM_topology"]
					keywords = ', '.join(sorted(families_summary[family]['function']['Keywords']))
					go_terms = sorted(families_summary[family]['function']['GO_terms'])
					function = families_summary[family]['function']['Function_description']

					if len(tm_type) > 0:
						tm_mode = 'Yes'
					else:
						tm_mode = 'No'
						
					f_go = '; '.join([go for go in go_terms if go.startswith('F:')])
					c_go = '; '.join([go for go in go_terms if go.startswith('C:')])
					p_go = '; '.join([go for go in go_terms if go.startswith('P:')])
			
				else:
					tm_mode = ''
					keywords = ''
					f_go = '' 
					c_go = '' 
					p_go = '' 
					function = ''
			else:
				tm_mode = 'n.a.'
				keywords = 'n.a.'
				f_go = 'n.a.' 
				c_go = 'n.a.' 
				p_go = 'n.a.'   
				function = 'n.a.'
			
			if 'structure' in families_summary[family]:
				structure = families_summary[family]['structure']
				if structure == '':
					structure = 'Click to model'
			else:
				structure = 'n.a.'
			
			data['UniprotKB'].append(uniprot_code)
			data['Keywords'].append(keywords)
			data['GO terms: Molecular function (F)'].append(f_go)
			data['GO terms: Cellular component (C)'].append(c_go)
			data['GO terms: Biological process (P)'].append(p_go)
			data['Function'].append(function)
			data['Predicted membrane protein'].append(tm_mode)
			data['Structure'].append(structure)
	
	columns = [TableColumn(field=i, title=i) if i not in ['Structure', 'UniprotKB']
		  else TableColumn(field=i, title=i, formatter=HTMLTemplateFormatter(template='<a href="https://swissmodel.expasy.org/repository/uniprot/<%= UniprotKB %>"><%= value %></a>')) if i=='Structure'
		  else TableColumn(field=i, title='Representative UniprotKB', formatter=HTMLTemplateFormatter(template='<a href="https://www.uniprot.org/uniprot/<%= value %>"><%= value %></a>'))  
		  for i in data.keys()]
	
	data = pd.DataFrame(data).sort_values(by='Frequency', ascending=False)
	
	return data, columns

#@Timer("create_families_frequency_table", text=lambda secs: f"Time to create families frequency table: {format_timespan(secs)}")
def create_families_frequency_table(operons, families_summary):
	
	t_data, t_columns = create_family_table_data(operons, families_summary)
	
	t = DataTable(source=ColumnDataSource(t_data), columns=t_columns, width=2250, height=300)
	
	div = Div(text="""<b>Represented families:</b>""") 
	
	return t, div


# 10. Routines to deal with Clans maps as input source of IDs

#@Timer("get_ncbicodes_order_in_clans", text=lambda secs: f"Time to get ncbicodes order in clans: {format_timespan(secs)}")
def get_ncbicodes_order_in_clans(clans_file):
   
	ncbids_ordered = []

	count = 0
	found_seq = False
	with open(clans_file, 'r') as clans:
		for line in clans:
			if '<seq>' in line:
				found_seq = True
			elif '</seq>' in line:
				found_seq = False
			elif found_seq and line.startswith('>'):
				line = line[1:]
				ncbi_code = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()
				ncbids_ordered.append(ncbi_code)
			   
	return ncbids_ordered

#@Timer("get_clusters_from_clans", text=lambda secs: f"Time to get clusters from clans: {format_timespan(secs)}")
def get_clusters_from_clans(clans_file, cluster_codes = None):

	ncbis_ordered = get_ncbicodes_order_in_clans(clans_file)

	clusters = {}
	
	if cluster_codes != None:
		found_seqgroup = False
		with open(clans_file, 'r') as in_clans:
			for line in in_clans:
				for cluster_code in cluster_codes:
					if '<seqgroup' in line:
						found_seqgroup = True
						found_allowed_cluster = False

					elif found_seqgroup:
						if 'name=' in line and cluster_code in line:
							current_cluster = line.split('=')[-1].strip()
							found_allowed_cluster = True
						   
						elif 'numbers=' in line and found_allowed_cluster:
							numbers = line.split('=')[-1].split(';')[:-1]
							numbers = [int(i) for i in numbers]
							numbers = [ncbis_ordered[i] for i in numbers]
							clusters[current_cluster] = numbers

							found_allowed_cluster = False

	else:
		label = clans_file.split('/')[-1].split('.')[0]
		clusters[label] = ncbis_ordered

	return clusters

# MAIN ROUTINES (corresponding to the major steps in the pipeline)

#@Timer("write_arguments_file", text=lambda secs: f"Time to write arguments file: {format_timespan(secs)}")
def write_arguments_file(args, out_label):

	out_file = '{}_input_arguments.log'.format(out_label)

	with open(out_file, 'w') as f:
		for key in args.__dict__:
			f.write('{}:\t{}\n'.format(key, args.__dict__[key]))


@Timer("parse_targets", text=lambda secs: f"Time to parse targets: {format_timespan(secs)}")
def parse_targets(targets, label = None, clans_patterns = None):

	targets_list = {}

	for target in targets:
		if os.path.isfile(target):

			if target.endswith('.clans'):
				targets_list = get_clusters_from_clans(target, cluster_codes = clans_patterns)

			else:
				curr_label = target.split('/')[-1].split('.')[0]
				if curr_label not in targets_list:
					targets_list[curr_label] = []

				is_fasta = False
				with open(target, 'r') as curr_infile:
					for line in curr_infile:
						if line.startswith('>'):
							is_fasta = True
							curr_target = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()

						elif not is_fasta:
							curr_target = line.strip().split()[0]

						if curr_target not in targets_list[curr_label]:
							targets_list[curr_label].append(curr_target)

		else:
			if label not in targets_list:
				targets_list[label] = []
		
			targets_list[label].append(target)

	print(' ... Found {} jobs to do: {}'.format(len(targets_list), list(targets_list.keys())))
	
	return targets_list

@Timer("download_and_parse_refseq_and_gb_databases", text=lambda secs: f"Time to download and parse refseq and gb databases: {format_timespan(secs)}")
def download_and_parse_refseq_and_gb_databases(databases = ['genbank', 'refseq']):
	
	database_assembly_mapping = {}

	for db in databases:
		print(' ... Taking care of {} summary table'.format(db))

		summary_table = '{}/assembly_summary_{}.txt'.format(os.getcwd(), db)

		if not os.path.isfile(summary_table):
			link = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/{}/assembly_summary_{}.txt'.format(db, db)
			print(' ... ... Downloading')
			t1 = Timer('db_download_single', text=lambda secs: f"Time for download {db}: {format_timespan(secs)}")
			t1.start()
			download_db_table = sp.Popen(['wget', link, '-O', summary_table], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
			stdout, stderr = download_db_table.communicate()
			t1.stop()
			print(' ... ... Done Downloading')

		else:
			print(' ... ... Summary table already exists (version from: {})'.format(time.ctime(os.path.getmtime(summary_table))))

		t2 = Timer('db_parsing_single', text=lambda secs: f"Time for parsing {db}: {format_timespan(secs)}")
		t2.start()
		with open(summary_table, 'r') as summary:
			print(' ... ... Parsing summary table')
			for line in summary:
				if not line.startswith('#'):
					data = line.strip().split('\t')
					database_assembly_mapping[data[0]] = data[19]
		t2.stop()
		print(' ... ... Done parsing')

	return database_assembly_mapping

#@Timer("get_genomic_context_information_for_ncbi_codes", text=lambda secs: f"Time to get genomic context information for ncbi codes: {format_timespan(secs)}")
def get_genomic_context_information_for_ncbi_codes(target_ncbi_codes, refseq_gb_assembly_map = None, n_flanking5 = None, n_flanking3 = None, exclude_partial = None, tmp_folder = None, threads = None):

	# Prepare all parallel jobs
	separated_jobs = chunk_list(target_ncbi_codes, threads)

	list_arguments = [i for i in zip(separated_jobs, [refseq_gb_assembly_map for job in separated_jobs], [n_flanking5 for job in separated_jobs], [n_flanking3 for job in separated_jobs], [exclude_partial for job in separated_jobs], [tmp_folder for job in separated_jobs], range(threads))]

	pool = ThreadPool(threads)
	results = pool.imap_unordered(get_genomic_contexts_for_ncbi_codes, list_arguments)

	all_syntenies = {key: dic[key] for dic in results for key in dic.keys()}
	
	return all_syntenies  

#@Timer("get_genomic_contexts_for_ncbi_codes", text=lambda secs: f"Time to get genomic contexts for ncbi codes: {format_timespan(secs)}")
def get_genomic_contexts_for_ncbi_codes(arguments):

	target_ncbi_codes = arguments[0]
	refseq_gb_assembly_map = arguments[1]
	n_flanking5 = arguments[2]
	n_flanking3 = arguments[3]
	exclude_partial = arguments[4]
	tmp_folder = arguments[5]
	thread_id = arguments[6]

	out_json = '{}_genomic_context_information.json'.format(thread_id)
	if os.path.isfile(out_json):
		with open(out_json, 'r') as fp:	 # allows to restart the searches without having to get all the info again in case it crashes
			all_syntenies = json.load(fp)
			starting_i = target_ncbi_codes.index(list(all_syntenies.keys())[-1])+1
			print(' ... Thread {}: A partial search was already performed for {} entrezIDs (from which, {} were valid). Will continue it.\n'.format(thread_id, starting_i+1, len(all_syntenies)))
	else:
		all_syntenies = {}
		starting_i = 0

	t1 = Timer('collecting_genomic_context', text=lambda secs: f"Time to collect all genomic context: {format_timespan(secs)}")
	t1.start()
	for i in range(starting_i, len(target_ncbi_codes)):
		curr_target_code = target_ncbi_codes[i].strip()
		curr_target_count = i+1

		t2 = Timer('find_ncbi_code_assembly',
				   text=lambda secs: f"Time to find assembly ncbi code for {curr_target_code}: {format_timespan(secs)}")
		t2.start()
		ncbi_code, assembly_id, assembly_link = find_ncbi_code_assembly(curr_target_code, refseq_gb_assembly_map)
		t2.stop()

		if assembly_id != 'nan':
			print(" ... {} belongs to assembly {} ({}/{})".format(curr_target_code, assembly_id, curr_target_count, len(target_ncbi_codes)))
			t3 = Timer('download_and_extract_assembly', text=lambda secs: f"Time to download and extract assembly {assembly_id} ({curr_target_code}): {format_timespan(secs)}")
			t3.start()
			assembly = download_and_extract_assembly(assembly_id, assembly_link, tmp_folder = tmp_folder, label = curr_target_code, get_chunk = True, chunk_size = ((n_flanking5+n_flanking3)*2)+1, target = ncbi_code)
			t3.stop()

			if assembly != 'nan':
				t4 = Timer('extract_flanking_genes',
						   text=lambda millisecs: f"Time to extract {n_flanking5 + n_flanking3} flanking genes for {curr_target_code} from {assembly_id}: {format_timespan(millisecs)}")
				t4.start()
				flanking_genes = get_n_flanking_genes(ncbi_code, assembly, n_5 = n_flanking5, n_3 = n_flanking3, exclude_partial = exclude_partial)
				t4.stop()
				if flanking_genes != 'nan':
					t5 = Timer('add_sequences_to_flanking_genes',
							   text=lambda
								   millisecs: f"Time to add_sequences_to_flanking_genes for {curr_target_code} from {assembly_id}: {format_timespan(millisecs)}")
					t5.start()
					flanking_genes = add_sequences_to_flanking_genes(flanking_genes, ncbi_code)
					t5.stop()

					print(' ... ... Species: {}'.format(flanking_genes['species'])) 

					if curr_target_code not in all_syntenies:
						all_syntenies[curr_target_code] = {}
					
					all_syntenies[curr_target_code]['flanking_genes'] = flanking_genes
					all_syntenies[curr_target_code]['assembly_id'] = [ncbi_code, assembly_id, assembly_link]
					all_syntenies[curr_target_code]['species'] = all_syntenies[curr_target_code]['flanking_genes']['species']
		else:
			print("\n ... > There is no assembly for {} ({}/{})\n".format(curr_target_code, curr_target_count, len(target_ncbi_codes)))

		with open(out_json, 'w') as fp:
			json.dump(all_syntenies, fp, indent=4)
	t1.stop()

	return all_syntenies  

#@Timer("find_and_add_protein_families", text=lambda secs: f"Time to find and add protein families: {format_timespan(secs)}")
def find_and_add_protein_families(in_syntenies, out_label = None, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, blast = None, mmseqs = None, min_coverage = None, default_base = None, tmp_folder = None, method = None):
	t1 = Timer(f"find_and_add_protein_families_{method}",
		   text=lambda secs: f"Time to find and add protein families using {method}: {format_timespan(secs)}")
	t1.start()
	if os.path.isfile('all_syntenies.json') and os.path.isfile('protein_families_summary.json'):

		in_syntenies = json.load(open('all_syntenies.json', 'r'))
		protein_families = json.load(open('protein_families_summary.json','r'))
		protein_families = {int(family): protein_families[family] for family in protein_families}

		print(' ... A job was previously run. Found {} families.'.format(len(protein_families)))

	else:
		print(' ... Doing all against all searches with {}'.format(method))

		t2 = Timer(f"compute_all_agains_all_distance_matrix_{method}",
				   text=lambda secs: f"Time to compute all agains all distance matrix using {method}: {format_timespan(secs)}")
		t2.start()
		distance_matrix, ordered_ncbi_codes = compute_all_agains_all_distance_matrix(in_syntenies, out_label = out_label, num_threads = num_threads, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, min_coverage = min_coverage, method = method, mmseqs = mmseqs, blast = blast, default_base = default_base, tmp_folder = tmp_folder)	
		t2.stop()

		t3 = Timer('find_clusters_in_distance_matrix',
				   text=lambda secs: f"Time to find clusters in distance matrix using {method}: {format_timespan(secs)}")
		t3.start()
		protein_clusters = find_clusters_in_distance_matrix(distance_matrix)
		t3.stop()

		t4 = Timer('mask_singleton_clusters',
				   text=lambda secs: f"Time to mask_singleton_clusters using {method}: {format_timespan(secs)}")
		t4.start()
		protein_clusters = mask_singleton_clusters(protein_clusters)
		t4.stop()

		print(' ... Assigning families')
		# define the correct numbering of the families so that the reference family is the largest, pseudogenes are > 10000 and all the others
		# are continuous
		t5 = Timer('assign_families',
				   text=lambda secs: f"Time to assign families: {format_timespan(secs)}")
		t5.start()

		curr_numbers = []
		for target in in_syntenies:
			in_syntenies[target]['flanking_genes']['families'] = []
			for i, ncbi_code in enumerate(in_syntenies[target]['flanking_genes']['ncbi_codes']):
				protein_name = in_syntenies[target]['flanking_genes']['names'][i]
				try:
					protein_family = protein_clusters[ordered_ncbi_codes.index(ncbi_code)]
				except:
					protein_family = 10000
				
				if protein_name == 'pseudogene':
					protein_family = 10000
				if ncbi_code == target:
					protein_family = max(protein_clusters)+1

				in_syntenies[target]['flanking_genes']['families'].append(protein_family)
				curr_numbers.append(protein_family)

				if ncbi_code == in_syntenies[target]['assembly_id'][0]:
					in_syntenies[target]['target_family'] = protein_family

		curr_numbers = sorted(list(set(curr_numbers)))
		for target in in_syntenies:
			for i, ncbi_code in enumerate(in_syntenies[target]['flanking_genes']['ncbi_codes']):
				protein_name = in_syntenies[target]['flanking_genes']['names'][i]
				protein_family = in_syntenies[target]['flanking_genes']['families'][i]

				if protein_family <= max(protein_clusters):
					if protein_family != range(min(protein_clusters), max(protein_clusters)+1)[curr_numbers.index(protein_family)]:
						protein_family = range(min(protein_clusters), max(protein_clusters)+1)[curr_numbers.index(protein_family)]
						in_syntenies[target]['flanking_genes']['families'][i] = protein_family

			in_syntenies[target]['target_family'] = in_syntenies[target]['flanking_genes']['families'][in_syntenies[target]['flanking_genes']['ncbi_codes'].index(in_syntenies[target]['assembly_id'][0])]
		t5.stop()

		t6 = Timer('get_protein_families_summary',
				   text=lambda secs: f"Time to get protein families summary: {format_timespan(secs)}")
		t6.start()
		protein_families = get_protein_families_summary(in_syntenies, write_to_file = True, out_label = out_label)
		t6.stop()

	t1.stop()

	return in_syntenies, protein_families

#@Timer("update_families_with_functions_and_structures", text=lambda secs: f"Time to update families with functions and structures: {format_timespan(secs)}")
def update_families_with_functions_and_structures(protein_families_summary, get_pdb = None, get_functional_annotations = None, threads = 1):

	# Prepare all parallel jobs
	separated_jobs = chunk_list(sorted(list(protein_families_summary.keys())), threads)

	list_arguments = [i for i in zip(separated_jobs, [protein_families_summary for job in separated_jobs], [get_pdb for job in separated_jobs], [get_functional_annotations for job in separated_jobs], range(threads))]

	pool = ThreadPool(threads)
	results = pool.imap_unordered(add_functions_and_structures_to_families, list_arguments)

	protein_families_summary = {key: dic[key] for dic in results for key in dic.keys()}

	json.dump(protein_families_summary, open('protein_families_summary.json', 'w'), indent = 4)

	return protein_families_summary

#@Timer("add_functions_and_structures_to_families", text=lambda secs: f"Time to add functions and structures to families: {format_timespan(secs)}")
def add_functions_and_structures_to_families(arguments):

	targets = arguments[0]
	protein_families_summary = arguments[1]
	get_pdb = arguments[2]
	get_functional_annotations = arguments[3]
	job_id = arguments[4]

	protein_families_summary = {key.item(): protein_families_summary[key.item()] for key in targets}

	all_members = [i for family in protein_families_summary for i in protein_families_summary[family]['members'] if family > 0 and family < 10000]

	print(' ... Thread {}: Mapping {} members'.format(job_id, len(all_members)))
	family_uniprot_codes = map_codes_to_uniprot(all_members)

	t1 = Timer('mapping_all_family_members', text=lambda secs: f"Time to map {len(all_members)} family members: {format_timespan(secs)}")
	t1.start()
	for family in sorted(list(protein_families_summary.keys())):
		t2 = Timer('map_one_familiy', text=lambda secs: f"Time to map one single family: {format_timespan(secs)}")
		t2.start()
		if 'function' not in protein_families_summary[family]:
			if family > 0 and family < 10000:

				print(' ... Thread {}: Family {} ({} members)'.format(job_id, family, len(protein_families_summary[family]['members'])))

				family_function = ''
				family_structure = ''
				family_uniprot_code = ''
				family_model_state = 'Model does not exist'

				for target in sorted(protein_families_summary[family]['members']):
					uniprot_code = family_uniprot_codes[target]

					if family_uniprot_code == '' or family_structure == '' or (family_function == '' or family_function['Function_description'] == ''):

						if get_pdb and family_structure == '':
							t3 = Timer('find_uniprot_in_swiss_model_repository',
									   text=lambda secs: f"Time to find {uniprot_code} in swiss model repository: {format_timespan(secs)}")
							t3.start()
							curr_pdb = find_uniprot_in_swiss_model_repository(uniprot_code)
							t3.stop()
							if 'nan' in curr_pdb:
								t4 = Timer('find_uniprot_in_alphafold_database',
										   text=lambda
											   secs: f"Time to find {uniprot_code} in alphafold database: {format_timespan(secs)}")
								t4.start()
								curr_pdb = find_uniprot_in_alphafold_database(uniprot_code)
								t4.stop()

							if 'nan' not in curr_pdb:
								family_structure = curr_pdb
								family_uniprot_code = uniprot_code
								family_model_state = 'Model exists'

							elif uniprot_code != 'nan' and curr_pdb != 'nan*':
								family_uniprot_code = uniprot_code

						if get_functional_annotations and (family_function == '' or family_function['Function_description'] == ''):
							t5 = Timer('get_uniprot_annotations',
									   text=lambda
										   secs: f"Time to get uniprot annotation for {uniprot_code}: {format_timespan(secs)}")
							t5.start()
							curr_uniprot_annotations = get_uniprot_annotations(uniprot_code, previous_annotations = family_function)
							t5.stop()

							if curr_uniprot_annotations != 'nan':
								family_function = curr_uniprot_annotations

							if uniprot_code != 'nan' and family_structure == '' and family_uniprot_code == '':
								family_uniprot_code = uniprot_code

				if family_uniprot_code == '':
					family_model_state = 'Not possible to map'

			else:
				family_function = 'n.a.'
				family_structure = 'n.a.'
				family_uniprot_code = 'n.a.'
				family_model_state = 'n.a.'

			protein_families_summary[family]['function'] = family_function
			protein_families_summary[family]['structure'] = family_structure
			protein_families_summary[family]['uniprot_code'] = family_uniprot_code
			protein_families_summary[family]['model_state'] = family_model_state

			json.dump(protein_families_summary, open('{}_protein_families_summary.json'.format(job_id), 'w'), indent = 4)
		t2.stop()
	t1.stop()
	return protein_families_summary

#@Timer("find_and_add_operon_types", text=lambda secs: f"Time to find and add operon types: {format_timespan(secs)}")
def find_and_add_operon_types(in_syntenies, protein_families_summary, label = None, advanced = False, min_family_freq = None, max_family_freq = None):

	print(' ... Using mode Advanced?', advanced)

	if len(in_syntenies) > 1:
		if advanced:
			t1 = Timer('find_operon_types_advanced',
					   text=lambda secs: f"Time to find operon types in advanced mode: {format_timespan(secs)}")
			t1.start()
			# get the clusters by excluding the most common families
			clean_coordinates, operon_clusters, ordered_ncbi_codes = find_operon_clusters_with_PaCMAP(in_syntenies, protein_families_summary, clean = True, min_freq = min_family_freq, max_freq = max_family_freq)
			# and now all PacMap coordinates by assuming all families. This will be later used for sorting the dendogram
			all_coordinates, ordered_ncbi_codes = find_operon_clusters_with_PaCMAP(in_syntenies, protein_families_summary, clean = False, coordinates_only = True)
			t1.stop()
		else:
			t2 = Timer('find_operon_types_standard',
					   text=lambda secs: f"Time to find operon types in standard mode: {format_timespan(secs)}")
			t2.start()
			distance_matrix, ordered_ncbi_codes = compute_operons_distance_matrix(in_syntenies, label = label)
			operon_clusters = find_clusters_in_distance_matrix(distance_matrix, t = 0.2)
			coordinates = [[np.nan, np.nan] for i in operon_clusters]
			t2.stop()
	else:
		operon_clusters = [1]
		ordered_ncbi_codes = list(in_syntenies.keys())
		coordinates = [[np.nan, np.nan] for i in operon_clusters]
	
	for i, target in enumerate(ordered_ncbi_codes):
		in_syntenies[target]['operon_type'] = int(operon_clusters[i])

		if advanced:
			in_syntenies[target]['operon_filtered_PaCMAP'] = list([float(a) for a in clean_coordinates[i]])
			in_syntenies[target]['operon_PaCMAP'] = list([float(a) for a in all_coordinates[i]])

	operon_types = get_operon_types_summary(in_syntenies, write_to_file = True, label = label)
	
	return in_syntenies, operon_types

#@Timer("annotate_TMs_in_all", text=lambda secs: f"Time to annotate TMs in all: {format_timespan(secs)}")
def annotate_TMs_in_all(in_syntenies, annotation_TM_mode, annotation_TM_file, label = None):

	out_dir = '{}/{}_TM_annotations'.format(os.getcwd(), label)
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	in_fasta, seqs_lens = write_flanking_sequences_to_fasta(in_syntenies, out_dir, label, exclude_pseudogenes = True)

	if annotation_TM_file == None:
		annotation_TM_file = run_TM_signal_peptide_annotation(in_fasta, annotation_TM_mode = annotation_TM_mode)
	else:
		annotation_TM_file = '{}/{}'.format(out_dir, annotation_TM_file.split('/')[-1])
	
	if os.path.isfile(annotation_TM_file):
		protein_annotations = parse_annotation_TM_file(annotation_TM_file, annotation_TM_mode)
		in_syntenies = add_TM_annotations_to_flanking_genes(in_syntenies, protein_annotations)
	else:
		print(' ... WARNING: It was not possible to annotate TM and signal peptides to collected protein sequences.\n')
		print('	 What to do now? There are 2 options:')
		print('	 (1) Check if {} is in your path'.format(annotation_TM_mode))
		print('	 (2) Use {} as input over the {} webserver (follow instructions bellow) and run GCsnap again using the flag -annotation_TM_file and the results as input in a txt file'.format(in_fasta, annotation_TM_mode))

		if annotation_TM_mode == 'phobius':
			print('		  Phobius webserver: http://phobius.sbc.su.se/')
			print('		  Run mode: "short"')
			print('		  Output file: copy-paste output text below the line to a .txt file')
			print('		  Save in folder: {}'.format(out_dir))

		elif annotation_TM_mode == 'tmhmm':
			print('		  TMHMM webserver: https://services.healthtech.dtu.dk/service.php?TMHMM-2.0')
			print('		  Run mode: "One line per protein"')
			print('		  Output file: copy-paste output text below the line to a .txt file')
			print('		  Save in folder: {}'.format(out_dir))
	return in_syntenies

#@Timer("make_genomic_context_figure", text=lambda secs: f"Time to make genomic context figure: {format_timespan(secs)}")
def make_genomic_context_figure(operons, most_populated_operon, all_syntenies, families_summary, cmap = None, label = None, out_format = None):

	# define the reference family as the one of the target in the most populated operon type
	reference_family = all_syntenies[operons[most_populated_operon]['target_members'][0]]['target_family']
	family_colors = define_family_colors(list(families_summary.keys()), reference_family, mode = 'matplotlib', cmap = cmap)

	draw_genomic_context(operons, all_syntenies, family_colors, reference_family, label = label, out_format = out_format)
	draw_genomic_context_legend(families_summary, family_colors, reference_family, label = label, out_format = out_format)

#@Timer("make_interactive_genomic_context_figure", text=lambda secs: f"Time to make interactive genomic context figure: {format_timespan(secs)}")
def make_interactive_genomic_context_figure(operons, all_syntenies, families_summary, taxonomy, most_populated_operon, tree = None, sort_mode = None, input_targets = None, gc_legend_mode = None, cmap = None, label = None, min_coocc = None, n_flanking5=None, n_flanking3=None, tree_format=None):

	# define the reference family as the one of the target in the most populated operon type
	reference_family = all_syntenies[operons[most_populated_operon]['target_members'][0]]['target_family']
	family_colors = define_family_colors(list(families_summary.keys()), reference_family, mode = 'bokeh', cmap = cmap)

	# output to static HTML file
	output_file("{}/{}_interactive_output.html".format(os.getcwd(), label))

	# Work on most conserved genomic context figure
	most_common_gc_figure = create_most_common_genomic_features_figure(operons, all_syntenies, families_summary, reference_family = reference_family, family_colors = family_colors, n_flanking5=n_flanking5, n_flanking3=n_flanking3)

	# Work on gene co-occurence figure
	coocurrence_figure, graph_coord = create_graph_figure(operons, reference_family, families_summary, family_colors, most_common_gc_figure, min_coocc = min_coocc, mode = 'coocurrence')
	adjacency_figure, graph_coord = create_graph_figure(operons, reference_family, families_summary, family_colors, most_common_gc_figure, min_coocc = min_coocc, mode = 'adjacency', graph_coord = graph_coord, previous_net = coocurrence_figure)

	# Work on dendogram for the genomic context block
	if tree != None:
		input_targets = [target for operon in operons for target in operons[operon]['target_members']]
	syn_dendogram, syn_den_data = make_dendogram_figure(taxonomy, operons, tree = tree, mode = sort_mode, input_targets = input_targets, show_leafs = False, height_factor = 25*1.2, tree_format = tree_format)

	# Work on the genomic context block
	genomic_context_figure = create_genomic_context_figure(operons, all_syntenies, family_colors, syn_den_data, syn_dendogram, most_common_gc_figure, reference_family, legend_mode = gc_legend_mode, height_factor = 25*1.2)

	# Work on the legend figure
	legend_figure = create_legend_figure(operons, families_summary, reference_family, family_colors, genomic_context_figure)

	# Make the tabs for the network and most common genomic context
	tab1 = Panel(child=coocurrence_figure, title='Gene co-occurrence network')
	tab2 = Panel(child=adjacency_figure, title='Gene adjacency network')
	tab3 = Panel(child=most_common_gc_figure, title='Most common genomic features')
	tabs = Tabs(tabs=[tab1, tab2, tab3])

	# Make a grid out of them
	grid = gridplot([[None, tabs, None], 
					 [syn_dendogram, genomic_context_figure, legend_figure]], merge_tools = True)

	save(grid)

#@Timer("make_advanced_interactive_genomic_context_figures", text=lambda secs: f"Time to make advanced interactive genomic context figures: {format_timespan(secs)}")
def make_advanced_interactive_genomic_context_figures(operons, all_syntenies, families_summary, taxonomy, most_populated_operon, tree = None, sort_mode = None, input_targets = None, gc_legend_mode = None, cmap = None, label = None, min_coocc = None, n_flanking5=None, n_flanking3=None, tree_format=None, max_family_freq=None, min_family_freq=None, min_family_freq_accross_contexts=None, clans_file=None, out_label = None, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, min_coverage = None, method = None, mmseqs = None, blast = None, default_base = None, tmp_folder = None):

	# define the reference family as the one of the target in the most populated operon type
	reference_family = all_syntenies[operons[most_populated_operon]['target_members'][0]]['target_family']
	family_colors = define_family_colors(list(families_summary.keys()), reference_family, mode = 'bokeh', cmap = cmap)

	# make first the summary page

	print(' ... Making summary page\n')

	# Create the dataframe with all data used to plot. it allows for cross interativity between plots
	clusters_colors = define_operons_colors(operons.keys(), mode = 'bokeh', cmap='gist_rainbow')
	scatter_data = create_data_structure(operons, clans_file, clusters_colors=clusters_colors, all_syntenies=all_syntenies, out_label = out_label, num_threads = num_threads, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, min_coverage = min_coverage, method = method, mmseqs = mmseqs, blast = blast, default_base = default_base, tmp_folder = tmp_folder)

	# Plot the scatter of the operons PaCMAP coordinates
	operons_scatter = create_operons_clusters_scatter(scatter_data)
	# Create network of operon types similarities (based on the centroids in PaCMAP space and the minimum distance between cluster members)
	operons_network, operons_distance_matrix = create_avg_operons_clusters_network(operons, operons_scatter, scatter_data, max_family_freq, min_family_freq)

	# Plot the CLANS map of the input target sequences if given
	clans_scatter = create_clans_map_scatter(scatter_data)
	scatter_row = gridplot([[operons_scatter, operons_network, clans_scatter]], merge_tools = True)

	# Now create two tabs
	all_tabs = []

	# 1. The tab of species (i.e. a table where each line corresponds to an input target and includes the taxonomy information as well as the cluster type it belongs to)
	targets_table = create_targets_summary_table(all_syntenies, taxonomy, clusters_colors)
	all_tabs.append(Panel(child = targets_table, title = 'Input targets summary'))

	# 2. The tab the shows the gene composition of each cluster type, sorted by similarity and connected by a cladogram build from the 
	# distance matrix used to build the network above

	# Make the dendogram
	oprn_dendogram, oprn_den_data = make_dendogram_figure(False, operons, show_leafs = False, height_factor = 25*1.2, mode = 'operon clusters', distance_matrix=(1-operons_distance_matrix), labels=sorted([i for i in operons if '-' not in i]), colors=clusters_colors)
	family_freq_figure = create_family_frequency_per_operon_figure(operons, families_summary, family_colors, oprn_dendogram, oprn_den_data, reference_family, height_factor = 25*1.2, min_freq = min_family_freq_accross_contexts/100)
	gc_row = gridplot([[oprn_dendogram, family_freq_figure]], merge_tools = True)
	all_tabs.append(Panel(child = gc_row, title = 'Genomic contexts clusters hierarchy'))

	all_tabs = Tabs(tabs=all_tabs)
	grid = gridplot([[scatter_row],[all_tabs]], merge_tools = False)

	# output to static HTML file
	output_file("{}/{}_advanced_operons_output_summary.html".format(os.getcwd(), label))
	save(grid)


	# Now, per operon type, plot each individual summary page

	print(' ... Making per operon type page')

	for operon_type in sorted(operons.keys()):
		curr_operon = {operon_type: operons[operon_type]}

		print(' ... ... {}'.format(operon_type))
		
		if len(curr_operon[operon_type]['target_members']) > 1:
			
			all_tabs = []

			div = Div(text="""<b>Detailed view of individual genomic context types:</b></br></br>
					  Below you find multiple <b>tabs</b> coresponding to each individual
					  cluster you see on the scatter plot above.</br>
					  <b>Click</b> on the tab to have a detailed view of the different genomic
					  context clusters, clustered based on the similarity of their family 
					  composition.</br></br>
					  The depiction is interactive, <b>hover</b> and <b>click</b> to get more information!</br></br>  """) 

			# Work on most conserved genomic context figure
			most_common_gc_figure = create_most_common_genomic_features_figure(curr_operon, all_syntenies, families_summary, reference_family = reference_family, family_colors = family_colors, n_flanking5=n_flanking5, n_flanking3=n_flanking3)

			# Make the dendogram
			if tree != None:
				input_targets = [target for operon in operons for target in curr_operon[operon]['target_members']]
			syn_dendogram, syn_den_data = make_dendogram_figure(taxonomy, curr_operon, mode = sort_mode, input_targets = input_targets, show_leafs = False, height_factor = 25*1.2, tree_format = tree_format)

			# Work on the genomic context block
			genomic_context_figure = create_genomic_context_figure(curr_operon, all_syntenies, family_colors, syn_den_data, syn_dendogram, most_common_gc_figure, reference_family, legend_mode = 'species', height_factor = 25*1.2)
			
			# Make the table of family frequencies
			family_table, table_div = create_families_frequency_table(curr_operon, families_summary)

			gc_row = row(syn_dendogram, genomic_context_figure)

			grid = gridplot([[gc_row],[table_div],[family_table]], merge_tools = True)
			curr_tab = Panel(child = grid, title = operon_type)
			all_tabs.append(curr_tab)
		
			all_tabs = Tabs(tabs=all_tabs)
			grid = gridplot([[div], [all_tabs]], merge_tools = True)

			output_file("{}/{}_advanced_operons_interactive_output_{}.html".format(os.getcwd(), label, operon_type))
			save(grid)


#@Timer("write_summary_table", text=lambda secs: f"Time to write summary table: {format_timespan(secs)}")
def write_summary_table(operons, all_syntenies, taxonomy, label = None):

	out_file = '{}_summary_table.tab'.format(label)
	
	with open(out_file, 'w') as outf:
		if 'TM_annotations' in all_syntenies[list(all_syntenies.keys())[0]]['flanking_genes']:
			outf.write('Operon type\tTarget\tAssemblyId\tGene direction\tGene start\tGene end\tRelative gene start\tRelative gene end\tProtein family code\tEntrzID\tProtein name\tTransmembrane/Signal peptide prediction\tSuperkingdom\tPhylum\tClass\tOrder\tGenus\tSpecies\n\n')
		else:
			outf.write('Operon type\tTarget\tAssemblyId\tGene direction\tGene start\tGene end\tRelative gene start\tRelative gene end\tProtein family code\tEntrzID\tProtein name\tSuperkingdom\tPhylum\tClass\tOrder\tGenus\tSpecies\n\n')
		for operon in operons:
			operon_type = operon.split()[-2]
			for target in operons[operon]['target_members']:
				outf.write('\n')
				for superkingdom in taxonomy.keys():
					for phylum in taxonomy[superkingdom].keys():
						for taxclass in taxonomy[superkingdom][phylum].keys():
							for order in taxonomy[superkingdom][phylum][taxclass].keys():
								for genus in taxonomy[superkingdom][phylum][taxclass][order].keys():
									for species in taxonomy[superkingdom][phylum][taxclass][order][genus].keys():
										if target in taxonomy[superkingdom][phylum][taxclass][order][genus][species]['target_members']:
											for i, prot_name in enumerate(all_syntenies[target]['flanking_genes']['names']):
												outf.write('{}\t{}\t{}\t{}\t{}\t{}\t'.format(operon_type, target, all_syntenies[target]['assembly_id'][1], all_syntenies[target]['flanking_genes']['directions'][i], all_syntenies[target]['flanking_genes']['starts'][i], all_syntenies[target]['flanking_genes']['ends'][i]))
												outf.write('{}\t{}\t{}\t{}\t{}\t'.format(all_syntenies[target]['flanking_genes']['relative_starts'][i], all_syntenies[target]['flanking_genes']['relative_ends'][i], all_syntenies[target]['flanking_genes']['families'][i], all_syntenies[target]['flanking_genes']['ncbi_codes'][i], all_syntenies[target]['flanking_genes']['names'][i]))
												
												if 'TM_annotations' in all_syntenies[target]['flanking_genes']:
													outf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(all_syntenies[target]['flanking_genes']['TM_annotations'][i], superkingdom, phylum, taxclass, order, genus, species))
												else:
													outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(superkingdom, phylum, taxclass, order, genus, species))
					
					
# MAIN CODE

@Timer("main", text=lambda secs: f"Time to run main: {format_timespan(secs)}")
def main():

	global timing_log_file  # Declare timing_log_file as global 

	# GET INPUTS
	parser = argparse.ArgumentParser(prog = 'GCsnap v1.0.17', usage = 'GCsnap -targets <targets> -user_email <user_email> [options]', 
									 description = 'GCsnap is a python-based, local tool that generates interactive snapshots\nof conserved protein-coding genomic contexts.',
									 epilog = 'Example: GCsnap -targets PHOL_ECOLI A0A0U4VKN7_9PSED A0A0S1Y445_9BORD -user_email <user_email')

	requiredNamed = parser.add_argument_group('Required arguments')
	optionalNamed = parser.add_argument_group('Options with defaults')

	# required inputs
	requiredNamed.add_argument('-targets', dest='targets', nargs='+', required=True, help='List of input targets. Can be a list of fasta files, a list of text files encompassing a list of protein sequence identifiers, a list of protein sequence identifiers, or a mix of them')
	# optional inputs
	optionalNamed.add_argument('-user_email', dest='user_email', type=str, default = None, help='Email address of the user. May be required to access NCBI databases and is not used for anything else (default: None)')
	optionalNamed.add_argument('-ncbi_api_key', dest='ncbi_api_key', default = None,type=str, help='The key for NCBI API, which allows for up to 10 queries per second to NCBI databases. Shall be obtained after obtaining an NCBI account (default: None)')
	optionalNamed.add_argument('-get_taxonomy', dest='get_taxonomy', default = 'True',type=str, help='Boolean statement to get and map taxonomy information (default: True)')
	optionalNamed.add_argument('-cpu', dest='n_cpu', default = 1,type=int, help='Number of cpus to use (default: 1)')
	optionalNamed.add_argument('-n_flanking', dest='n_flanking', default = 4,type=int, help='Number of flanking sequences (to each side) to take (default: 4)')
	optionalNamed.add_argument('-n_flanking5', dest='n_flanking5', default = 4,type=int, help="Number of flanking sequences to take on the 5' (default: 4)")
	optionalNamed.add_argument('-n_flanking3', dest='n_flanking3', default = 4,type=int, help="Number of flanking sequences to take on the 3' (default: 4)")
	optionalNamed.add_argument('-exclude_partial', dest='exclude_partial', default = False,type=bool, help='Boolean statement to exclude partial operon/genomic_context blocks (default: False)\nIf turned off, partial cases will still be ignored to get the most common genomic features')
	optionalNamed.add_argument('-out_label', dest='out_label', default = 'default',type=str, help='The label to append to the out folder (default: "default"). Important when the input list corresponds to raw sequence identifiers.')
	optionalNamed.add_argument('-out_label_suffix', dest='out_label_suffix', default = '',type=str, help='A suffix to add to the out_label (default: "").')
	optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp',type=str, help='The temporary folder (default: /tmp). May be changed so that intermediary files (e.g., assembly files) are saved somewhere else.')
	optionalNamed.add_argument('-collect_only', dest='collect_only', default = False,type=bool, help='Boolean statement to make GCsnap collect genomic contexts only, without comparing them (default: False).')
	# operon clustering
	optionalNamed.add_argument('-n_max_operons', dest='n_max', default = 30,type=int, help='Maximum number of top most populated operon/genomic_context block types (default: 30)')
	optionalNamed.add_argument('-operon_cluster_advanced', dest='operon_cluster_advanced', default = False,type=bool, help='Boolean statement to use the operon clustering advanced mode (using PacMAP) (default: False)')
	optionalNamed.add_argument('-max_family_freq', dest='max_family_freq', default = 20,type=int, help='Maximum frequency of a family in the set of genomic cotexts found to be considered for advanced operon clustering (default: 20)')
	optionalNamed.add_argument('-min_family_freq', dest='min_family_freq', default = 2,type=int, help='Minimum frequency of a family in the set of genomic cotexts found to be considered for advanced operon clustering (default: 20)')
	#protein family identification
	optionalNamed.add_argument('-n_iterations', dest='num_iterations', default = 1,type=int, help='Number of iterations for all-against-all searches (default: 1). Required to define protein families.')
	optionalNamed.add_argument('-evalue', dest='max_evalue', default = 1e-3,type=float, help='Max e-value at which two sequences are considered to be homologous (default: 1e-3). Required to define protein families.')
	optionalNamed.add_argument('-coverage', dest='min_coverage', default = 70,type=float, help='Minimum coverage of target and subject a match needs to be so that two sequences are considered to be homologous (default: 70). Required to define protein families.')
	optionalNamed.add_argument('-base', dest='default_base', default = 10,type=int, help='Artificial distance value for two sequences that do not match with an E-value better than -evalue (default: 10).')
	optionalNamed.add_argument('-all-against-all_method', dest='clustering_method', default = 'psiblast',type=str, choices=['mmseqs', 'psiblast'], help='Method for clustering (default: psiblast)')
	optionalNamed.add_argument('-psiblast_location', dest='blast', default = 'psiblast',type=str, help='Location of psiBLAST (if not in path) (default: psiblast)')
	optionalNamed.add_argument('-mmseqs_location', dest='mmseqs', default = 'mmseqs',type=str, help='Location of MMseqs (if not in path) (default: mmseqs)')
	# figure making optional inputs
	optionalNamed.add_argument('-genomic_context_cmap', dest='genomic_context_cmap', default = 'Spectral',type=str, help='Color map (as of matplotlib) to assign colors to and plot the syntenic blocks (default: Spectral)')
	optionalNamed.add_argument('-out_format', dest='out_format', default = 'png',type=str, help='Output format of the core figures (default: png)')
	optionalNamed.add_argument('-print_color_summary', dest='print_color_summary', default = False, type=bool, help='Boolean statement to print the RGBA codes of the colors defined (default: False)')
	# annotation optional inputs
	optionalNamed.add_argument('-get_pdb', dest='get_pdbs', default = 'True', type=str, help='Boolean statement to get PDB information for representatives of the families found (default: True)\nTurn off to make it faster.')
	optionalNamed.add_argument('-get_functional_annotations', dest='get_functional_annotations', default = 'True' ,type=str, help='Boolean statement to find functional annotations for representatives of the families found (default: True)\nTurn off to make it faster.')
	optionalNamed.add_argument('-annotate_TM', dest='annotate_TM', default = False, type=bool, help='Boolean statement to find sequence features in the flanking genes (default: False)')
	optionalNamed.add_argument('-annotation_TM_mode', dest='annotation_TM_mode', default = 'uniprot', type=str, choices=['phobius', 'tmhmm', 'uniprot'], help='Method to use to find transmembrane segments (default: uniprot)')
	optionalNamed.add_argument('-annotation_TM_file', dest='annotation_TM_file', default = None, type=str, help='File with pre-computed transmembrane features. Only use when the targets correspond to a single project (no multiple fasta or text files) (default: None)')
	# interactive optional inputs
	optionalNamed.add_argument('-interactive', dest='interactive', default = 'True',type=str, help='Boolean statement to make the interactive html output (default: True). WARNING: It requires the Bokeh python package. It will check if it is installed')
	optionalNamed.add_argument('-gc_legend_mode', dest='gc_legend_mode', default = 'species',type=str, choices=['species', 'ncbi_code'], help='Mode of the genomic context legend (default: species)')
	optionalNamed.add_argument('-min_coocc', dest='min_coocc', default = 0.30,type=float,  help='Minimum maximum co-occurrence of two genes to be connected in the graphs (default: 0.30)')
	optionalNamed.add_argument('-min_freq_accross_contexts', dest='min_family_freq_accross_contexts', default = 30,type=float,  help='Minimum frequency of a family in a conserved genomic context type to be considered as a member (default: 30)')
	optionalNamed.add_argument('-sort_mode', dest='sort_mode', default = 'taxonomy',type=str, choices=['taxonomy', 'as_input', 'tree', 'operon'], help='Mode to sort the genomic contexts (default: taxonomy)')
	optionalNamed.add_argument('-in_tree', dest='in_tree', default = None, type=str, help='Input phylogenetic tree. Only use when the targets correspond to a single project (no multiple fasta or text files) (default: None)')
	optionalNamed.add_argument('-in_tree_format', dest='in_tree_format', default = "newick", type=str, help='Format of the input phylogenetic tree (default: newick)')
	# clans map optional inputs
	optionalNamed.add_argument('-clans_patterns', dest='clans_patterns', default = None,type=str, nargs='+', help='Patterns to identify the clusters to analyse. They will be used to select the individual clusters in the clans map to analyse (default: None).')
	optionalNamed.add_argument('-clans_file', dest='clans_file', default = None,type=str, help='Clans file. Used only for advanced interactive output representation (default: None or input clans file if the input is a clans file and -operon_cluster_advanced is set to True).')


	# print help message if GCsnap is run also without any arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)

	# Define inputs
	args = parser.parse_args()

	targets = args.targets
	Entrez.email = args.user_email

	n_flanking = args.n_flanking
	n_flanking5 = args.n_flanking5
	n_flanking3 = args.n_flanking3

	if n_flanking3 == n_flanking5 and n_flanking != n_flanking3:
		if n_flanking != 4:
			n_flanking3 = n_flanking
			n_flanking5 = n_flanking

	exclude_partial = args.exclude_partial
	collect_only = args.collect_only
	n_cpus = args.n_cpu
	operon_cluster_advanced = args.operon_cluster_advanced
	max_family_freq = args.max_family_freq
	min_family_freq = args.min_family_freq
	min_family_freq_accross_contexts = args.min_family_freq_accross_contexts
	method = args.clustering_method
	n_max = args.n_max
	num_alignments = 50000
	max_evalue = args.max_evalue
	min_coverage = args.min_coverage
	num_iterations = args.num_iterations
	default_base = args.default_base
	out_label = args.out_label
	out_label_suffix = args.out_label_suffix
	genomic_context_cmap = args.genomic_context_cmap
	min_coocc = args.min_coocc
	print_color_summary = args.print_color_summary
	tmp_folder = args.tmp_folder

	gc_legend_mode = args.gc_legend_mode

	annotate_TM = args.annotate_TM
	annotation_TM_mode = args.annotation_TM_mode
	annotation_TM_file = args.annotation_TM_file

	ncbi_api_key = args.ncbi_api_key
	if ncbi_api_key is not None:
		Entrez.api_key = ncbi_api_key

	if annotation_TM_file is not None:
		annotate_TM = True

	get_pdb = args.get_pdbs
	if get_pdb == 'True':
		get_pdb = True
	else:
		get_pdb = False

	get_functional_annotations = args.get_functional_annotations
	if get_functional_annotations == 'True':
		get_functional_annotations = True
	else:
		get_functional_annotations = False

	get_taxonomy = args.get_taxonomy
	if get_taxonomy == 'True':
		get_taxonomy = True
	else:
		get_taxonomy = False

	interactive_output = args.interactive
	if interactive_output == 'True':
		interactive_output = True
	else:
		interactive_output = False

	in_tree = args.in_tree
	in_tree_format = args.in_tree_format
	sort_mode = args.sort_mode
	if in_tree is not None:
		sort_mode = 'tree'

	clans_patterns = args.clans_patterns
	clans_file = args.clans_file
	if len(targets) == 1 and targets[0].endswith('.clans') and clans_file is None:
		clans_file = targets[0]

	if clans_file is not None:
		clans_file = os.path.abspath(clans_file)

	# define programs location
	blast = args.blast
	mmseqs = args.mmseqs

	# install cache to make it faster
	if sys.version_info[0] == 3 and sys.version_info[1] < 8:
		requests_cache.install_cache()

	# set starting directory
	starting_directory = os.getcwd()
	print(os.getcwd())

	## START PIPELINE HERE

	try:

		# Download and parse RefSeq and Genbank databases
		print('\nDownloading and parsing RefSeq and Genbank summary tables\n')
		refseq_gb_assembly_map = download_and_parse_refseq_and_gb_databases()

		t_all_steps = Timer('All steps 0-10', text=lambda secs: f"Time to run steps 1-9: {format_timespan(secs)}")
		t_all_steps.start()

		t0 = Timer('Step 0: Parse Targets', text=lambda secs: f"Time to run step 1: {format_timespan(secs)}")
		t0.start()
		# Get the list of target codes
		print('\nParsing targets\n')
		targets = parse_targets(targets, label = out_label, clans_patterns = clans_patterns)
		t0.stop()

		for out_label in targets:

			start = time.time()
			# Define job and create working directory
			print("\nWorking on job '{}'".format(out_label))
			if len(out_label_suffix) > 0:
				working_dir = '{}/{}_{}'.format(starting_directory, out_label, out_label_suffix)
			else:
				working_dir = '{}/{}'.format(starting_directory, out_label)

			if not os.path.isdir(working_dir):
				os.mkdir(working_dir)
			os.chdir(working_dir)

			# Write arguments file, which summarises all arguments used
			write_arguments_file(args, out_label)

			# Collect the genomic_context of all target ncbi codes
			curr_targets = targets[out_label]

			if len(curr_targets) > 1:
				print("\n 1. Collecting the genomic contexts of {} unique input entrezIDs (may take some time)\n".format(len(curr_targets)))
				t1 = Timer('Step 1: Collecting the genomic contexts', text=lambda secs: f"Time to run step 1: {format_timespan(secs)}")
				t1.start()
				all_syntenies = get_genomic_context_information_for_ncbi_codes(curr_targets, refseq_gb_assembly_map = refseq_gb_assembly_map, n_flanking5 = n_flanking5, n_flanking3 = n_flanking3, exclude_partial = exclude_partial, tmp_folder = tmp_folder, threads = n_cpus)
				t1.stop()

				if len(all_syntenies) > 1:
					if not collect_only:
						# Find shared protein families by running all-against-all blast searches for all proteins collected
						print("\n 2. Finding protein families (may take some time depending on the number of flanking sequences taken)\n")
						t2 = Timer('Step 2: Finding protein families', text=lambda secs: f"Time run step 2: {format_timespan(secs)}")
						t2.start()
						all_syntenies, protein_families_summary = find_and_add_protein_families(all_syntenies, out_label = out_label, num_threads = n_cpus, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, blast = blast, mmseqs = mmseqs, min_coverage = min_coverage, default_base = default_base, tmp_folder = tmp_folder, method = method)
						t2.stop()

						# Search for functional information and pdb structures (experimental or homology-models, in Swiss-repository) for representatives of the families found
						if get_pdb or get_functional_annotations:

							print("\n 3. Annotating functions and/or finding structures for the protein families found\n")
							t3 = Timer('Step 3: Annotating functions and structures', text=lambda secs: f"Time run step 3: {format_timespan(secs)}")
							t3.start()
							protein_families_summary = update_families_with_functions_and_structures(protein_families_summary, get_pdb = get_pdb, get_functional_annotations = get_functional_annotations, threads = n_cpus)
							t3.stop()
						else:
							print("\n 3. Neither functions will be annotated, neither structures will be searched\n")

						# Find operon/genomic_context types by clustering them by similarity
						print("\n 4. Finding operon/genomic_context types\n")
						t4 = Timer('Step 4: Finding operon/genomic_context', text=lambda secs: f"Time run step 4: {format_timespan(secs)}")
						t4.start()
						all_syntenies, operon_types_summary = find_and_add_operon_types(all_syntenies, protein_families_summary, label = out_label, advanced = operon_cluster_advanced, min_family_freq = min_family_freq, max_family_freq = max_family_freq)
						t4.stop()

						# Select top N most populated operon/genomic_context types
						print("\n 5. Selecting top {} most common operon/genomic_context types\n".format(n_max))
						t5 = Timer('Step 5: Selecting most common operon/genomic context', text=lambda secs: f"Time run step 5: {format_timespan(secs)}")
						t5.start()
						selected_operons, most_populated_operon = find_most_populated_operon_types(operon_types_summary, nmax = n_max)
						json.dump(selected_operons, open('selected_operons.json', 'w'), indent = 4)
						t5.stop()

						# get taxonomy information
						if get_taxonomy:
							# Map taxonomy to the input targets. Load if already precomputed
							print("\n 6. Mapping taxonomy (may take some time)\n")
							t6 = Timer('Step 6: Mapping taxonomy', text=lambda secs: f"Time run step 6: {format_timespan(secs)}")
							t6.start()
							taxonomy = map_taxonomy_to_targets(all_syntenies, threads = n_cpus)
							json.dump(taxonomy, open('taxonomy.json', 'w'), indent = 4)
							t6.stop()

						else:
							print("\n 6. Taxonomy will not be mapped\n")
							t6b = Timer('Step 6b: No mapping', text=lambda secs: f"Time run step 6b: {format_timespan(secs)}")
							t6b.start()
							taxonomy = map_taxonomy_to_targets(all_syntenies, mode = 'as_input')
							t6b.stop()

						# Annotate transmembrane segments
						if annotate_TM:
							# Try to annotate transmembrane segments
							print("\n 7. Finding ALL proteins with transmembrane segments and signal peptides\n")
							t7 = Timer('Step 7: Finding ALL proteins with transmembrane segments', text=lambda secs: f"Time run step 7 using {annotation_TM_mode}: {format_timespan(secs)}")
							t7.start()
							all_syntenies = annotate_TMs_in_all(all_syntenies, annotation_TM_mode, annotation_TM_file, label = out_label)
							t7.stop()

						else:
							print("\n 7. Transmembrane segments and signal peptides will not be searched\n")

						# Make operon/genomic_context conservation figure
						print("\n 8. Making operon/genomic_context blocks figure\n")
						try:
							t8 = Timer('Step 8: Making genomic context figure', text=lambda secs: f"Time run step 8: {format_timespan(secs)}")
							t8.start()
							make_genomic_context_figure(selected_operons, most_populated_operon, all_syntenies, protein_families_summary, cmap = genomic_context_cmap, label = out_label, out_format = args.out_format)
							t8.stop()
						except:
							print(' ... Images not created due to minor errors (likely they are too big)')
						# Make interactive HTML output
						if interactive_output:
							print("\n 9. Making interactive html output file\n")
							#check if Bokeh is available
							# try:
							t9 = Timer('Step 9: Making interactive html output', text=lambda secs: f"Time run step 9 (operon cluster advanced mode: {operon_cluster_advanced}): {format_timespan(secs)}")
							t9.start()
							if operon_cluster_advanced:
								make_advanced_interactive_genomic_context_figures(selected_operons, all_syntenies, protein_families_summary, taxonomy, most_populated_operon, input_targets = curr_targets, tree = in_tree, gc_legend_mode = gc_legend_mode, cmap = genomic_context_cmap, label = out_label, sort_mode = sort_mode, min_coocc = min_coocc, n_flanking5=n_flanking5, n_flanking3=n_flanking3, tree_format = in_tree_format, max_family_freq=max_family_freq, min_family_freq=min_family_freq, min_family_freq_accross_contexts=min_family_freq_accross_contexts, clans_file=clans_file, out_label = out_label, num_threads = n_cpus, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, blast = blast, mmseqs = mmseqs, min_coverage = min_coverage, default_base = default_base, tmp_folder = tmp_folder, method = method)
							else:
								make_interactive_genomic_context_figure(selected_operons, all_syntenies, protein_families_summary, taxonomy, most_populated_operon, input_targets = curr_targets, tree = in_tree, gc_legend_mode = gc_legend_mode, cmap = genomic_context_cmap, label = out_label, sort_mode = sort_mode, min_coocc = min_coocc, n_flanking5=n_flanking5, n_flanking3=n_flanking3, tree_format = in_tree_format)
							t9.stop()

							# except:
							#	 print(sys.exc_info())
							#	 print(' --> ERROR: Not possible to generate the interactive output. If the error is about Bokeh, check if the correct version is installed and run again. You can install it with "pip install bokeh==1.3.4"')
							#	 print("\n ... Interactive html output file will not be generated\n")

						else:
							print("\n 9. Interactive html output file will not be generated\n")

						# Write summary table
						t10 = Timer('Step 10: Write output', text=lambda secs: f"Time run step 9 (operon cluster advanced mode: {operon_cluster_advanced}): {format_timespan(secs)}")
						t10.start()
						print("\n Finished {}: Writting summary table\n".format(out_label))
						write_summary_table(selected_operons, all_syntenies, taxonomy, label = out_label)
						json.dump(protein_families_summary, open('protein_families_summary.json', 'w'), indent = 4)

					else:
						print(' GCsnap was asked to collect genomic context only. Will not proceed further.')
						t10 = Timer('Step 10: Write output', text=lambda secs: f"Time run step 9 (operon cluster advanced mode: {operon_cluster_advanced}): {format_timespan(secs)}")
						t10.start()

					json.dump(all_syntenies, open('all_syntenies.json', 'w'), indent = 4)
					t10.stop()

					end = time.time()
					numb_seconds = end - start
					print("\n#### Finished {} after: {} \n".format(out_label, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

				else:
					print('\n --> WARNING: Found genomic contexts for less than 2 targets. Job {} will not be further analysed.'.format(out_label))
			else:
				print('\n --> WARNING: Job {} has less than 2 targets. It will not be analysed.'.format(out_label))

			t_all_steps.stop()

		with open(timing_log_file, 'a+') as f:
			f.write(f'{out_label}\n')
			f.write(f'nr_of_targets: {len(targets[list(targets.keys())[0]])}\n')
			f.write(f'\tTotal_time\tCount\tMax\tMin\tMean\tStd\tMedian\n')

			for timing in Timer.timers:
				line_to_write = [timing, '\t', str(Timer.timers.total(timing)), '\t', str(Timer.timers.count(timing)),
								 '\t', str(Timer.timers.max(timing)), '\t', str(Timer.timers.min(timing)),
								 '\t', str(Timer.timers.mean(timing)), '\t', str(Timer.timers.stdev(timing)),
								 '\t', str(Timer.timers.median(timing)), '\n']
				f.writelines(line_to_write)

	except Exception as e:
		print(e)
		with open(timing_log_file, 'a+') as f:
			f.write(f'{out_label}\n')
			f.write(f'nr_of_targets: {len(targets[list(targets.keys())[0]])}\n')
			f.write(f'\tTotal_time\tCount\tMax\tMin\tMean\tStd\tMedian\n')

			for timing in Timer.timers:
				line_to_write = [timing, '\t', str(Timer.timers.total(timing)), '\t', str(Timer.timers.count(timing)),
								 '\t', str(Timer.timers.max(timing)), '\t', str(Timer.timers.min(timing)),
								 '\t', str(Timer.timers.mean(timing)), '\t', str(Timer.timers.stdev(timing)),
								 '\t', str(Timer.timers.median(timing)), '\n']
				f.writelines(line_to_write)


if __name__ == '__main__':
	# now = datetime.now()
	# current_time = now.strftime("%d-%m-%Y_%H-%M-%S")
	# timing_log_file = f'times_{current_time}.tsv'

	main_start = time.time()
	main()
	main_end = time.time()
	main_numb_seconds = main_end - main_start

	# add time for running main() to timing log file
	with open(timing_log_file, 'a+') as f:
		line_to_write = ['main', '\t', str(Timer.timers.total('main')), '\t', str(Timer.timers.count('main')),
						 '\t', str(Timer.timers.max('main')), '\t', str(Timer.timers.min('main')),
						 '\t', str(Timer.timers.mean('main')), '\t', str(Timer.timers.stdev('main')),
						 '\t', str(Timer.timers.median('main')), '\n']
		f.writelines(line_to_write)
