#!/usr/bin/env python3

import configparser
import argparse
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d
import numpy as np
import pandas as pd
import sys
import time
import gc
import itertools
import psutil #wrapper
import os #wrapper
from functools import partial
from joblib import Parallel, delayed

def get_config():
	config=configparser.ConfigParser()
	config.read("sim_config.ini")
	return config["I"] #INSERTION or ROI

def get_arguments():
	parser = argparse.ArgumentParser(description='Simulation of insertions (I) or regions-of-interests (ROI).')
	parser.add_argument('--reference_genome_path', type=str, help='Path to reference genome')
	parser.add_argument('--vector_sequence_path', type=str, help='Path to vector sequence/ not really used currently')
	parser.add_argument('--sequenced_data_path', type=str, help='Path to sequenced data for custom read length distribution')
	parser.add_argument('--output_path', type=str, help='Output path')
	parser.add_argument('--experiment_name', type=str, help='Experiment name')
	parser.add_argument('--output_path_plots', type=str, help='Output path for plots')
	parser.add_argument('--insertion_probability', type=float, help='Insertion probability: [0-1]')
	parser.add_argument('--chr_restriction', type=str, help='Chromosome restriction. If "unrestricted", all chromosomes from the reference are used and "_" are changed to "-" in the gene names".')
	parser.add_argument('--bedpath', type=str, help='Path to BED file for insertions. Insertions are randomly spread into the respective intervals. Interval-size is considered. Special: If num insertions and num of bed entries match, each entry will be chosen once, disregaridng its length. This allows for fixed insertions!')
	parser.add_argument('--barcode_weights', type=str, help='Barcode weights. E.g. {"Barcode_0: 5, "Barcode_1: 2"}: A read will be be five times more likely to be from Barcode 0 than from Barcode 2,3,...n. ')
	parser.add_argument('--chromosome_weights', type=str, help='Chromosome weights. E.g. {"chr2":0}: There will be no reads drawn from chromosome 2. Resembles monosomy option with 0.')
	parser.add_argument('--insertion_numbers', type=int, help='Number of insertions')
	parser.add_argument('--n_barcodes', type=int, help='Number of barcodes, i.e. number of genomes')
	parser.add_argument('--iterations', type=int, help='Number of iterations')
	parser.add_argument('--scaling', type=bool, help='defines whether there is genome scaling necessary. Currently, onle True/False are implemented. If true, 50%% of the findings will be discarded to account for the difference of haplotype (reference genome) and the physiologicla state (diploid)')
	parser.add_argument('--parallel_jobs', type=int, help='Number of parallel jobs')
	parser.add_argument('--mode', type=str, help='Mode. Either I for Insertion moder or ROI for region of interest mode. ROI needs a BED file with ROI genomic locations based on the referece genome.')
	parser.add_argument('--coverages', type=str, help='Coverages. N times the reference genome.')
	parser.add_argument('--mean_read_lengths', type=str, help='Mean read lengths. Only if no custom read length distirbution is used.')
	parser.add_argument('--roi_bedpath', type=str, help='Path to ROI BED file')
	parser.add_argument('--blocked_regions_bedpath', type=str, help='Path to blocked regions BED file')
	parser.add_argument('--barcodes_to_check_blocked_regions', type=str, help='Barcodes to check blocked regions. Only these Barcodes will be affected by the blocking.')
	parser.add_argument('--monosomie', type=str, help='Monosomy. Same as chromosome weigths set to 0.')
	
	return parser.parse_args()


def parse_comma_separated_list(value, value_type):
	if value is None or value.lower() == 'none':
		return None
	return [value_type(item) for item in value.split(',')]

def parse_stringified_dict(value):
	if value is None or value.lower() == 'none':
		return None
	try:
		return eval(value)
	except (SyntaxError, NameError):
		raise ValueError("Invalid dictionary format for argument")

def main():
	config = get_config()
	args = get_arguments()

	reference_genome_path = args.reference_genome_path or config.get('reference_genome_path')
	vector_sequence_path = args.vector_sequence_path or config.get('vector_sequence_path')
	sequenced_data_path = args.sequenced_data_path or config.get('sequenced_data_path')
	output_path = args.output_path or config.get('output_path')
	experiment_name = args.experiment_name or config.get('experiment_name')
	output_path_plots = args.output_path_plots or config.get('output_path_plots')
	insertion_probability = args.insertion_probability or float(config.get('insertion_probability'))
	chr_restriction = args.chr_restriction or config.get('chr_restriction')
	bedpath = args.bedpath or config.get('bedpath')
	barcode_weights = parse_stringified_dict(args.barcode_weights) if args.barcode_weights else parse_stringified_dict(config.get('barcode_weights'))
	chromosome_weights = parse_stringified_dict(args.chromosome_weights) if args.chromosome_weights else parse_stringified_dict(config.get('chromosome_weights'))
	insertion_numbers = args.insertion_numbers or int(config.get('insertion_numbers'))
	n_barcodes = args.n_barcodes or int(config.get('n_barcodes'))
	iterations = args.iterations or int(config.get('iterations'))
	scaling = args.scaling if args.scaling is not None else config.getboolean('scaling')
	parallel_jobs = args.parallel_jobs or int(config.get('parallel_jobs'))
	mode = args.mode or config.get('mode')
	coverages = parse_comma_separated_list(args.coverages, float) if args.coverages else parse_comma_separated_list(config.get('coverages'), float)
	mean_read_lengths = parse_comma_separated_list(args.mean_read_lengths, int) if args.mean_read_lengths else parse_comma_separated_list(config.get('mean_read_lengths'), int)
	roi_bedpath = args.roi_bedpath or config.get('roi_bedpath')
	blocked_regions_bedpath = args.blocked_regions_bedpath or config.get('blocked_regions_bedpath')
	barcodes_to_check_blocked_regions = parse_comma_separated_list(args.barcodes_to_check_blocked_regions, str) if args.barcodes_to_check_blocked_regions else parse_comma_separated_list(config.get('barcodes_to_check_blocked_regions'), str)
	monosomie = parse_comma_separated_list(args.monosomie, str) if args.monosomie else parse_comma_separated_list(config.get('monosomie'), str)

	
	#combinations
	combinations = itertools.product(mean_read_lengths, coverages)

	### WRAPPER START
	'''
	Wrapper to monitor the resources used by key steps.
	'''
	def elapsed_since(start):
		return time.strftime("%H:%M:%S", time.gmtime(time.time() - start))


	def get_process_memory():
		process = psutil.Process(os.getpid())
		return process.memory_info().rss


	def profile(func):
		def wrapper(*args, **kwargs):
			mem_before = get_process_memory()
			start = time.time()
			result = func(*args, **kwargs)
			elapsed_time = elapsed_since(start)
			mem_after = get_process_memory()
			print("{}: memory before: {:,}, after: {:,}, consumed: {:,}; exec time: {}".format(
				func.__name__,
				mem_before, mem_after, mem_after - mem_before,
				elapsed_time))
			return result
		return wrapper
	#### WRAPPER END

	def check_barcoding(n_barcodes):
		'''
		Evaluates barcoding condition.
		'''
		if n_barcodes >= 1:
			return True
		return False

	def readbed(bedpath, list_of_chromosomes_in_reference, barcoding=False):
		'''
		Reads the bed file and only keeps entries from chromosomes in the reference.
		If the weights column is missing, it adds a weights column with 0 for each row.
		Barcoding only becomes true for the bed-guided insertion placement.
		'''
		try:
			# Attempt to read the bed file with the weight column
			bed = pd.read_csv(bedpath, sep='\t', header=None, usecols=[0, 1, 2, 3, 4], names=['chrom', 'start', 'end', 'ID', 'weight'])
		except ValueError as e:
			# If the weight column is missing, read without it and add a weight column with 0
			if "Too many columns specified: expected 5 and found 4" in str(e):
				print("Weight column not found but ID exists. Adding a default weight of 0.")
				bed = pd.read_csv(bedpath, sep='\t', header=None, usecols=[0, 1, 2, 3],names=['chrom', 'start', 'end', 'ID'])
				bed['weight'] = 0  # Add the weight column with 0 as default
			else:
				raise e  # Raise other exceptions
		except:
			try: 
				print("Try reading BED minimally...")
				bed = pd.read_csv(bedpath, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
				bed['ID'] = ''  # Add a placeholder ID if missing
				bed['weight'] = 0  # Add the weight column with 0 as default
			except Exception as e:
				print("Could not read the BED or no BED provided:", e)
				return None 

		# Handling barcoding if required
		if barcoding:
			print('Barcoding selected: Transforming the names in the bed... Adding the following barcode...')
			print('_'.join(list(list_of_chromosomes_in_reference)[0].split("_")[:-1]))
			barcode = '_'.join(list(list_of_chromosomes_in_reference)[0].split("_")[:-1])
			bed['chrom'] = barcode + '_' + bed['chrom'].astype(str)
			bed = bed[bed["chrom"].isin(list_of_chromosomes_in_reference)]
		else:
			print('Only keeping the following chromosomes...')
			print(list_of_chromosomes_in_reference)
			bed = bed[bed["chrom"].isin(list_of_chromosomes_in_reference)]
		
		return bed

	def update_coordinates(df, input_chromosome_dict, include_weights=False):
		'''
		Uses a df (bed-like) and a dict of chromosome broders in global format and returns the global cooridnates of the bed entries.
		'''
		try:
			updated_coordinates = {}
			for index, row in df.iterrows():
				ID=row["ID"]
				chrom = row['chrom']
				start = row['start'] + input_chromosome_dict[chrom][0]
				end = row['end'] + input_chromosome_dict[chrom][0]
				if include_weights: #adds weights if there are any /// a weight of 1 means that the corresponding coordinates will be blocked with a probability of 1 
					updated_coordinates[ID] = {row['weight'], start, end}
				else:
					updated_coordinates[ID] = {start, end}
			return updated_coordinates
		except:
			return None

	def convert_to_normal_bed(chromosome_dict, items_dict):
		'''
		Converts items dictionary to DataFrame in BED format using chromosome dictionary.
		'''
		bed_data = []
		for item, coordinates in items_dict.items():
			global_start, global_end = coordinates
			for chrom, dimensions in chromosome_dict.items():
				chr_start, chr_end = dimensions
				if chr_start < global_start < chr_end:
					bed_data.append([chrom, global_start - chr_start, global_end - chr_start, item])

		bed_df = pd.DataFrame(bed_data, columns=['chrom', 'start', 'end', 'item'])
		return bed_df

	def add_monosomy_regions(monosomie, input_chromosome_dict, masked_regions=None, chromosome_weights=None):
		'''
		Adds monosomy regions to the masked regions with optional weights.
		'''
		if monosomie is None:
			return masked_regions

		if not masked_regions:
			masked_regions = {}
		try:
			for chromosome in monosomie:
				masked_regions[chromosome] = input_chromosome_dict[chromosome].copy()
				if chromosome_weights and chromosome in chromosome_weights:
					masked_regions[chromosome].insert(0, chromosome_weights[chromosome])
				else:
					print("No weight provided for %s. Added 0 probability." %chromosome)
					masked_regions[chromosome].insert(0, 0)
		except KeyError:
			print("Chromosome names in monosomie list and reference genome do not overlap")

		return masked_regions


	def collapse_fasta(path_to_fasta):
		'''
		Collapses fasta file into single fasta string
		'''
		with open(path_to_fasta, 'r') as fasta_file:
			seqList=[]
			for record in SeqIO.parse(fasta_file, 'fasta'): 
				seqList.append(record)
		return ''.join(seqList)

	def pseudo_fasta_coordinates(path_to_fasta):
		'''
		Takes in FASTA with separate chr entries and returns a table with the pseudo entry borders of each entry after the collapsing of the entries. Currently filters atypical chromosomes.
		'''
		with open(path_to_fasta, 'r') as fasta_file:
			entries={}
			seqList=[]
			updated_length=0
			for record in SeqIO.parse(fasta_file, 'fasta'):
				if chr_restriction == "unrestricted":
					#this part makes sure now downstream errors occur while doing splits by _
					record.id = record.id.replace("_", "-")
					entries[record.id] = [updated_length, updated_length + len(record.seq)]
					updated_length = updated_length + len(record.seq)
					seqList.append(str(record.seq))
				else:   
					if '_' not in record.id and 'M' not in record.id:
							entries[record.id] = [updated_length, updated_length + len(record.seq)]
							updated_length = updated_length + len(record.seq)
							seqList.append(str(record.seq))

		return ''.join(seqList), entries

	def barcode_genome(chromosome_dir, barcode):
		'''
		Adds prefix to chromosome dict
		'''
		chromosome_dir = {f'Barcode_{barcode}_{k}': v for k, v in chromosome_dir.items()}
		return chromosome_dir

	def roi_barcoding(roi_dict, n_barcodes):
		'''
		Add _i suffix to each key in the dictionary for each index i in the specified range.
		'''
		new_dict = {}
		for key, value in roi_dict.items():
			for i in range(n_barcodes):
				new_key = f"{key}_{i}"
				new_dict[new_key] = value
		return new_dict


	def get_chromosome(insert_position, chromosome_dir):
		'''
		Checks within which chromosome the random insertion landed.
		'''
		for chromosome, coordinates in chromosome_dir.items():
			start, end = coordinates
			if start <= insert_position <= end:
				return chromosome
		return None

	def insertions_per_chromosome(chromosome_dir, insertion_dir):
		'''
		Creates a table that lists the number of insertions per chromosome.
		'''
		
		# Convert chromosome_dir and insertion_dir to sets for faster lookups
		chromosome_set = {chromosome: tuple(coordinates) for chromosome, coordinates in chromosome_dir.items()}
		insertion_set = {insertion: tuple(coordinates) for insertion, coordinates in insertion_dir.items()}
		chromosome_lengths = {}
		
		# Iterate through chromosome_dir to get chromosome lengths
		for chromosome, coordinates in chromosome_set.items():
			chromosome_lengths[chromosome] = coordinates[1] - coordinates[0]

		# Initialize a dictionary to count insertions per chromosome
		insertions_per_chromosome = {chromosome: 0 for chromosome in chromosome_dir}

		# Iterate through insertions
		for insertion, insertion_coordinates in insertion_set.items():
			# Check if the insertion falls within any chromosome
			for chromosome, chromosome_coordinates in chromosome_set.items():
				if chromosome_coordinates[0] <= insertion_coordinates[0] < chromosome_coordinates[1] or \
				   chromosome_coordinates[0] < insertion_coordinates[1] <= chromosome_coordinates[1]:
					insertions_per_chromosome[chromosome] += 1

		# Convert the dictionary to a DataFrame
		insertions_df = pd.DataFrame.from_dict(insertions_per_chromosome, orient='index', columns=['Num_Insertions'])
		insertions_df.index.name = 'Chromosome'
		# Add a column for chromosome lengths
		insertions_df['Chromosome_Length'] = insertions_df.index.map(chromosome_lengths)

		return insertions_df

	@profile
	def add_insertions_to_genome_sequence_with_bed(reference_sequence, insertion_sequence, num_insertions, chromosome_dir, bed_df=None): #costs 3 billion per barcode!
		'''
		Randomly add insertion sequence into the reference genome or within specified regions.

		In principle, for each insertion, a random position on the reference sequence is chosen and its coordinates are stored. 
		For all follwoing insertions, the coordinates of their location are updated if needed (i.e. if they happend to insert at an earlier position)
		In case of a bed-guided insertion, each region in the file is assigned a probability according to its length.
		'''
		position = {}
		if bed_df is not None:
			print("BED guided insertion pattern...")
			# Step 1: Calculate probabilities based on region lengths
			print("Calculating insertion probabilities (region length / sum of all regions lengths)...")
			region_lengths = bed_df['end'] - bed_df['start']
			region_probabilities = region_lengths / region_lengths.sum()
			updated_reference_sequence = reference_sequence

			for i in range(num_insertions):
				if random.random() < insertion_probability: 
					# Step 2: Randomly select insertion regions #so that each region is selected once!
					if len(bed_df.index) == num_insertions:
						selected_region = bed_df.iloc[i]
					else:
						selected_region_index = np.random.choice(bed_df.index, p=region_probabilities)
						selected_region = bed_df.iloc[selected_region_index]

					# Step 3: Perform insertions within selected regions
					chromosome = selected_region['chrom']
					chromosome_range = chromosome_dir[chromosome]

					insert_position = random.randint(selected_region['start'], selected_region['end'])
					# Adjust insertion position to the global genomic coordinates
					global_insert_position = chromosome_range[0] + insert_position

					# Insert the insertion sequence into the reference sequence
					updated_reference_sequence = (
						updated_reference_sequence[:global_insert_position] +
						insertion_sequence +
						updated_reference_sequence[global_insert_position:]
					)

					# Update position table
					for key, value in position.items():
					# If the insertion is after the current position, update the position
						if value[0] >= global_insert_position:
							position[key][0] += len(insertion_sequence)
							position[key][1] += len(insertion_sequence)

					# Add the new insertion position and (barcoded) name
					insertion_name = chromosome.split('chr')[0] + "insertion_%s" %i
					position[insertion_name] = [global_insert_position, global_insert_position + len(insertion_sequence)]


			return updated_reference_sequence, position
		
		#if no bed is provided
		for i in range(num_insertions):
			if random.random() < insertion_probability:
				# Choose a random position to insert the smaller sequence
				insert_position = random.randint(0, len(reference_sequence))
				
				#check in which chr it landed
				chromosome = get_chromosome(insert_position, chromosome_dir)
				
				# Insert the insertion sequence at the chosen position
				updated_reference_sequence = (
					reference_sequence[:insert_position] +
					insertion_sequence +
					reference_sequence[insert_position:]
				)
				# Update position table
				for key, value in position.items():
				# If the insertion is after the current position, update the position
					if value[0] >= insert_position:
						position[key][0] += len(insertion_sequence)
						position[key][1] += len(insertion_sequence)

				# Add the new insertion position and (barcoded) name
				insertion_name = chromosome.split('chr')[0] + "insertion_%s" %i
				position[insertion_name]= [insert_position, insert_position + len(insertion_sequence)]

		return updated_reference_sequence, position

	@profile
	def create_barcoded_insertion_genome(reference_genome_path, bedpath, insertion_fasta, n_barcodes):
		'''
		Pre-processment step of the insertion mode.

		1.) The reference genome cooridnates are transformed into a string-like format (One single FASTA string) and the chromsome borders are stored
		2.) Based on user input, masked regions are defined and weighted
		3.) For each barcode, the chromosome borders get a prefix and inserttion positions are randomly chosen. If an insertion-bed is defined, the insertions are only placed within these regions accoridng to their length-based weights.
		4.) The length of the reference genome (for coverage estimation), the insertion cooridnates, the masked regions, and the chromosome border dict are returned
		''' 
		print("Create Barcoded Genome...")
		insertion_dict={}
		#1
		fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
		
		#2 #create blocked regions file
		blocked_bed = readbed(blocked_regions_bedpath, chromosome_dir.keys()) #barcoding should be default false
		masked_regions = update_coordinates(blocked_bed, chromosome_dir, include_weights=True)
		masked_regions = add_monosomy_regions(monosomie, chromosome_dir, masked_regions=masked_regions, chromosome_weights=chromosome_weights)
		#3
		for i in range(n_barcodes):
			barcoded_chromosome_dir = barcode_genome(chromosome_dir, i)
			#optional bed-guided insertion
			bed_df = readbed(bedpath, barcoded_chromosome_dir.keys(), barcoding=check_barcoding(n_barcodes))
			mod_fasta, insertion_dir = add_insertions_to_genome_sequence_with_bed(fasta, insertion_fasta, insertion_numbers, barcoded_chromosome_dir, bed_df)
			insertion_dict.update(insertion_dir)
			del barcoded_chromosome_dir, bed_df, insertion_dir
		#4
		return len(mod_fasta), insertion_dict, masked_regions, chromosome_dir

	def get_read_length_distribution_from_real_data(path_to_sequenced_fasta):
		'''
		Uses an input fasta and generates a custom read length distribution.
		'''
		lengths=[]
		for record in SeqIO.parse(path_to_sequenced_fasta, 'fasta'):
				lengths.append(len(record.seq))
		return lengths

	def generate_reads(fasta, read_length_distribution, num_reads):
		'''
		Generate reads based on custom read length distribution and pulled from fasta ref.
		'''
		reads = []
		for read_length in read_length_distribution:
			start_position = random.randint(0, len(fasta) - read_length)
			read = fasta[start_position:start_position + read_length]
			reads.append(read)
		return reads

	def check_if_blocked_region(start_position, read_length, blocked_region_coordinates=None):
		'''
		Checks whether the start position is inside the coordinates of a masked/blocked region.
		Blocked can mean that no reads are obtained from there or only with a certain probability.
		This is useful to simulate chromosomal aberrations or targeted sequencing
		Probability of 0 means all reads with a start or end position in the respective interval are discarded! (0,5: Half of them are discarded, 1: None are discarded)
		'''
		if blocked_region_coordinates:
			for values in blocked_region_coordinates.values():
				probability, start, stop = values
				if (start < start_position < stop) or (start < start_position + read_length < stop): #checks if start lays in the blocked region or if the end lays in a blocked region
					# The start position falls within a blocked region
					if random.random() >= probability:
						return True  # position is blocked 

	@profile
	def generate_reads_based_on_coverage(fasta, read_length_distribution, coverage, PRECOMPUTE_RANDOM):
		'''
		Randomly pulls a read of size X derived from the read length distribution from the fasta until the fasta is N times covered (coverage).
		
		If barcoding:
		Each pulled read will be assigned to a barcode by the barcodes weight. This accounts for the fact that we have a heterogenous mixture of cells for sequencing.

		If masked regions are provided:
		Each pulled read will be checked for its location. If its position lays within a masked region, it is discarded accoridng to the defined probability.
		(Only specific barcodes can have masked regions as well. This is implemented via "barcodes_to_check_blocked_regions" (default: All barcodes are used))

		'''
		print("Coverage: " + str(coverage))
		print("Pulling reads...")
		covered_length = 0
		total_length = fasta
		#reads = [] #only a goodf idea if we are not testing high coverages, otherwise memory is floated
		read_coordinates = {}
		barcode_names = ["Barcode_" + str(i) for i in range(n_barcodes)]

		#Reads pulled until desired coverage is reached
		while covered_length < coverage * total_length:
			random_barcode = random.choices(barcode_names, [get_weighted_probabilities(i, n_barcodes, barcode_weights) for i in barcode_names])[0] #chooses one of the barcodes based on weighted probability
			read_length = PRECOMPUTE_RANDOM.pop()
			start_position = random.randint(0, total_length - read_length)
			# Check if the random barcode is in the list of barcodes that require checking for blocked regions
			random_barcode_number = random_barcode.split('_')[-1] #otherwise user input needs to be weird

			if barcodes_to_check_blocked_regions and random_barcode_number in barcodes_to_check_blocked_regions:
				# Check if the start position falls within a blocked region
				if check_if_blocked_region(start_position, read_length, masked_regions):
					# If the start position is in a blocked region, continue to the next iteration
					continue
			
			# Record the read coordinates
			read_coordinates[f"Read_{len(read_coordinates)}_{random_barcode}"] = (start_position, start_position + read_length)
			# Update the total covered length
			covered_length += read_length
		return read_coordinates, covered_length



	def generate_read_length_distribution(num_reads, mean_read_length, distribution='lognormal', **kwargs):
		'''
		Generate a list of read lengths with a customizable mean read length and distribution.
		'''
		if distribution == 'normal':
			# Generate read lengths from a normal distribution
			read_lengths = np.random.normal(mean_read_length, mean_read_length / 10, num_reads)
		elif distribution == 'lognormal':
			# Generate read lengths from a log-normal distribution
			#adjust the mean so it matches the lognormal case
			read_lengths = np.random.lognormal(mean=np.log(mean_read_length) - 0.5, sigma=1.0, size=num_reads)
		else:
			raise ValueError("Unsupported distribution. Supported options: 'normal', 'lognormal'.")
		
		# Filter out zero-length reads
		read_lengths = np.round(read_lengths).astype(int)
		read_lengths = read_lengths[read_lengths > 0]
		
		return read_lengths #not returning list


	def get_weighted_probabilities(insertion_name,n_barcodes, weights_dict):
		'''
		Uses a target name and checks for a key in the weights  and returns the weighting factor.
		'''
		if weights_dict is not None:
			# Calculate the common denominator
			common_denominator = (sum(weights_dict.values()) + n_barcodes - len(weights_dict)) * n_barcodes
			#Weight provided: Weighted share of the barcode of the common denominator used
			for key in weights_dict: 
				if any(key == part for part in insertion_name.split("_")) or key == insertion_name: #added the part after the or
					#print((weights_dict[key] * n_barcodes) / common_denominator)
					return (weights_dict[key] * n_barcodes) / common_denominator
			
			#No weight provided for this barcode, using the barcodes share of the common denominator 
			return n_barcodes / common_denominator
		else:
			#print("No weights provided, using equal weights.") # = probability that insertion/roi is from the right genome is 1/number of genomes 
			return 1 / n_barcodes


	#Part 3: Simulate how many insertions can be found using our sequencing approach with different parameter settings
	@profile
	def count_insertions(insertion_dir, n_barcodes, read_dir):
		'''
		Counts the number of full-length and partial insertions for each insertion/roi. Genome scale factor due to diploidy of the human genome. 
		'''
		data = []
		print("Counting...")
		if scaling == True:
			genome_scale_facor = 0.5 
		else:
			genome_scale_facor = 1

		for insertion, insertion_positions in insertion_dir.items():
			full_length_count = 0
			partial_count = 0
			insertion_start, insertion_end = insertion_positions
			insertion_data = {'Insertion': insertion}

			for read, read_positions in read_dir.items():
				if read.split("_")[-1] == insertion.split("_")[1]: #if barcodes match
					read_start, read_end = read_positions

					if read_start <= insertion_start and insertion_end <= read_end:  # Full-length insertion
						if random.random() <= genome_scale_facor:
							full_length_count += 1
					elif (insertion_start < read_start and insertion_end > read_start) or \
						 (insertion_start < read_end and insertion_end > read_end):  # Partial insertion
						if random.random() <= genome_scale_facor:
							partial_count += 1
					
			insertion_data['full_matches'] = full_length_count
			insertion_data['partial_matches'] = partial_count
			data.append(insertion_data)
		return data



	def count_barcode_occurrences(dictionary):
		'''
		Counts how many reads are from which barcode. Sanity checks the weighted barcoding.
		'''
		suffix_count = {}
		for key in dictionary.keys():
			suffix = key.split("_")[-1]  # Get the last part of the key after "_"
			suffix_count[suffix] = suffix_count.get(suffix, 0) + 1  # Increment count for the suffix
		return suffix_count

	def normalize_ROI_by_length(roi_input_bed, roi_counted_insertions, scaling_factor=1):
		'''
		Normalizes found ROI matches by their length and normalizes by number of reads
		'''
		results_full = []
		results_partial = []
		results_combined=[]
		lengths=[]
		for _, row_counted in roi_counted_insertions.iterrows():
			# Initialize length to None
			length = None
			
			# Iterate through each row in df1 to find partial match
			for _, row_initial in roi_input_bed.iterrows():
				if row_initial.ID in row_counted.Insertion:
					# Calculate the length
					length = row_initial.end - row_initial.start
					lengths.append(length)
					break
			
			if length is not None:
				# Calculate the metric
				metric_full = ((row_counted['full_matches'] / length) * scaling_factor) / row_counted.Total_Reads #Reads per base (sf=1, kb sf=1000, or mb)
				results_full.append(metric_full)

				metric_partial = ((row_counted['partial_matches'] / length) * scaling_factor) / row_counted.Total_Reads #Reads per base (sf=1, kb sf=1000, or mb)
				results_partial.append(metric_partial)

				#combined ratio
				metric_combined = (((row_counted['partial_matches'] + row_counted['full_matches']) / length) * scaling_factor) / row_counted.Total_Reads #Reads per base (sf=1, kb sf=1000, or mb)
				results_combined.append(metric_combined)

			else:
				# If no partial match found, append None
				results_full.append(0)
				results_partial.append(0)
		
		# Add the results as a new column to DF2
		roi_counted_insertions['Length'] = lengths
		roi_counted_insertions['Full_Ratio'] = results_full
		roi_counted_insertions['Partial_Ratio'] = results_partial
		roi_counted_insertions['Combined_Ratio'] = results_combined
		
		return roi_counted_insertions

	def save_histogram(data, bins, mean_read_length, coverage):
		"""
		Generates a histogram from a list of numbers and saves it as a PNG file.
		"""
		mean_value = np.mean(data)
		median_value = np.median(data)
		nreads = len(data)
		# Create the histogram
		plt.figure(figsize=(10, 6))
		plt.hist(data, bins=bins, edgecolor='black', alpha=0.7)
		
		# Add a vertical line for the mean
		plt.axvline(mean_value, color='red', linestyle='dashed', linewidth=1)
		plt.text(mean_value, plt.ylim()[1]*0.9, f'Mean: {mean_value:.2f}', color='red')
		plt.axvline(median_value, color='blue', linestyle='dashed', linewidth=1)
		plt.text(median_value, plt.ylim()[1]*0.8, f'Median: {median_value:.2f}', color='blue')
		plt.text(plt.xlim()[1]*0.8, plt.ylim()[1]*0.8, f'N reads: {nreads}', color='blue')
		# Add labels and title
		plt.xlabel('Value')
		plt.ylabel('Frequency')
		plt.title('Histogram with Mean and Median')
		
		# Save the histogram as a PNG file
		output_file = f"{output_path_plots}/{experiment_name}_{mean_read_length}_{coverage}_histogram.png"
		plt.savefig(output_file, format='png', dpi=300)
		plt.close()

	def bin_coverage(coverage, bin_size):
		"""
		Aggregate coverage data into bins of a specified size.
		"""
		num_bins = int(np.ceil(len(coverage) / bin_size))
		binned_coverage = np.zeros(num_bins)
		
		for i in range(num_bins):
			start = i * bin_size
			end = start + bin_size
			binned_coverage[i] = np.mean(coverage[start:end])
		
		return binned_coverage

	@profile
	def plot_reads_coverage(reference_genome_length,bin_size, reads_dict, mean_read_length, current_coverage, insertion_dict):
		"""
		Plots a coverage-like plot using the reference genome and reads information.
		"""
		ref_length = reference_genome_length
		smooth_sigma = 3
		# Initialize coverage array
		coverage = np.zeros(ref_length)
		
		# Extract unique suffixes and assign colors
		unique_suffixes = set(read_id.split('_')[-1] for read_id in reads_dict.keys())
		colors = list(mcolors.TABLEAU_COLORS.keys())
		suffix_color_map = {suffix: colors[i % len(colors)] for i, suffix in enumerate(unique_suffixes)}
		
		# Populate coverage array based on reads
		for read_id, (start, stop) in reads_dict.items():
			suffix = read_id.split('_')[-1]
			coverage[start:stop] += 1
		
		# Bin the coverage
		binned_coverage = bin_coverage(coverage, bin_size)
		bin_positions = np.arange(0, ref_length, bin_size)
		
		# Smooth the binned coverage using a Gaussian filter
		smoothed_binned_coverage = gaussian_filter1d(binned_coverage, sigma=smooth_sigma)

		# Plot the binned coverage
		plt.figure(figsize=(20, 6))
		plt.plot(bin_positions[:len(smoothed_binned_coverage)], smoothed_binned_coverage, drawstyle='steps-pre', color='gray', alpha=0.5)
		#plt.plot(bin_positions[:len(binned_coverage)], binned_coverage, color='r', marker='o')
		# Plot individual reads with different colors based on suffix
		for suffix in unique_suffixes:
			coverage_suffix = np.zeros(ref_length)
			for read_id, (start, stop) in reads_dict.items():
				if read_id.endswith(suffix):
					coverage_suffix[start:stop] += 1
			binned_coverage_suffix = bin_coverage(coverage_suffix, bin_size)
			smoothed_binned_coverage_suffix = gaussian_filter1d(binned_coverage_suffix, sigma=smooth_sigma)
			plt.plot(bin_positions[:len(smoothed_binned_coverage_suffix)], smoothed_binned_coverage_suffix, drawstyle='steps-pre', color=suffix_color_map[suffix], label=suffix, alpha=0.5)

		# Add vertical lines at specified positions
		for key, positions in insertion_dict.items():
			suffix = key.split('_')[1]
			for pos in positions:
				plt.axvline(x=pos, color=suffix_color_map[suffix], linestyle='--')
		
		# Create custom legend
		legend_handles = [Line2D([0], [0], color=suffix_color_map[suffix], lw=2, label=suffix) for suffix in unique_suffixes]
		plt.legend(handles=legend_handles, title='Barcode', bbox_to_anchor=(1.05, 1), loc='upper left')
		
		plt.xlabel('Position on "one-string" Reference Genome (1e6 binned)')
		plt.ylabel('Read Coverage')
		plt.title('Read Coverage Plot')
		
		# Save the plot
		output_file = f"{output_path_plots}/{experiment_name}_{mean_read_length}_{current_coverage}_coverage.png"
		plt.savefig(output_file)
		plt.close()
		
		print(f"Plot saved as {output_file}")


	def process_combination(mean_read_length, coverage, length_mod_fasta, insertion_dir, iteration):
		'''
		This is the heart of the calculations, since it simulates the sequencing procedure. In brief:

		Creates a read length distribution based on mean read length and draws artificial reads from it until desired coverage is reached. 
		Then it compares the coordinates of the target regions (ROI or I) with the coordinates of artifical reads of each barcode and checks whether they are partially or fully contained.
		Masking is optional and is performed during the read generation step.
		'''
		#to prevent memory overfloat

		try:
			custom_read_length_distribution = get_read_length_distribution_from_real_data(sequenced_data_path) #for experimental data
			print("Custom FASTA data provided.")
			save_histogram(custom_read_length_distribution, 20, mean_read_length, coverage)
		except:
			print("No custom read length distribution provided... Generating artificial one...")
			custom_read_length_distribution = generate_read_length_distribution(num_reads=1000000, mean_read_length=mean_read_length, distribution='lognormal')
			print(custom_read_length_distribution)
			save_histogram(custom_read_length_distribution, 20, mean_read_length, coverage)
			print('Calculating for:' + str(mean_read_length))
		
		PRECOMPUTE_RANDOM = [random.choice(custom_read_length_distribution) for _ in range(100000000)]
		custom_cov_coordinates, covered_length = generate_reads_based_on_coverage(length_mod_fasta, custom_read_length_distribution, coverage, PRECOMPUTE_RANDOM)
		
		#plot coverage
		#plot_reads_coverage(length_mod_fasta, 1000000, custom_cov_coordinates, mean_read_length, coverage, insertion_dir) #might need to be commented out for high coverage runs! otherwise, the plot cannot be drawn!
		# Sanity check for barcoding
		barcode_distribution = count_barcode_occurrences(custom_cov_coordinates)
		barcode_distribution["coverage"] = coverage
		barcode_distribution["mean_read_length"] = mean_read_length
		barcode_distribution["iteration"] = iteration   
		barcode_distribution["N_bases"] = covered_length  #new, check if this works for roi too, although I'd be very surprised if not
		
		#Detections per ROI/Insertion
		detected = count_insertions(insertion_dir, n_barcodes, custom_cov_coordinates)
		suffix = "_" + str(iteration)
		for insertion_data in detected:
			insertion_data["mean_read_length"] = mean_read_length
			insertion_data["coverage"] = coverage
			insertion_data["Insertion"] += suffix

		return detected, barcode_distribution
		
		#save memory from overload
		del custom_read_length_distribution
		del custom_cov_coordinates
		gc.collect()

	def parallel_replicates(n_iterations): #x100 down to ~5.5h
			'''
			This snippet is th wrapper for the parallelization over the iterations. Each Iteration will be performed as a single job. 
			'''
			results=[]
			barcode_distributions=[]
			print("Iteration number: %i" % n_iterations)
			for mean_read_length, coverage in combinations:
				out, barcode_distribution = process_combination(mean_read_length, coverage, genome_size, target_regions, n_iterations)
				results.append(out)
				barcode_distributions.append(barcode_distribution)
			return results, barcode_distributions


	#insertion mode
	if mode == "I":
		'''
		I mode needs a quite complex pre-treatment, which is performed within the create_barcoded_insertion_genome function. In brief:
		1.) The reference genome cooridnates are transformed into a string-like format (One single FASTA string) and the chromsome borders are stored
		2.) Based on user input, masked regions are defined and weighted
		3.) For each barcode, the chromosome borders get a prefix and inserttion positions are randomly chosen. If an insertion-bed is defined, the insertions are only placed within these regions accoridng to their length-based weights.
		4.) The length of the reference genome (for coverage estimation), the insertion cooridnates, the masked regions, and the chromosome border dict are returned
		'''
		print("Insertion mode selected...")
		t0 = time.time()
		#run once
		#1 Creates the vector in the right format
		#insertion_fasta = collapse_fasta(vector_sequence_path)
		insertion_fasta = 'X' * 5000 #artificial insertion

		length_mod_fasta, insertion_dict, masked_regions, chromosome_dir = create_barcoded_insertion_genome(reference_genome_path=reference_genome_path, bedpath=bedpath, insertion_fasta=insertion_fasta, n_barcodes=n_barcodes)
		print("Number of insertions:")
		print(len(insertion_dict.keys()))
		print(insertion_dict)
		print(chromosome_dir)
		if masked_regions:
			masked_regions = {k: v for k, v in masked_regions.items() if len(v) == 3}
			print("Masked regions:")
			print(masked_regions)
		#get locations of the insertions
		insertion_locations_bed = convert_to_normal_bed(chromosome_dir, insertion_dict)
		print(insertion_locations_bed)

		#input variables
		genome_size = length_mod_fasta
		target_regions = insertion_dict

	#region of interest mode    
	elif mode == "ROI":
		'''
		ROI mode needs a different pre-treatment of the reference genome, as insertions are not spread randomly but have a fixed placement as defined in the bed file.
		'''
		print("Region of interest mode selected...")
		t0 = time.time()
		
		#create global cooridnates from bed based on provided genome ref to adjust ROIs to string-like genome
		fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
		bed = readbed(roi_bedpath, chromosome_dir.keys())
		print(bed.head())
		roi_dict = update_coordinates(bed, chromosome_dir)
		print(roi_dict)
		roi_dict = roi_barcoding(roi_dict, n_barcodes)
		
		#create blocked regions file
		if blocked_regions_bedpath is not None:
			blocked_bed = readbed(blocked_regions_bedpath, chromosome_dir.keys()) #barcoding should be default false
			masked_regions = update_coordinates(blocked_bed, chromosome_dir, include_weights=True)
		
		masked_regions = add_monosomy_regions(monosomie, chromosome_dir, masked_regions=masked_regions, chromosome_weights=chromosome_weights)
		print("masking...")
		print(masked_regions)
		#input variables
		genome_size = len(fasta)
		target_regions = roi_dict


	else:
		print("no valid mode selected.")
		sys.exit()  


	#Shared processing: The input files differ between ROI and I, but the downstream handling is the same!
	parallel_results= Parallel(n_jobs=parallel_jobs)(delayed(parallel_replicates)(i) for i in range(iterations))

	#the next part unpacks the nested structure of the parallel results
	results_list = []
	barcode_distributions_list = []

	# Iterate over parallel_results to extract results and barcode distributions
	for iteration_results, barcode_distributions in parallel_results:
		results_list.extend(iteration_results)
		barcode_distributions_list.append(barcode_distributions)

	# Flatten and covnert to df
	flattened_results_list = [item for sublist in results_list for item in sublist]
	results_df = pd.DataFrame(flattened_results_list)

	# Flatten and covnert to df
	flattened_barcode_distributions_list = list(itertools.chain(*barcode_distributions_list))
	barcode_distributions_df = pd.DataFrame(flattened_barcode_distributions_list)

	#ROI-sepcific transformations
	if mode == "ROI":
		'''
		For this mode, the files need some additonal modifications, namely the normalization of the detected ROIs by their length and total number of reads.
		The output dataframe is thereby much bigger than the insertion output. 
		'''

		# Split the 'Insertion' column in df1 to extract the 'iteration' number
		results_df['iteration'] = results_df['Insertion'].str.split('_').str[-1].astype(int)
		barcode_distributions_df['Total_Reads'] = barcode_distributions_df.iloc[:, :n_barcodes].sum(axis=1)
		#barcode_distributions_df = barcode_distributions_df.drop(columns=barcode_distributions_df.columns[:n_barcodes])
		results_df = pd.merge(results_df, barcode_distributions_df, on=['coverage', 'mean_read_length', 'iteration'], how='inner')
		results_df = normalize_ROI_by_length(bed, results_df, scaling_factor=0) #ROIPKB


	t1 = time.time()
	total = t1-t0
	print(total)
	print("Done.")

	#output file writing: Insertions have one mor eoutput file: The location bed of the random insertions
	try:
		insertion_locations_bed.to_csv(output_path+experiment_name+"_"+"insertion_locations.bed", sep='\t', header=True, index=False)
		print("Insertion mode: Writing files...")
		print(insertion_locations_bed.head())
		#### I pre processing
		results_df["Barcode"] = results_df["Insertion"].str.split("_").str[0] #changed, maybe causes problems with insertions?
		results_df["Iteration"] = results_df["Insertion"].str.split("_").str[4]
		combined_mean_df = results_df.groupby(['coverage', 'mean_read_length', 'Barcode']).agg({'full_matches': 'mean', 'partial_matches': 'mean'}).reset_index() #mean over iterations
		print("mean overview:")
		print(combined_mean_df)
	except:
		print("ROI mode: Writing files...")

	barcode_distributions_df.to_csv(output_path+experiment_name+"_"+"barcode_distribution_table.csv", sep='\t', header=True, index=False)
	results_df.to_csv(output_path+experiment_name+"_"+"matches_table.csv", sep='\t', header=True, index=False)

	print(f"Reference Genome Path: {reference_genome_path}")
	print(f"Vector Sequence Path: {vector_sequence_path}")
	print(f"Sequenced Data Path: {sequenced_data_path}")
	print(f"Output Path: {output_path}")
	print(f"Experiment Name: {experiment_name}")
	print(f"Output Path Plots: {output_path_plots}")
	print(f"Insertion Probability: {insertion_probability}")
	print(f"Chromosome Restriction: {chr_restriction}")
	print(f"Bed Path: {bedpath}")
	print(f"Barcode Weights: {barcode_weights}")
	print(f"Chromosome Weights: {chromosome_weights}")
	print(f"Insertion Numbers: {insertion_numbers}")
	print(f"Number of Barcodes: {n_barcodes}")
	print(f"Iterations: {iterations}")
	print(f"Scaling: {scaling}")
	print(f"Parallel Jobs: {parallel_jobs}")
	print(f"Mode: {mode}")
	print(f"Coverages: {coverages}")
	print(f"Mean Read Lengths: {mean_read_lengths}")
	print(f"ROI Bed Path: {roi_bedpath}")
	print(f"Blocked Regions Bed Path: {blocked_regions_bedpath}")
	print(f"Barcodes to Check Blocked Regions: {barcodes_to_check_blocked_regions}")
	print(f"Monosomie: {monosomie}")

if __name__ == "__main__":
	main()