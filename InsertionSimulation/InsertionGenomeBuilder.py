#!/usr/bin/env python3

from Bio import SeqIO
import random
import matplotlib.pyplot as plt
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

#will at some point be changes into argparse or sth similar
reference_genome_path = "/home/weichan/permanent/Projects/VIS/dev/VIS_Magdeburg_withBasecalling/hg38.fa" #reads will be created based on this reference
vector_sequence_path = "/home/weichan/permanent/Projects/VIS/dev/VIS_Magdeburg_withBasecalling/pSLCAR-CD19-28z.fasta"#vector #currently 8866 - 42000 (not observed in data): 5kb long should be enough!
sequenced_data_path = "/home/weichan/permanent/Projects/VIS/VIS_integration_site/Results/FullRunAfterModulaization_BUFFERMODE100_CD19_cd247_Vector_integration_site/FASTA/Full_MK025_GFP+.fa"
output_path = "./out/Homogeneous_NotWeigthed_20_Iterations_SummaryTable.csv"
bedpath = None #"/home/weichan/permanent/Projects/VIS/dev/UCSC/UCSC_GENCODE_V44_Introns_04_24" #default setting to None
weights_dict = {}#{"Barcode_0": 10, "Barcode_1": 5}
insertion_numbers=5
n_barcodes=1 #add function to set the default to 1 if barcoding = FALSE #doesn't work: barcoding is either tgrue and > 1 or false
iterations=20
parallel_jobs=10
#Part 1: Create an insertion-infiltrated chromosome and check where the insertions happen.

### WRAPPER START
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
	if n_barcodes >= 1:
		return True
	return False

def readbed(bedpath, list_of_chromosomes_in_reference, barcoding=False):
	try:
		bed = pd.read_csv(bedpath, sep='\t', header=None, usecols=[0,1,2], names=['chrom', 'start', 'end'])
		if barcoding:
			print('Barcoding selected: Transforming the names in the bed... Adding the following barcode...')
			print('_'.join(list(list_of_chromosomes_in_reference)[0].split("_")[:-1]))
			barcode='_'.join(list(list_of_chromosomes_in_reference)[0].split("_")[:-1])
			bed['chrom'] = barcode + '_'+ bed['chrom'].astype(str)
			bed = bed[bed["chrom"].isin(list_of_chromosomes_in_reference)]
		else:
			print('Only keeping the following chromosomes...')
			print(list_of_chromosomes_in_reference)
			bed = bed[bed["chrom"].isin(list_of_chromosomes_in_reference)]
		return bed
	except:
		return None

def collapse_fasta(path_to_fasta):
	with open(path_to_fasta, 'r') as fasta_file:
		# extracting multiple data in single fasta file using biopython
		seqList=[]
		for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
			seqList.append(str(record.seq))
	return ''.join(seqList)

def pseudo_fasta_coordinates(path_to_fasta):
	'''
	Takes in FASTA with separate entries and returns a table with the pseudo entry borders of each entry after the collapsing of the entries.
	'''
	with open(path_to_fasta, 'r') as fasta_file:
		entries={}
		seqList=[]
		updated_length=0
		for record in SeqIO.parse(fasta_file, 'fasta'):
			if '_' not in record.id and 'M' not in record.id:# and '2' in record.id: #the 2 modification will be removed later but is good for assessing the performance!
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

def get_chromosome(insert_position, chromosome_dir):
	'Checks where the random insertion landed in the artifical one string full genome'
	for chromosome, (start, end) in chromosome_dir.items():
		if start <= insert_position <= end:
			return chromosome
	return None

def insertions_per_chromosome(chromosome_dir, insertion_dir):
	"""
	Creates a table that lists the number of insertions per chromosome.
	"""
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

#new code for bed consideration option

#This is the basic functionality that needs to be implemented to only have bed file specific insertions based on the weighted coordinates
@profile
def add_insertions_to_genome_sequence_with_bed(reference_sequence, insertion_sequence, num_insertions, chromosome_dir, bed_df=None): #costs 3 billion per barcode!
	"""Randomly add insertion sequence into the reference genome or within specified regions."""
	position = {}
	if bed_df is not None:
		# Step 1: Calculate probabilities based on region lengths
		print("Calculating insertion probabilities (region length / sum of all regions lengths)...")
		region_lengths = bed_df['end'] - bed_df['start']
		region_probabilities = region_lengths / region_lengths.sum()
		updated_reference_sequence = reference_sequence

		for i in range(num_insertions):
			# Step 2: Randomly select insertion regions
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
	#return reference_sequence, position
	#updated_reference_sequence, position = add_insertions_to_genome_sequence(reference_sequence=reference_sequence, insertion_sequence=insertion_sequence, num_insertions=num_insertions)
	return updated_reference_sequence, position


#1 Creates the vector in the right format
#insertion_fasta = collapse_fasta(vector_sequence_path)
insertion_fasta = 'X' * 5000
'''
#2 Creates a single-string reference genome, while storing the chromosome borders
fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path, barcoding=barcoding, barcode=5)
bed_df = readbed(bedpath, chromosome_dir.keys(), barcoding=barcoding)
#3 Randomly inserts insertio sequence into the single-string reference x times
mod_fasta, insertion_dir = add_insertions_to_genome_sequence_with_bed(fasta, insertion_fasta, insertion_numbers, chromosome_dir, bed_df)
#4 Creates an overview table of the insertion distribution
print(insertion_dir)
insertions_df = insertions_per_chromosome(chromosome_dir, insertion_dir)
print(insertions_df)
'''
@profile
def create_barcoded_insertion_genome(reference_genome_path, bedpath, insertion_fasta, n_barcodes):
	"""
	Calls all functions n times and combines their results in one barcoded output file. This allows us to simulate a mixture of different cells.
	""" 
	print("Create Barcoded Genome...")
	#mod_fasta_dict={} #list of fastas
	insertion_dict={} #list of dirs
	
	fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
	
	for i in range(n_barcodes):
		barcoded_chromosome_dir = barcode_genome(chromosome_dir, i)
		bed_df = readbed(bedpath, barcoded_chromosome_dir.keys(), barcoding=check_barcoding(n_barcodes))
		mod_fasta, insertion_dir = add_insertions_to_genome_sequence_with_bed(fasta, insertion_fasta, insertion_numbers, barcoded_chromosome_dir, bed_df)
		insertion_dict.update(insertion_dir)
		#mod_fasta_dict.update({str(i): mod_fasta})
		del barcoded_chromosome_dir, bed_df, insertion_dir
	
	return len(mod_fasta), insertion_dict #''.join( #_dict #now returns just one mod fasta


#Part 2: Create artifical reads
def get_read_length_distribution_from_real_data(path_to_sequenced_fasta):
	lengths=[]
	for record in SeqIO.parse(path_to_sequenced_fasta, 'fasta'):
			lengths.append(len(record.seq))
	return lengths

def generate_reads(fasta, read_length_distribution, num_reads):
	"""Generate reads based on custom read length distribution and pulled from fasta ref."""
	reads = []
	for read_length in read_length_distribution:
		start_position = random.randint(0, len(fasta) - read_length)
		read = fasta[start_position:start_position + read_length]
		reads.append(read)
	return reads

@profile
def generate_reads_based_on_coverage(fasta, read_length_distribution, coverage, PRECOMPUTE_RANDOM):
	'''
	Randomly pulls a read of size X derived from the read length distribution from the fasta until the fasta is N times covered (coverage). 
	'''
	print("Coverage: " + str(coverage))
	print("Pulling reads...")
	covered_length = 0
	total_length = fasta
	#reads = [] #only a goodf idea if we are not testing high coverages, otherwise memory is floated
	read_coordinates = {}
	while covered_length < coverage * total_length: #human genome size now #total_length: #6 seconds per 3x10^9 string length for mean read length of 1k
		#read_length = random.choice(read_length_distribution)
		read_length = PRECOMPUTE_RANDOM.pop()
		start_position = random.randint(0, total_length - read_length)
		#read = fasta[start_position:start_position + read_length]
		#reads.append(read)
		covered_length += read_length
		read_coordinates["Read_%s" % len(read_coordinates.keys())] =[start_position, start_position + read_length]
	return read_coordinates #reads optional



#@timeit
def generate_read_length_distribution(num_reads, mean_read_length, distribution='lognormal', **kwargs):
	"""
	Generate a list of read lengths with a customizable mean read length and distribution.
	"""
	if distribution == 'normal':
		# Generate read lengths from a normal distribution
		read_lengths = np.random.normal(mean_read_length, mean_read_length / 10, num_reads)
	elif distribution == 'lognormal':
		# Generate read lengths from a log-normal distribution
		read_lengths = np.random.lognormal(mean=np.log(mean_read_length), sigma=1.0, size=num_reads)
	else:
		raise ValueError("Unsupported distribution. Supported options: 'normal', 'lognormal'.")
	
	# Ensure that all read lengths are positive integers
	read_lengths = np.round(np.abs(read_lengths)).astype(int)
	# Filter out zero-length reads
	read_lengths = read_lengths[read_lengths > 0]
	
	return read_lengths #not returning list


def get_weighted_probabilities(insertion_name,n_barcodes, weights_dict):
	"""
	Checks in the insertion name for a key of weights dir and returns the weighting factor.
	This is then used to increase the probability
	"""
	if weights_dict is not None:
		#print("Add weights to barcodes for ...")
		#print(insertion_name)
		# Calculate the common denominator
		common_denominator = (sum(weights_dict.values()) + n_barcodes - len(weights_dict)) * n_barcodes
		# Check if insertion_name contains any key in weights_dict
		for key in weights_dict: 
			if key in insertion_name:
				#print(key)
				#print((weights_dict[key] * n_barcodes) / common_denominator)
				return (weights_dict[key] * n_barcodes) / common_denominator
		# No weight provided, using 1/Number of barcodes for this barcode
		#print("No weight provided. Using standard weight for this barcode. Standard weigth:")
		#print(n_barcodes / common_denominator)
		return n_barcodes / common_denominator
	else:
		#print("No weights provided, using equal weights.")
		return 1 / n_barcodes

#Part 3: Simulate how many insertions can be found using our sequencing approach with different parameter settings
@profile
def count_insertions(insertion_dir, n_barcodes,weights_dict, read_dir):
	"""
	Count the number of full-length and partial insertions for each insertion.
	"""
	data = []
	print("Counting...")

	for insertion, insertion_positions in insertion_dir.items():
		full_length_count = 0
		partial_count = 0
		insertion_start, insertion_end = insertion_positions
		insertion_data = {'Insertion': insertion}

		weight = get_weighted_probabilities(insertion, n_barcodes, weights_dict) #get probability of insertion/genome: This probability will then be used for eac read!
		

		'''
		The weights are supposed to illustrate the follwoing probability: 
		If I choose a random read from my collection, how probable is it that this read contains an insertion and both of them are from the same genome.
		'''

		 #we have a certain probability that the read we pulled is from the correct genome. This probability is defined in the weighting function accoridng to the assigned weigths.

		for read, read_positions in read_dir.items():
			if weight >= random.random():
				read_start, read_end = read_positions

				if read_start <= insertion_start and insertion_end <= read_end:  # Full-length insertion
					full_length_count += 1
				elif (insertion_start < read_start and insertion_end > read_start) or \
					 (insertion_start < read_end and insertion_end > read_end):  # Partial insertion
					partial_count += 1
				
		insertion_data['full_matches'] = full_length_count
		insertion_data['partial_matches'] = partial_count
		data.append(insertion_data)
		

	return data

'''
### Calling everything in the right order
#1 Creates the vector in the right format
insertion_fasta = collapse_fasta(vector_sequence_path)
#2 Creates a single-string reference genome, while storing the chromosome borders
fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
#3 Randomly inserts insertio sequence into the single-string reference x times
mod_fasta, insertion_dir = add_insertions_to_genome_sequence(fasta, insertion_fasta, insertion_numbers)
#4 Creates an overview table of the insertion distribution
insertions_df = insertions_per_chromosome(chromosome_dir, insertion_dir)
print(insertions_df)
#5 creates histogram from actual sequenced distribution and the randomly pulled reads
sequence_read_length_distribution = get_read_length_distribution_from_real_data(sequenced_data_path)
plt.hist(sequence_read_length_distribution, bins=400)
reads = generate_reads(mod_fasta, sequence_read_length_distribution, num_reads=len(sequence_read_length_distribution))
plt.hist([len(i) for i in reads], bins=400)

cov_reads, cov_coordinates = generate_reads_based_on_coverage(mod_fasta, sequence_read_length_distribution, coverage=coverage)
plt.hist([len(i) for i in cov_reads], bins=400)
plt.savefig('sequence_based_read_length_distribution.pdf')
plt.close()

#6 creates histogram for artificially created distribution and the randomly pulled reads
custom_read_length_distribution = generate_read_length_distribution(500000, mean_read_length = mean_read_length, distribution = 'lognormal') #currently only normal and lognormal work
plt.hist(custom_read_length_distribution, bins=400)

custom_cov_reads, custom_cov_coordinates = generate_reads_based_on_coverage(mod_fasta, custom_read_length_distribution, coverage=coverage)
plt.hist([len(i) for i in custom_cov_reads], bins=400)
plt.savefig('calculated_read_length_distribution.pdf')
#7 creates final table with detected insertions
detected, full, partial = count_insertions(insertion_dir, custom_cov_coordinates)
print("Total Reads: " + str(len(custom_cov_reads)))
print(full)
print(partial)
'''
t0 = time.time()
coverages = [1, 5, 10, 15, 20] 
mean_read_lengths = [1000, 2000, 3000, 4000, 5000, 6000,7000,8000,9000,10000,15000,20000]
#mean_read_lengths = [5000, 8000, 12000]
#coverages=[1,5,10,15] #* 10 #,5,10] #* 10 #coverage with some influence on the runtime
combinations = itertools.product(mean_read_lengths, coverages)
#coverages=[1,2]
#mean_read_lengths=[5000, 6000]

#run once
#1 Creates the vector in the right format
#insertion_fasta = collapse_fasta(vector_sequence_path)
insertion_fasta = 'X' * 5000 #artificial insertion

'''
#2 Creates a single-string reference genome, while storing the chromosome borders
fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
#3 import bed with possible insertion positions 
bed_df = readbed(bedpath, chromosome_dir.keys())
#3 Randomly inserts insertio sequence into the single-string reference x times
mod_fasta, insertion_dir = add_insertions_to_genome_sequence_with_bed(fasta, insertion_fasta, insertion_numbers, chromosome_dir, bed_df)

del fasta
#del insertion_fasta
del bed_df
'''
# Define a function to process each combination
def process_combination(mean_read_length, coverage, insertion_numbers, weights_dict, length_mod_fasta, insertion_dir, iteration):
	
	#Creates a read length distribution based mean read length and draws artificial reads from it based on coverage settings. Then it compares the coordinates of the previously created
	#insertions in the reference fasta and checks whether they are partially or fully contained within the artificial reads.
	
	#custom_read_length_distribution = get_read_length_distribution_from_real_data(sequenced_data_path) #for experimental data
	print('Calculating for:' + str(mean_read_length))
	custom_read_length_distribution = generate_read_length_distribution(1000000, mean_read_length=mean_read_length, distribution='lognormal')
	PRECOMPUTE_RANDOM = [random.choice(custom_read_length_distribution) for _ in range(100000000)]
	#custom_cov_coordinates = generate_reads_based_on_coverage_parallel(mod_fasta, custom_read_length_distribution, coverage=coverage, num_processes=10)
	custom_cov_coordinates = generate_reads_based_on_coverage(length_mod_fasta, custom_read_length_distribution, coverage, PRECOMPUTE_RANDOM)
	#detected, full, partial = count_insertions(insertion_dir, custom_cov_coordinates)
	detected = count_insertions(insertion_dir, n_barcodes, weights_dict, custom_cov_coordinates)
	suffix = "_" + str(iteration)
	for insertion_data in detected:
		insertion_data["mean_read_length"] = mean_read_length
		insertion_data["coverage"] = coverage
		insertion_data["Insertion"] += suffix
	#return [mean_read_length, coverage, insertion_numbers, full, partial]
	return detected
	#save memory from overload
	del custom_read_length_distribution
	del custom_cov_coordinates
	gc.collect()

#barcode test
length_mod_fasta, insertion_dict = create_barcoded_insertion_genome(reference_genome_path=reference_genome_path, bedpath=bedpath, insertion_fasta=insertion_fasta, n_barcodes=n_barcodes)
print(insertion_dict)

def parallel_replicates(n_interations): #x100 down to ~5.5h
	results=[]
	print("Iteration number: %i" % n_interations)
	for mean_read_length, coverage in combinations:
		out = process_combination(mean_read_length, coverage, insertion_numbers, weights_dict, length_mod_fasta, insertion_dict, n_interations)
		results.append(out)
	return results

parallel_results= Parallel(n_jobs=parallel_jobs)(delayed(parallel_replicates)(i) for i in range(iterations))
finaloutput = list(itertools.chain(*parallel_results))

#results=[]

#for mean_read_length, coverage in combinations:
#    out = process_combination(mean_read_length, coverage, insertion_numbers, big_fasta, big_insertion_dict)
#    results.append(out)

#print(results)
#sys.exit()
# Convert finaloutput to DataFrame
if any(isinstance(i, list) for i in finaloutput):
	finaloutput_df = pd.DataFrame(sum(finaloutput, []))
#finaloutput_df = pd.DataFrame(finaloutput, columns=['mean_read_length', 'coverage', 'insertion_number', 'full_matches', 'partial_matches'])
t1 = time.time()
total = t1-t0
print(total)
print(finaloutput_df)
# Write DataFrame to CSV
finaloutput_df.to_csv(output_path, sep='\t', header=True, index=False)
print("Done.")
sys.exit()