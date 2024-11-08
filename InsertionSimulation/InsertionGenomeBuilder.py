#%%
#!/usr/bin/env python3

from Bio import SeqIO
import random
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d
import math
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

# class-based
from config_handler import ConfigHandler


def main():
	
	# Set the desired section dynamically
	section = 'I'  # or 'I' depending on your use case

	# Initialize ConfigHandler with the chosen section
	config_handler = ConfigHandler(config_file="sim_config.ini", section=section)
	param_dictionary = config_handler.parse_config()

	mean_read_lengths = param_dictionary.get('mean_read_lengths')
	coverages = param_dictionary.get('coverages')
	
	if type(mean_read_lengths) is not list:
		mean_read_lengths = [mean_read_lengths]
	if type(coverages) is not list:
		coverages = [coverages]

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

	def bedcheck(bed):
		"""
		Check if the first column of a df contains the string "chr" and if the second and thirs column contain an integer.  
		"""
		if bed.iloc[:, 0].str.startswith("chr").all() and \
			bed.iloc[:, 1].apply(lambda x: isinstance(x, int)).all() and \
			bed.iloc[:, 2].apply(lambda x: isinstance(x, int)).all():
			return True
		
		else:
			print("The file does not seem to be in BED format. Make sure the data looks like: chrN integer integer")
			return False

	def readbed(bedpath, list_of_chromosomes_in_reference, barcoding=False):
		"""
		Reads bed files in a flexible manner. Adds None fr columns that are non-existent. 
		"""
		bed = None
		print(bedpath)

		try:
			# Read only the available columns initially
			bed = pd.read_csv(bedpath, sep='\t', header=0)
			print(bed.head())
			if not bedcheck(bed):
				return None
			
			all_cols = ["chrom", "start", "end", "ID", "Barcode", "weight"]
			missing_cols = list(set(all_cols) - set(bed.columns))
			print("Your BED is missing the following columns. Trying to fill in default values where possible...")
			print(missing_cols)
			for missing_col in missing_cols: 
				if missing_col == "weight":
					bed[missing_col] = 1 
				elif missing_col == "ID":
					bed[missing_col] = bed["chrom"].astype(str) + "_" + bed["start"].astype(str) + "_" + bed["end"].astype(str) + "_" + bed["weight"] .astype(str)
				else:
					bed[missing_col] = [[]] * len(bed)
			
			# Fill empty or NaN values in the 'weight' column with the default value of 1
			bed['weight'] = bed['weight'].apply(lambda x: 1 if x in ["", None] or (isinstance(x, float) and math.isnan(x)) else x)
		
			# Fill NaN values in the 'Barcode' column with an empty list
			bed['Barcode'] = bed['Barcode'].apply(lambda x: [] if x in ["", None] or (isinstance(x, float) and math.isnan(x)) else x)
			
			
			# Apply barcoding transformation if required
			if barcoding:
				print('Barcoding selected: Transforming the chromosome names in the bed...')
				barcode = '_'.join(list(list_of_chromosomes_in_reference)[0].split("_")[:-1])
				bed['chrom'] = barcode + '_' + bed['chrom'].astype(str)

			# Filter out chromosomes not in the reference list
			print('Only keeping the following chromosomes:', list_of_chromosomes_in_reference)
			bed = bed[bed["chrom"].isin(list_of_chromosomes_in_reference)]
					
		except Exception as e:
			print("Error reading the BED file:", e)
			return None
		
		return bed

	def update_coordinates(df, input_chromosome_dict):
		'''
		Uses a df (bed-like) and a dict of chromosome broders in global format and returns the global cooridnates of the bed entries.
		Checks if start and stop are 0s and if so, sets the coordinates to the full chromosome length for this row. 
		'''
		try:
			
			updated_coordinates = {}
			barcode_exists = 'Barcode' in df.columns
			weight_exists = 'weight' in df.columns
			
			for index, row in df.iterrows():
				ID=row["ID"]
				chrom = row['chrom']
				
				if row["start"] == row["end"] == 0: #full chromosome mode
					print(f"Whole {chrom} will be used for {ID}.")
					start = input_chromosome_dict[chrom][0]
					end = input_chromosome_dict[chrom][1]
				else:
					start = row['start'] + input_chromosome_dict[chrom][0]
					end = row['end'] + input_chromosome_dict[chrom][0]

				# Prepare the entry for updated_coordinates based on conditions
				entry = {'start': start, 'end': end, 'weight': row['weight'], 'Barcode': row['Barcode']}
				updated_coordinates[ID] = entry
			return updated_coordinates
		except:
			print("BED coordinates could not be updated.")
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
				if param_dictionary.get('chr_restriction') == "unrestricted":
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

	def insertions_per_chromosome(chromosome_dir, insertion_dict):
		'''
		Creates a table that lists the number of insertions per chromosome.
		'''
		
		# Convert chromosome_dir and insertion_dict to sets for faster lookups
		chromosome_set = {chromosome: tuple(coordinates) for chromosome, coordinates in chromosome_dir.items()}
		insertion_set = {insertion: tuple(coordinates) for insertion, coordinates in insertion_dict.items()}
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
		
		print(f"insertion_number_distribution: {param_dictionary.get('insertion_number_distribution')}")

		if param_dictionary.get('insertion_number_distribution') == 'poisson': 
			num_insertions = np.random.poisson(num_insertions)
			print(f"Number of insertions drawn from Poisson distribution: {num_insertions}")
		else:
			print(f"Using exactly {num_insertions}.")

		if bed_df is not None:
			print("BED guided insertion pattern...")
			print(bed_df.head())
			# Step 1: Calculate probabilities based on region lengths
			print("Calculating insertion probabilities (region length / sum of all regions lengths)...")
			region_lengths = bed_df['end'] - bed_df['start']
			region_probabilities = region_lengths / region_lengths.sum()
			updated_reference_sequence = reference_sequence

			for i in range(num_insertions):
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

		print(updated_reference_sequence)
		print(position)
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
		collected_insertion_dict={}
		#1
		fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
		
		#2 #create blocked regions file
		if param_dictionary.get('blocked_regions_bedpath') is not None:
			blocked_bed = readbed(param_dictionary.get('blocked_regions_bedpath'), chromosome_dir.keys())
			masked_regions = update_coordinates(blocked_bed, chromosome_dir)
		else:
			masked_regions = None
		
		print(masked_regions)

		#3
		for i in range(n_barcodes):
			barcoded_chromosome_dir = barcode_genome(chromosome_dir, i)
			#optional bed-guided insertion
			bed_df = readbed(bedpath, barcoded_chromosome_dir.keys(), barcoding=check_barcoding(n_barcodes))
			mod_fasta, insertion_dict = add_insertions_to_genome_sequence_with_bed(fasta, insertion_fasta, param_dictionary.get('insertion_numbers'), barcoded_chromosome_dir, bed_df)
			collected_insertion_dict.update(insertion_dict)
			del barcoded_chromosome_dir, bed_df, insertion_dict
		#4
		return len(mod_fasta), collected_insertion_dict, masked_regions, chromosome_dir

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

	def check_if_blocked_region(random_barcode_number, start_position, read_length, masked_regions=None):
		'''
		Checks whether the start position is inside the coordinates of a masked/blocked region.
		Blocked can mean that no reads are obtained from there or only with a certain probability.
		'''
		for key, values in masked_regions.items():
			# Convert Barcode string to list if needed (e.g., '[]' to [])
			barcodes = eval(values["Barcode"]) if isinstance(values["Barcode"], str) else values["Barcode"]
			if barcodes is None:
				barcodes = []
			if int(random_barcode_number) in barcodes or not barcodes:
				start = values["start"]
				stop = values["end"]
				weight = values["weight"] #full blockage if weight not provided
				
				if (start < start_position < stop) or (start < start_position + read_length < stop): #checks if start lays in the blocked region or if the end lays in a blocked region
					# The start position falls within a blocked region
					if random.random() < weight:
						return True  # position is blocked 

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
		barcode_names = ["Barcode_" + str(i) for i in range(param_dictionary.get('n_barcodes'))]

		#Reads pulled until desired coverage is reached
		while covered_length < coverage * total_length:
			random_barcode = random.choices(barcode_names, [get_weighted_probabilities(i, param_dictionary.get('n_barcodes'), param_dictionary.get('barcode_weights')) for i in barcode_names])[0] #chooses one of the barcodes based on weighted probability
			read_length = PRECOMPUTE_RANDOM.pop()
			start_position = random.randint(0, total_length - read_length)
			# Check if the random barcode is in the list of barcodes that require checking for blocked regions
			random_barcode_number = random_barcode.split('_')[-1] #otherwise user input needs to be weird


			# Check if the start position falls within a blocked region
			if masked_regions and check_if_blocked_region(random_barcode_number, start_position, read_length, masked_regions):
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


	@profile
	def count_insertions(insertion_dict, read_dir):
		'''
		Counts the number of full-length and partial insertions for each insertion/roi.
		Requires at least 'min_overlap' between the insertion and read intervals
		Genome scale factor due to diploidy of the human genome.
		'''
		data = []
		print("Counting...")
		if param_dictionary.get('scaling') == True:
			genome_scale_facor = 0.5 
		else:
			genome_scale_facor = 1

		for key, values in insertion_dict.items():
			full_length_count = 0
			partial_count = 0
			#necessary since Insertions and ROIs have a different structure
			if type(values) == list: 
				start = values[0]
				end = values[1]
			else:
				start = values["start"]
				end = values["end"]
			countdata = {'Insertion': key}
			overlaps=[]
			for read, read_positions in read_dir.items():
				if read.split("_")[-1] == key.split("_")[1]:  # If barcodes match
					read_start, read_end = read_positions

					# Check if there is overlap first
					if max(start, read_start) < min(end, read_end):
						
						# Calculate overlap between the insertion and the read
						overlap = min(end, read_end) - max(start, read_start)
						
						# Full-length insertion: check if the read fully covers the insertion
						if read_start <= start and end <= read_end:
							if overlap >= param_dictionary.get('min_overlap_for_detection'):  # Ensure the overlap is at least 'x'
								if random.random() <= genome_scale_facor:
									full_length_count += 1
									overlaps.append(overlap)


						# Partial insertion: check if the read partially overlaps with the insertion
						elif (start < read_start and end > read_start) or \
							 (start < read_end and end > read_end):
							if overlap >= param_dictionary.get('min_overlap_for_detection'):  # Ensure the overlap is at least 'x'
								if random.random() <= genome_scale_facor:
									partial_count += 1
									overlaps.append(overlap)


			countdata['full_matches'] = full_length_count
			countdata['partial_matches'] = partial_count
			countdata['overlap'] = overlaps
			data.append(countdata)

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
		
		#plt.show()
		# Save the histogram as a PNG file
		output_file = f"{param_dictionary.get('output_path_plots')}/{param_dictionary.get('experiment_name')}_{mean_read_length}_{coverage}_histogram.png"
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
		plt.show()
		# Save the plot
		output_file = f"{param_dictionary.get('output_path_plots')}/{param_dictionary.get('experiment_name')}_{mean_read_length}_{current_coverage}_coverage.png"
		plt.savefig(output_file)
		plt.close()
		
		print(f"Plot saved as {output_file}")


	def process_combination(mean_read_length, coverage, length_mod_fasta, insertion_dict, iteration):
		'''
		This is the heart of the calculations, since it simulates the sequencing procedure. In brief:

		Creates a read length distribution based on mean read length and draws artificial reads from it until desired coverage is reached. 
		Then it compares the coordinates of the target regions (ROI or I) with the coordinates of artifical reads of each barcode and checks whether they are partially or fully contained.
		Masking is optional and is performed during the read generation step.
		'''
		#to prevent memory overfloat

		try:
			custom_read_length_distribution = get_read_length_distribution_from_real_data(param_dictionary.get('sequenced_data_path')) #for experimental data
			print("Custom FASTA data provided.")
			#save_histogram(custom_read_length_distribution, 20, mean_read_length, coverage)
		except:
			print("No custom read length distribution provided... Generating artificial one...")
			custom_read_length_distribution = generate_read_length_distribution(num_reads=1000000, mean_read_length=mean_read_length, distribution='lognormal')
			#save_histogram(custom_read_length_distribution, 20, mean_read_length, coverage)
			print('Calculating for:' + str(mean_read_length))
		
		PRECOMPUTE_RANDOM = [random.choice(custom_read_length_distribution) for _ in range(100000000)]
		custom_cov_coordinates, covered_length = generate_reads_based_on_coverage(length_mod_fasta, custom_read_length_distribution, coverage, PRECOMPUTE_RANDOM)
		#plot coverage
		#plot_reads_coverage(length_mod_fasta, 1000000, custom_cov_coordinates, mean_read_length, coverage, insertion_dict) #might need to be commented out for high coverage runs! otherwise, the plot cannot be drawn!
		# Sanity check for barcoding
		barcode_distribution = count_barcode_occurrences(custom_cov_coordinates)
		barcode_distribution["coverage"] = coverage
		barcode_distribution["mean_read_length"] = mean_read_length
		barcode_distribution["iteration"] = iteration   
		barcode_distribution["N_bases"] = covered_length  #new, check if this works for roi too, although I'd be very surprised if not
		
		#Detections per ROI/Insertion
		detected = count_insertions(insertion_dict, custom_cov_coordinates)
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

	def parallel_replicates(n_iterations): #x100 faster
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
	if param_dictionary.get('mode') == "I":
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
		insertion_fasta = 'X' * param_dictionary.get('insertion_length')

		length_mod_fasta, insertion_dict, masked_regions, chromosome_dir = create_barcoded_insertion_genome(reference_genome_path=param_dictionary.get('reference_genome_path'), bedpath=param_dictionary.get('bedpath'), insertion_fasta=insertion_fasta, n_barcodes=param_dictionary.get('n_barcodes'))
		print("Number of insertions:")
		print(len(insertion_dict.keys()))
		print(insertion_dict)
		print(chromosome_dir)
		if masked_regions:
			print("Masked regions:")
			print(masked_regions)
		#get locations of the insertions
		insertion_locations_bed = convert_to_normal_bed(chromosome_dir, insertion_dict)
		print(insertion_locations_bed)

		#input variables
		genome_size = length_mod_fasta
		target_regions = insertion_dict

	#region of interest mode    
	elif param_dictionary.get('mode') == "ROI":
		'''
		ROI mode needs a different pre-treatment of the reference genome, as ROIs are not spread randomly but have a fixed placement as defined in the bed file.
		'''
		print("Region of interest mode selected...")
		t0 = time.time()
		
		#create global cooridnates from bed based on provided genome ref to adjust ROIs to string-like genome
		fasta, chromosome_dir = pseudo_fasta_coordinates(param_dictionary.get('reference_genome_path'))
		bed = readbed(param_dictionary.get('roi_bedpath'), chromosome_dir.keys())
		roi_dict = update_coordinates(bed, chromosome_dir)
		roi_dict = roi_barcoding(roi_dict, param_dictionary.get('n_barcodes'))
		
		#create blocked regions file
		if param_dictionary.get('blocked_regions_bedpath') is not None:
			blocked_bed = readbed(param_dictionary.get('blocked_regions_bedpath'), chromosome_dir.keys()) #barcoding should be default false
			masked_regions = update_coordinates(blocked_bed, chromosome_dir)
			print("Masked regions:")
			print(masked_regions)
		
		#input variables
		genome_size = len(fasta)
		target_regions = roi_dict


	else:
		print("no valid mode selected.")
		sys.exit()  


	#Shared processing: The input files differ between ROI and I, but the downstream handling is the same!
	parallel_results= Parallel(n_jobs=param_dictionary.get('parallel_jobs'))(delayed(parallel_replicates)(i) for i in range(param_dictionary.get('iterations')))

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
	if param_dictionary.get('mode') == "ROI":
		'''
		For this mode, the files need some additonal modifications, namely the normalization of the detected ROIs by their length and total number of reads.
		The output dataframe is thereby much bigger than the insertion output. 
		'''

		# Split the 'Insertion' column in df1 to extract the 'iteration' number
		results_df['iteration'] = results_df['Insertion'].str.split('_').str[-1].astype(int)
		barcode_distributions_df['Total_Reads'] = barcode_distributions_df.iloc[:, :param_dictionary.get('n_barcodes')].sum(axis=1)
		#barcode_distributions_df = barcode_distributions_df.drop(columns=barcode_distributions_df.columns[:n_barcodes])
		results_df = pd.merge(results_df, barcode_distributions_df, on=['coverage', 'mean_read_length', 'iteration'], how='inner')
		results_df = normalize_ROI_by_length(bed, results_df, scaling_factor=0) #ROIPKB


	t1 = time.time()
	total = t1-t0
	print(total)
	print("Done.")

	#output file writing: Insertions have one mor eoutput file: The location bed of the random insertions
	try:
		insertion_locations_bed.to_csv(param_dictionary.get('output_path')+param_dictionary.get('experiment_name')+"_"+"insertion_locations.bed", sep='\t', header=True, index=False)
		print("Insertion mode: Writing files...")
		print(insertion_locations_bed.head())
		#### I pre processing
		results_df["Barcode"] = results_df["Insertion"].str.split("_").str[0] #changed, maybe causes problems with insertions?
		results_df["Iteration"] = results_df["Insertion"].str.split("_").str[4]
		print(results_df.tail())
	except:
		print("ROI mode: Writing files...")
		print(results_df.head())

	barcode_distributions_df.to_csv(param_dictionary.get('output_path')+param_dictionary.get('experiment_name')+"_"+"barcode_distribution_table.csv", sep='\t', header=True, index=False)
	results_df.to_csv(param_dictionary.get('output_path')+param_dictionary.get('experiment_name')+"_"+"matches_table.csv", sep='\t', header=True, index=False)
	
	print("Input parameters for the simulation:")
	
	# Now param_dictionary is accessible with correctly typed values
	for param, value in param_dictionary.items():
		print(f"{param} = {value} ({type(value)})")
	
if __name__ == "__main__":
	try:
		main()
	except ValueError as e:
		print("Configuration error:", e)
		sys.exit(1)
# %%


