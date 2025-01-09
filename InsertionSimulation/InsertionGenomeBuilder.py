#%%	
#!/usr/bin/env python3
from Bio import SeqIO
import random
import numpy as np
import pandas as pd
import sys
import time
import gc
import itertools
from functools import partial
from joblib import Parallel, delayed

# custom config reader
from config_handler import parse_config

#custom functions
#from create_insertion_genome import add_insertions_to_genome_sequence_with_bed
#from utils import profile, check_barcoding, roi_barcoding, barcode_genome
from bed_operations import global_to_chromosome_coordinates
from counting import normalize_ROI_by_length
from genome_generation import create_barcoded_insertion_genome, create_barcoded_roi_genome
from combined_calculations import process_combination

def main():
	# config
	config_file = "sim_config_roi.ini"
	try:
		param_dictionary = parse_config(config_file)
		print(param_dictionary)
	except Exception as e:
		print(f"Unexpected error: {e}")
	
	def parallel_replicates(n_iterations): #x100 faster
			'''
			This snippet is th wrapper for the parallelization over the iterations. Each Iteration will be performed as a single job. 
			'''
			results=[]
			barcode_distributions=[]
			print("Iteration number: %i" % n_iterations)
			for mean_read_length, coverage in param_dictionary.get('combinations'):
				out, barcode_distribution = process_combination(mean_read_length, coverage, genome_size, 
													target_regions, n_iterations, 
													param_dictionary.get('sequenced_data_path'), 
													param_dictionary.get('n_barcodes'), 
													param_dictionary.get('barcode_weights'), 
													masked_regions, 
													param_dictionary.get('output_path'),
													param_dictionary.get('scaling'), 
													param_dictionary.get('min_overlap_for_detection'))
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
		genome_size, insertion_dict, masked_regions, chromosome_dir = create_barcoded_insertion_genome(param_dictionary.get('reference_genome_path'), 
																								 param_dictionary.get('bedpath'),
																								 param_dictionary.get('blocked_regions_bedpath'),
																								 param_dictionary.get('chr_restriction'), 
																								 param_dictionary.get('insertion_length'),
																								 param_dictionary.get('insertion_numbers'),
																								 param_dictionary.get('insertion_number_distribution'),
																								 param_dictionary.get('n_barcodes'))
		print("Number of insertions:")
		print(len(insertion_dict.keys()))
		print(insertion_dict)
		print(chromosome_dir)
		if masked_regions:
			print("Masked regions:")
			print(masked_regions)
		#get locations of the insertions
		insertion_locations_bed = global_to_chromosome_coordinates(chromosome_dir, insertion_dict)
		print(insertion_locations_bed)

		#input variables
		target_regions = insertion_dict

	#region of interest mode    
	elif param_dictionary.get('mode') == "ROI":
		'''
		ROI mode needs a different pre-treatment of the reference genome, as ROIs are not spread randomly but have a fixed placement as defined in the bed file.
		'''
		print("Region of interest mode selected...")
		t0 = time.time()
		
		target_regions, bed, masked_regions, genome_size, chromosome_dir = create_barcoded_roi_genome(param_dictionary.get('reference_genome_path'), 
																					 param_dictionary.get('chr_restriction'), 
																					 param_dictionary.get('roi_bedpath'), 
																					 param_dictionary.get('n_barcodes'), 
																					 param_dictionary.get('blocked_regions_bedpath'))
	


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
		results_df = normalize_ROI_by_length(bed, results_df, scaling_factor=1000) #ROIPKB


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
