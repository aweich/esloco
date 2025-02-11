#%%	
#!/usr/bin/env python3
from Bio import SeqIO
import os
import random
import numpy as np
import pandas as pd
import sys
import time
import logging
import itertools
from functools import partial
from joblib import Parallel, delayed
from tqdm import tqdm

# custom config reader
from config_handler import parse_config

#custom functions
#from create_insertion_genome import add_insertions_to_genome_sequence_with_bed
#from utils import profile, check_barcoding, roi_barcoding, barcode_genome
from bed_operations import global_to_chromosome_coordinates
from counting import normalize_ROI_by_length
from genome_generation import create_barcoded_insertion_genome, create_barcoded_roi_genome
from combined_calculations import run_simulation_iteration

def main():
    """ Main function to execute the entire simulation. """

    # Load configuration
    if len(sys.argv) != 2:
        print("Usage: python InsertionGenomeBuilder.py <config_file>")
        sys.exit(1)
    
    config_file = sys.argv[1]
    
    try:
        param_dictionary = parse_config(config_file)
    except Exception as e:
        print(f"Error parsing config: {e}")
        sys.exit(1)

    # Extract key parameters once
    mode = param_dictionary.get("mode")
    output_path = param_dictionary.get("output_path")
    experiment_name = param_dictionary.get("experiment_name")
    num_iterations = param_dictionary.get("iterations")
    parallel_jobs = param_dictionary.get("parallel_jobs")

	# Create output dir
    if not os.path.exists(output_path):
         os.makedirs(output_path)
    
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
    filename=os.path.join(output_path, 'tcr_simple_log.log'),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filemode='w'
	)
    
    logging.info(f"Running in mode: {mode}")

    start_time = time.time()
    genome_size, target_regions, masked_regions, chromosome_dir = None, None, None, None

    # Process based on mode
    if mode == "I":
        logging.info("Processing Insertion Mode...")
        genome_size, insertion_dict, masked_regions, chromosome_dir = create_barcoded_insertion_genome(
            param_dictionary['reference_genome_path'], param_dictionary['bedpath'],
            param_dictionary['blocked_regions_bedpath'], param_dictionary['chr_restriction'],
            param_dictionary['insertion_length'], param_dictionary['insertion_numbers'],
            param_dictionary['insertion_number_distribution'], param_dictionary['n_barcodes']
        )
        logging.info(f"Number of insertions: {len(insertion_dict)}")
        target_regions = insertion_dict
        insertion_locations_bed = global_to_chromosome_coordinates(chromosome_dir, insertion_dict)

    elif mode == "ROI":
        logging.info("Processing Region of Interest Mode...")
        target_regions, bed, masked_regions, genome_size, chromosome_dir = create_barcoded_roi_genome(
            param_dictionary['reference_genome_path'], param_dictionary['chr_restriction'],
            param_dictionary['roi_bedpath'], param_dictionary['n_barcodes'],
            param_dictionary['blocked_regions_bedpath']
        )

    else:
        logging.error("Error: Invalid mode selected.")
        sys.exit(1)

    # Parallel execution
    parallel_results = Parallel(n_jobs=parallel_jobs)(
        delayed(run_simulation_iteration)(i, param_dictionary, genome_size, target_regions, masked_regions)
        for i in tqdm(range(num_iterations))
    )

    # Unpack nested structure
    results_list, barcode_distributions_list = [], []
    for iteration_results, barcode_distributions in parallel_results:
        results_list.extend(iteration_results)
        barcode_distributions_list.append(barcode_distributions)

    # Convert to DataFrame
    results_df = pd.DataFrame([item for sublist in results_list for item in sublist])
    barcode_distributions_df = pd.DataFrame(itertools.chain(*barcode_distributions_list))

    # ROI-specific normalization
    if mode == "ROI":
        results_df['iteration'] = results_df['Insertion'].str.split('_').str[-1].astype(int)
        barcode_distributions_df['Total_Reads'] = barcode_distributions_df.iloc[:, :param_dictionary['n_barcodes']].sum(axis=1)
        results_df = results_df.merge(barcode_distributions_df, on=['coverage', 'mean_read_length', 'iteration'], how='inner')
        results_df = normalize_ROI_by_length(bed, results_df, scaling_factor=1000)

    # Save results
    try:
        if mode == "I":
            insertion_locations_bed.to_csv(f"{output_path}{experiment_name}_insertion_locations.bed",
                                           sep='\t', header=True, index=False)
            logging.info(f"Insertion locations saved.")

        barcode_distributions_df.to_csv(f"{output_path}{experiment_name}_barcode_distribution_table.csv",
                                        sep='\t', header=True, index=False)
        results_df.to_csv(f"{output_path}{experiment_name}_matches_table.csv",
                          sep='\t', header=True, index=False)

    except Exception as e:
        logging.error(f"Error writing output files: {e}")

    #  Configuration for reference
    logging.info("Simulation Parameters:")
    for param, value in param_dictionary.items():
        logging.info(f"{param} = {value} ({type(value)})")

    logging.info(f"Total execution time: {time.time() - start_time:.2f} seconds.")
    logging.info("Simulation complete!")



if __name__ == "__main__":
	try:
		main()
	except ValueError as e:
		print("Configuration error:", e)
		sys.exit(1)

# %%
