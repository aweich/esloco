#!/usr/bin/env python3
#%%
import sys
import time
import logging
import itertools
import os
import numpy as np
import pandas as pd
from joblib import delayed
from tqdm_joblib import ParallelPbar

# custom config reader
from esloco.config_handler import parse_config

#custom functions
#from create_insertion_genome import add_insertions_to_genome_sequence_with_bed
from esloco.utils import track_usage, setup_logging
from esloco.bed_operations import global_to_chromosome_coordinates
from esloco.genome_generation import create_barcoded_insertion_genome, create_barcoded_roi_genome
from esloco.combined_calculations import run_simulation_iteration

def print_help():
    print("")
    print("Usage: ")
    print("         esloco --config <config_file>")
    print("")
    print("Run a simulation for local coverage estimation.")
    print("")
    print("Options:")
    print("     --config  <config_file>    Path to the configuration file.")
    print("     --help                     Show this help message.")
    print("")
    sys.exit(0)

def main():
    """ Main body to execute the entire simulation. """

    # Load configuration
    if len(sys.argv) != 3 or sys.argv[1] == "--help" or sys.argv[1] != "--config":
        print_help()
        sys.exit(1)

    if sys.argv[1] == "--config":
        config_file = sys.argv[2]

    if not os.path.exists(config_file):
        print(f"Configuration file '{config_file}' does not exist.")
        sys.exit(1)

    # immediate CLI feedback after starting the simulation
    print(f"Checking config...")

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

    # Further CLI startup feedback
    print(f"Starting simulation: {experiment_name}")
    print(f"Running in mode: {mode}")

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Setup logfile
    log_file = os.path.join(output_path, f"{experiment_name}_log.log")

    # Check if log file already exists
    if os.path.exists(log_file):
        logging.warning(f"Log file already exists. Continuing simulation. New info will be appended to '{log_file}'.")

    setup_logging(log_file)

    # prints startup message to log after initializing of the log file
    logging.info(f"Starting simulation: {experiment_name}")
    logging.info(f"Running in mode: {mode}")

    # Setup monitoring
    start_time = time.time()
    track_usage("Before")

    # Assign Seed
    np.random.seed(param_dictionary['seed'])
    logging.info(f"Seed for reproducibility: {param_dictionary['seed']}")
    genome_size, target_regions, masked_regions, chromosome_dir = None, None, None, None

    # Process based on mode
    if mode == "I":
        logging.info("Processing Insertion Mode...")
        genome_size, insertion_dict, masked_regions, chromosome_dir = create_barcoded_insertion_genome(
            parallel_jobs,
            param_dictionary['reference_genome_path'],
            param_dictionary['bedpath'],
            param_dictionary['blocked_regions_bedpath'],
            param_dictionary['chr_restriction'],
            param_dictionary['insertion_length'],
            param_dictionary['insertion_numbers'],
            param_dictionary['insertion_number_distribution'],
            param_dictionary['n_barcodes'],
        )
        logging.info(f"Number of insertions: {len(insertion_dict)}")
        target_regions = insertion_dict
        insertion_locations_bed = global_to_chromosome_coordinates(chromosome_dir, insertion_dict)
    elif mode == "ROI":
        logging.info("Processing Region of Interest Mode...")
        target_regions, bed, masked_regions, genome_size, chromosome_dir = create_barcoded_roi_genome(
            param_dictionary['reference_genome_path'],
            param_dictionary['chr_restriction'],
            param_dictionary['roi_bedpath'],
            param_dictionary['n_barcodes'],
            param_dictionary['blocked_regions_bedpath'])
        logging.info(f"Chromosome borders: {chromosome_dir}")
        logging.info(f"Target regions: {target_regions}")
        logging.info(f"Original Target region coordinates: {bed}")
    else:
        logging.error("Error: Invalid mode selected.")
        sys.exit(1)

    # Parallel execution
    parallel_results = ParallelPbar("Iterations...")(n_jobs=parallel_jobs)(
        delayed(run_simulation_iteration)(i, param_dictionary, genome_size, target_regions, masked_regions, log_file)
        for i in range(num_iterations)
    )

    # Unpack nested structure
    results_list, barcode_distributions_list = [], []
    for iteration_results, barcode_distributions in parallel_results:
        results_list.extend(iteration_results)
        barcode_distributions_list.append(barcode_distributions)

    # Convert to DataFrame
    results_df = pd.DataFrame([item for sublist in results_list for item in sublist])
    barcode_distributions_df = pd.DataFrame(itertools.chain(*barcode_distributions_list))

    # output formatting
    if mode == "ROI":
        # works with {ROI}_{barcode-number}_{iteration-number}
        results_df[["target", "barcode", "iteration"]] = results_df["target_region"].str.rsplit("_", n=2, expand=True)
        barcode_distributions_df['total_Reads'] = barcode_distributions_df.iloc[:, :param_dictionary['n_barcodes']].sum(axis=1)

    else:
        # works with Barcode_{barcode-number}_insertion_{insertion-number}_{iteration-number}
        # Split the target_region column into its components
        target_split = results_df["target_region"].str.rsplit("_", n=4, expand=True)
        results_df["target"] = target_split.iloc[:, :4].agg("_".join, axis=1)
        results_df["barcode"] = target_split.iloc[:, 1]
        results_df["insertion"] = target_split.iloc[:, 2] + "_" + target_split.iloc[:, 3]
        results_df["iteration"] = target_split.iloc[:, 4]

        # Calculate total reads for barcode distributions
        barcode_distributions_df['total_Reads'] = barcode_distributions_df.iloc[:, :param_dictionary['n_barcodes']].sum(axis=1)

    # simple results
    simple_results = results_df.groupby(["target", "mean_read_length", "coverage"], as_index=False).agg(
        mean_full_matches=('full_matches', 'mean'),
        sd_full_matches=('full_matches', 'std'),
        mean_partial_matches=('partial_matches', 'mean'),
        sd_partial_matches=('partial_matches', 'std'),
        mean_bases_on_target=('on_target_bases', 'mean'),
        sd_bases_on_target=('on_target_bases', 'std')
    )

    # Save results
    try:
        if mode == "I":
            insertion_locations_bed.to_csv(f"{output_path}{experiment_name}_insertion_locations.bed",
                                           sep='\t',
                                           header=True,
                                           index=False)
            logging.info(f"Insertion locations saved.")
        barcode_distributions_df.to_csv(f"{output_path}{experiment_name}_barcode_distribution_table.csv",
                                        sep='\t',
                                        header=True,
                                        index=False)
        results_df.to_csv(f"{output_path}{experiment_name}_matches_table.csv",
                          sep='\t',
                          header=True,
                          index=False)
        simple_results.to_csv(f"{output_path}{experiment_name}_summary.csv",
                              sep='\t',
                              header=True,
                              index=False)

    except Exception as e:
        logging.error(f"Error writing output files: {e}")

    #  Configuration for reference
    logging.info("Simulation Parameters:")
    for param, value in param_dictionary.items():
        logging.info(f"{param} = {value} ({type(value)})")

    track_usage("After")
    logging.info(f"Total execution time: {time.time() - start_time:.2f} seconds.")
    logging.info("Simulation complete!")
    print("Simulation complete!")
    sys.exit()

if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        print("Configuration error:", e)
        sys.exit(1)

# %%
