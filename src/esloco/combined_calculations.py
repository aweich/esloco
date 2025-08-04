import gc
import logging
import numpy as np

from esloco.config_handler import seq_read_data
from esloco.coverage import plot_reads_coverage
from esloco.counting import count_matches, count_barcode_occurrences
from esloco.read_operations import generate_read_length_distribution, generate_reads_based_on_coverage
from esloco.utils import profile_iteration, setup_logging, track_usage

def process_combination(
        mean_read_length,
        coverage,
        genome_size,
        target_regions,
        iteration,
        sequenced_data_path,
        min_read_length,
        n_barcodes,
        barcode_weights,
        masked_regions,
        output_path,
        scaling,
        min_overlap_for_detection,
        no_cov_plots,
        ):
    '''
    Creates a read length distribution based on mean read length and draws artificial reads from it until desired coverage is reached. 

    Then it compares the coordinates of the target regions (ROI or I) with the coordinates of artifical reads of each barcode and checks whether they are partially or fully contained.
    Masking is optional and is performed during the read generation step. 
    '''

    try:
        custom_read_length_distribution = seq_read_data(sequenced_data_path,
                                                        distribution=True,
                                                        min_read_length=min_read_length) #for sequencing data
        logging.info("Custom FASTA/FASTQ data provided.")

    except Exception:
        logging.info("No custom read length distribution provided... Generating artificial one...")
        custom_read_length_distribution = generate_read_length_distribution(num_reads=1000000,
                                                                            mean_read_length=mean_read_length,
                                                                            min_read_length = min_read_length,
                                                                            distribution='lognormal')
        logging.info(f'Calculating for: {mean_read_length}')

        if iteration == 1:
            logging.info("Saving generated read length distribution to numpy file...")
            np.save(f"{output_path}/{coverage}_{mean_read_length}_dsitribution.npy", custom_read_length_distribution)

    reads_to_simulate = int(coverage * genome_size / mean_read_length)
    precomputed_lengths = np.random.choice(custom_read_length_distribution, size=reads_to_simulate*10) #x10 to avoid running out of reads #check if this is faster?
    custom_cov_coordinates, covered_length = generate_reads_based_on_coverage(genome_size,
                                                                              coverage,
                                                                              precomputed_lengths,
                                                                              n_barcodes,
                                                                              barcode_weights,
                                                                              masked_regions)

    #plot coverage only for first iteration of each parameter combination
    if iteration < 1 and not no_cov_plots:
        plot_reads_coverage(genome_size,
                            1000000,
                            custom_cov_coordinates,
                            mean_read_length,
                            coverage,
                            target_regions,
                            output_path)

    # Sanity check for barcoding
    barcode_distribution = count_barcode_occurrences(custom_cov_coordinates)
    barcode_distribution["coverage"] = coverage
    barcode_distribution["mean_read_length"] = mean_read_length
    barcode_distribution["iteration"] = iteration
    barcode_distribution["N_bases"] = covered_length

    #Detections per ROI/Insertion
    detected = count_matches(target_regions, custom_cov_coordinates, scaling, min_overlap_for_detection)
    suffix = "_" + str(iteration)
    for insertion_data in detected:
        insertion_data["mean_read_length"] = mean_read_length
        insertion_data["coverage"] = coverage
        insertion_data["target_region"] += suffix

    #save memory from overload
    del custom_read_length_distribution
    del custom_cov_coordinates
    gc.collect()
    return detected, barcode_distribution

@profile_iteration
def run_simulation_iteration(iteration_id, param_dictionary, genome_size, target_regions, masked_regions, log_file):
    """
    Executes a single iteration of the simulation.
    """
    setup_logging(log_file)
    results, barcode_distributions = [], []
    for mean_read_length, coverage in param_dictionary['combinations']:
        out, barcode_distribution = process_combination(
            mean_read_length,
            coverage,
            genome_size,
            target_regions,
            iteration_id,
            param_dictionary['sequenced_data_path'],
            param_dictionary['min_read_length'],
            param_dictionary['n_barcodes'],
            param_dictionary['barcode_weights'],
            masked_regions,
            param_dictionary['output_path_plots'],
            param_dictionary['scaling'],
            param_dictionary['min_overlap_for_detection'],
            param_dictionary['no_cov_plots'],
        )
        results.append(out)
        barcode_distributions.append(barcode_distribution)
    track_usage(f"Iteration_{iteration_id}")
    return results, barcode_distributions
