import logging
from tqdm_joblib import ParallelPbar
from joblib import delayed
from tqdm import tqdm

from esloco.create_insertion_genome import add_insertions_to_genome_sequence_with_bed
from esloco.utils import check_barcoding, barcode_genome, track_usage
from esloco.fasta_operations import pseudo_fasta_coordinates
from esloco.bed_operations import readbed, chromosome_to_global_coordinates

def parallel_barcoded_insertion_genome(i,chromosome_dir, bedpath, n_barcodes, ref_genome_size, insertion_length, insertion_numbers, insertion_number_distribution):
    barcoded_chromosome_dir = barcode_genome(chromosome_dir, i)
    if not bedpath or bedpath.lower() == "none":
        logging.info("Insertions will be placed randomly...")
        bed_df = None
    else:
        logging.info(f"Guided insertion placement within regions defined in {bedpath}.")
        bed_df = readbed(bedpath, barcoded_chromosome_dir.keys(), barcoding=check_barcoding(n_barcodes))
    insertion_dict = add_insertions_to_genome_sequence_with_bed(
        ref_genome_size, insertion_length, insertion_numbers, barcoded_chromosome_dir, insertion_number_distribution, bed_df
    )
    return insertion_dict

def create_barcoded_insertion_genome(parallel_jobs, reference_genome_path, bedpath, blocked_regions_bedpath, restriction, insertion_length, insertion_numbers, insertion_number_distribution, n_barcodes):
    '''
    Pre-processment step of the insertion mode.
    The reference genome cooridnates are transformed into a string-like format (One single FASTA string) and the chromsome borders are stored
    Based on user input, masked regions are defined and weighted
    For each barcode, the chromosome borders get a prefix and inserttion positions are randomly chosen. 
    If an insertion-bed is defined, the insertions are only placed within these regions accoridng to their length-based weights.
    The length of the reference genome (for coverage estimation), the insertion cooridnates, the masked regions, and the chromosome border dict are returned
    '''
    logging.info("Create Barcoded Insertion Genome...")
    collected_insertion_dict = {}
    #1
    ref_genome_size, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path, restriction)
    #2 #create blocked regions file
    if not blocked_regions_bedpath or blocked_regions_bedpath.lower() == "none":
        logging.info("No regions provided for masking...")
        masked_regions = None
    else:
        logging.info(f"Masking regions as defined in {blocked_regions_bedpath}")
        blocked_bed = readbed(blocked_regions_bedpath, chromosome_dir.keys())
        masked_regions = chromosome_to_global_coordinates(blocked_bed, chromosome_dir)
    #3 Parallelized step
    parallel_results = ParallelPbar(f"Creating {n_barcodes} I Genome(s)...")(n_jobs=parallel_jobs)(
        delayed(parallel_barcoded_insertion_genome)(i,chromosome_dir, bedpath, n_barcodes, ref_genome_size, insertion_length, insertion_numbers, insertion_number_distribution)
        for i in range(n_barcodes)
    )
    #4 Unpack
    for insertion_dict in parallel_results:
        for key, value in insertion_dict.items():
            if key not in collected_insertion_dict:
                collected_insertion_dict[key] = []
            collected_insertion_dict[key].extend(value)
    track_usage("create_barcoded_insertion_genome")
    return ref_genome_size, collected_insertion_dict, masked_regions, chromosome_dir

def create_barcoded_roi_genome(reference_genome_path, restriction, roi_bedpath, n_barcodes, blocked_regions_bedpath):
    '''
    Pre-processment step of the roi mode.
    The reference genome cooridnates are transformed into a string-like format (One single FASTA string) and the chromsome borders are stored
    Based on user input, masked regions are defined and weighted
    The pre-defined ROIs are barcoded.
    The length of the reference genome (for coverage estimation), the insertion cooridnates, the masked regions, and the chromosome border dict are returned
    '''
    logging.info("Create Barcoded ROI Genome...")
    #create global cooridnates from bed based on provided genome ref to adjust ROIs to string-like genome
    genome_size, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path, restriction)
    bed = readbed(roi_bedpath, chromosome_dir.keys())
    roi_dict = chromosome_to_global_coordinates(bed, chromosome_dir)
    logging.info("Barcoding ROIs...")
    barcoded_roi_dict = {}
    for i in tqdm(range(n_barcodes),  desc=f"Creating {n_barcodes} ROI Genome(s)..."):
        for key, value in roi_dict.items():
            new_key = f"{key}_{i}"
            barcoded_roi_dict[new_key] = value
    #create blocked regions file
    if not blocked_regions_bedpath or blocked_regions_bedpath.lower() == "none":
        logging.info("No regions provided for masking...")
        masked_regions = None
    else:
        logging.info(f"Masking regions as defined in {blocked_regions_bedpath}.")
        blocked_bed = readbed(blocked_regions_bedpath, chromosome_dir.keys())
        masked_regions = chromosome_to_global_coordinates(blocked_bed, chromosome_dir)
    track_usage("create_barcoded_roi_genome")
    return barcoded_roi_dict, bed, masked_regions, genome_size, chromosome_dir
