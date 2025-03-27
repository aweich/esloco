# config_handler.py
#config
import logging
import configparser
import ast
import json
import itertools
import sys

#mean read length calculation
import numpy as np
from Bio import SeqIO
import gzip

def mean_read_length(fasta_file):
    """
    Efficient mean read length calculation based on path to FASTA file.
    Supports both plain and gzipped FASTA files.
    """
    # Open the file, handling gzipped files if necessary
    open_func = gzip.open if fasta_file.endswith(".gz") else open
    with open_func(fasta_file, "rt") as handle:
        # Extract lengths and convert them to a numpy array
        lengths = np.array([len(record.seq) for record in SeqIO.parse(handle, "fasta")])
    return np.mean(lengths) if lengths.size > 0 else 0

def parse_config(config_file):
    """
    Parses the given configuration file and returns a dictionary of parameters.
    
    Args:
        config_file: Path to the configuration file.
    
    Returns:
        dict: A dictionary containing all parsed configuration parameters.
    """
    try:
        config = configparser.ConfigParser()
        config.read(config_file)

        # Required parameters (must be present in config)
        mode = config.get("COMMON", "mode", fallback=None)
        reference_genome_path = config.get("COMMON", "reference_genome_path", fallback=None)

        # Validate required parameters
        if mode is None:
            logging.error("Missing required parameter: mode")
            raise ValueError("Missing required parameter: mode")
        
        if reference_genome_path is None:
            logging.error("Missing required parameter: reference_genome_path")
            raise ValueError("Missing required parameter: reference_genome_path")

        # Optional parameters with default values
        sequenced_data_path = config.get("COMMON", "sequenced_data_path", fallback=None)
        output_path = config.get("COMMON", "output_path", fallback="./output/")
        experiment_name = config.get("COMMON", "experiment_name", fallback="default_experiment")
        output_path_plots = config.get("COMMON", "output_path_plots", fallback=output_path)
        min_overlap_for_detection = config.getint("COMMON", "min_overlap_for_detection", fallback=1)
        chr_restriction = config.get("COMMON", "chr_restriction", fallback=None)
        barcode_weights = config.get("COMMON", "barcode_weights", fallback="None")
        barcode_weights = ast.literal_eval(barcode_weights) if barcode_weights != "None" else None
        n_barcodes = config.getint("COMMON", "n_barcodes", fallback=1)
        iterations = config.getint("COMMON", "iterations", fallback=1)
        scaling = config.getfloat("COMMON", "scaling", fallback=1.0)
        parallel_jobs = config.getint("COMMON", "parallel_jobs", fallback=1)

        # Handle coverages and mean_read_lengths safely
        coverages = json.loads(config.get("COMMON", "coverages", fallback="[1]"))
        if not isinstance(coverages, list):
            coverages = [coverages]

        mean_read_lengths = json.loads(config.get("COMMON", "mean_read_lengths", fallback="[1000]"))
        if not isinstance(mean_read_lengths, list):
            mean_read_lengths = [mean_read_lengths]

        if sequenced_data_path:
            logging.info("Sequencing data provided, calculating mean read length...")
            mrl = int(mean_read_length(sequenced_data_path))
            logging.info(f"Mean read length set to: {mrl}")
            mean_read_lengths = [mrl]

        blocked_regions_bedpath = config.get("COMMON", "blocked_regions_bedpath", fallback=None)

        # Mode-specific parameters (handled safely)
        if mode == "ROI":
            roi_bedpath = config.get("ROI", "roi_bedpath", fallback=None)
        elif mode == "I":
            insertion_length = config.getint("I", "insertion_length", fallback=1000)
            insertion_number_distribution = config.get("I", "insertion_number_distribution", fallback="poisson")
            bedpath = config.get("I", "bedpath", fallback=None)
            insertion_numbers = config.getint("I", "insertion_numbers", fallback=5)
        else:
            logging.error("Invalid mode selected. Exiting.")
            raise ValueError("Invalid mode selected. Allowed values: 'ROI' or 'I'.")

        # Generate combinations of mean_read_lengths and coverages
        combinations = list(itertools.product(mean_read_lengths, coverages))

        # Store all configuration parameters in a dictionary
        param_dictionary = {key: value for key, value in locals().items() if key != "config" and not key.startswith("_")}
        return param_dictionary

    except Exception as e:
        logging.error(f"Configuration parsing failed: {e}")
        sys.exit(1)