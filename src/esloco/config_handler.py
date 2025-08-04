# config_handler.py
#config
import logging
import configparser
import ast
import json
import itertools
import sys
import os

#mean read length calculation
import gzip
import numpy as np
from Bio import SeqIO

def seq_read_data(fasta_file, distribution=False, min_read_length=0):
    """
    Efficient mean read length calculation based on path to FASTA or FASTQ file.
    Supports both plain and gzipped files.
    """
    # Open the file, handling gzipped files if necessary
    open_func = gzip.open if fasta_file.endswith(".gz") else open
    with open_func(fasta_file, "rt") as handle:

        # Determine file format based on extension
        fastq_exts = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
        fasta_exts = (".fasta", ".fasta.gz", ".fa", ".fa.gz")

        if fasta_file.endswith(fastq_exts):
            file_format = "fastq"
        elif fasta_file.endswith(fasta_exts):
            file_format = "fasta"
        else:
            raise ValueError(f"Unsupported file format for {fasta_file}. Only FASTA and FASTQ are allowed.")

        # Extract lengths and convert them to a numpy array
        lengths = []
        for record in SeqIO.parse(handle, file_format):
            if file_format == "fastq":
                # Ensure quality scores are present for FASTQ
                if not record.letter_annotations.get("phred_quality"):
                    raise ValueError(f"This is not a correct FASTQ file {fasta_file}. Provide FASTA or FASTQ.")
            # Filter sequences based on minimum read length
            if len(record.seq) > min_read_length:
                lengths.append(len(record.seq))
        lengths = np.array(lengths)
    if lengths.size == 0:
        raise ValueError(f"No reads found in {fasta_file}")
    if distribution:
        return lengths
    return np.mean(lengths)

def random_seed():
    """
    Generates a random seed if not provided.
    """
    return np.random.randint(0, 2**32 - 1)

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
            print("Missing required parameter: mode")
            raise ValueError("Missing required parameter: mode")

        if reference_genome_path is None:
            print("Missing required parameter: reference_genome_path")
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
        min_read_length = config.getint("COMMON", "min_read_length", fallback=1)
        no_cov_plots = config.getboolean("COMMON", "no_cov_plots", fallback=False)
        parallel_jobs = config.getint("COMMON", "parallel_jobs", fallback=1)
        seed = config.getint("COMMON", "seed", fallback=random_seed())

        #make sure output paths exist
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        if not os.path.exists(output_path_plots):
            os.makedirs(output_path_plots)

        # Handle coverages and mean_read_lengths safely
        coverages = json.loads(config.get("COMMON", "coverages", fallback="[1]"))
        if not isinstance(coverages, list):
            coverages = [coverages]

        mean_read_lengths = json.loads(config.get("COMMON", "mean_read_lengths", fallback="[1000]"))

        if not isinstance(mean_read_lengths, list):
            mean_read_lengths = [mean_read_lengths]

        if sequenced_data_path:
            print("Sequencing data provided. Calculating mean read length might take a while...")
            logging.info("Sequencing data provided, calculating mean read length...")
            mrl = int(seq_read_data(sequenced_data_path, min_read_length=min_read_length))
            logging.info("Mean read length set to: {mrl}")
            mean_read_lengths = [mrl]

        blocked_regions_bedpath = config.get("COMMON", "blocked_regions_bedpath", fallback=None)

        # Mode-specific parameters (handled safely)
        if mode == "ROI":
            roi_bedpath = config.get("ROI", "roi_bedpath", fallback=None)
        elif mode == "I":
            insertion_length = config.getint("I", "insertion_length", fallback=1000)
            insertion_number_distribution = config.get("I", "insertion_number_distribution", fallback=None)
            bedpath = config.get("I", "bedpath", fallback=None)
            insertion_numbers = config.getfloat("I", "insertion_numbers", fallback=5.0)
        else:
            print("Invalid mode selected. Exiting.")
            raise ValueError("Invalid mode selected. Allowed values: 'ROI' or 'I'.")

        # Generate combinations of mean_read_lengths and coverages
        combinations = list(itertools.product(mean_read_lengths, coverages))

        # Store all configuration parameters in a dictionary
        param_dictionary = {key: value for key, value in locals().items() if key != "config" and not key.startswith("_")}
        return param_dictionary

    except Exception as e:
        print("Configuration parsing failed: %s.", e)
        sys.exit(1)
