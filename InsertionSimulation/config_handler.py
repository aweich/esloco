# config_handler.py
#config
import configparser
import ast
import json
import itertools
import sys

#mean read length calculation
import numpy as np
from Bio import SeqIO

def mean_read_length(fasta_file):
    """
    Efficient mean read length calculation based on path to FASTA file.
    """
    # Extract lengths and convert them to a numpy array
    lengths = np.array([len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")])
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

        # Common parameters
        mode = config['COMMON']['mode']
        reference_genome_path = config['COMMON']['reference_genome_path']
        sequenced_data_path = config['COMMON']['sequenced_data_path']
        output_path = config['COMMON']['output_path']
        experiment_name = config['COMMON']['experiment_name']
        output_path_plots = config['COMMON']['output_path_plots']
        min_overlap_for_detection = int(config['COMMON']['min_overlap_for_detection'])
        chr_restriction = config['COMMON']['chr_restriction']
        barcode_weights = ast.literal_eval(config['COMMON']['barcode_weights'])
        n_barcodes = int(config['COMMON']['n_barcodes'])
        iterations = int(config['COMMON']['iterations'])
        scaling = config['COMMON']['scaling']
        parallel_jobs = int(config['COMMON']['parallel_jobs'])

        coverages = json.loads(config.get("COMMON", "coverages"))
        if type(coverages) is not list:
            coverages = [coverages]

        mean_read_lengths = json.loads(config.get("COMMON", "mean_read_lengths"))
        if type(mean_read_lengths) is not list:
            mean_read_lengths = [mean_read_lengths]
        
        if sequenced_data_path:
            print("Sequencing data and a mean read length were provided...")
            print("Only sequencing data will be used for the read length distribution...")
            mrl = int(mean_read_length(sequenced_data_path))
            print(f"Mean read length is set to your data specific mean: {mrl}")
            mean_read_lengths = [mrl]
            del mrl

        blocked_regions_bedpath = config['COMMON']['blocked_regions_bedpath']

        # Mode-specific parameters
        if mode == "ROI":
            roi_bedpath = config['ROI']['roi_bedpath']
        elif mode == "I":
            insertion_length = int(config['I']['insertion_length'])
            insertion_number_distribution = config['I']['insertion_number_distribution']
            bedpath = config['I']['bedpath']
            insertion_numbers = int(config['I']['insertion_numbers'])
        else:
            print("No valid running mode selected.")
            sys.exit()

        # Generate combinations of mean_read_lengths and coverages
        combinations = list(itertools.product(mean_read_lengths, coverages))

        # Store all configuration parameters in a dictionary
        param_dictionary = {key: value for key, value in locals().items() if key != "config" and not key.startswith("_")}
        return param_dictionary

    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)