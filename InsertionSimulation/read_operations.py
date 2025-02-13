import random
from Bio import SeqIO
import numpy as np
import logging
import sys

from utils import track_usage

def get_read_length_distribution_from_real_data(path_to_sequenced_fasta):
    '''
    Uses an input fasta and generates a custom read length distribution.
    '''
    lengths=[]
    for record in SeqIO.parse(path_to_sequenced_fasta, 'fasta'):
            lengths.append(len(record.seq))
    return lengths

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
        logging.error("Unsupported distribution. Supported options: 'normal', 'lognormal'.")
        sys.exit(1)
    
    # Filter out zero-length reads
    read_lengths = np.round(read_lengths).astype(int)
    read_lengths = read_lengths[read_lengths > 0]
    
    return read_lengths #not returning list

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
					return (weights_dict[key] * n_barcodes) / common_denominator
			
			#No weight provided for this barcode, using the barcodes share of the common denominator 
			return n_barcodes / common_denominator
		else:
			return 1 / n_barcodes

def generate_reads_based_on_coverage(genome_size, coverage, precomputed_lengths, n_barcodes, barcode_weights, masked_regions=None):
    '''
    Randomly pulls a read of size X derived from the read length distribution from the fasta until the fasta is N times covered (coverage).
    '''
    logging.info("Coverage: " + str(coverage))
    logging.info("Pulling reads...")
    covered_length = 0
    read_coordinates = {}
    barcode_names = ["Barcode_" + str(i) for i in range(n_barcodes)]

    #Reads pulled until desired coverage is reached
    while covered_length < coverage * genome_size:
        random_barcode = random.choices(barcode_names, [get_weighted_probabilities(i, n_barcodes, barcode_weights) for i in barcode_names])[0] #chooses one of the barcodes based on weighted probability
        read_length = precomputed_lengths.pop()
        start_position = random.randint(0, genome_size - read_length)
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
    track_usage("generate_reads_based_on_coverage")
    return read_coordinates, covered_length

