import logging
import numpy as np
import sys
from utils import track_usage

def calc_bias(read_length, target_length, min_overlap):    
    '''
    Calculates the bias correction factor based on the read length and target length.
    '''
    potential_start_full = read_length - target_length
    potential_start_partial = target_length - 1 + read_length - 1 - min_overlap - potential_start_full

    # Calculate the bias correction factor
    bias = potential_start_full / potential_start_partial
    return bias


def count_matches(insertion_dict, read_dir, scaling, min_overlap):
    '''
    Counts the number of full-length and partial insertions for each insertion/roi.
    Requires at least 'min_overlap' between the insertion and read intervals
    Genome scale factor due to diploidy of the human genome.
    '''
    data = []
    logging.info("Counting...")
    
    if not isinstance(scaling, (int, float)):
        scaling = 1 

    for key, values in insertion_dict.items():
        full_length_count = 0
        partial_count = 0
        
        #necessary since Insertions and ROIs have a different structure
        if type(values) == list: 
            start = int(values[0])
            end = int(values[1])
        else:
            try:
                start = int(values["start"])
                end = int(values["end"])
            except ValueError as e:
                logging.error(f"Invalid start or end value in key {key}: {values}. Error: {e}")
                continue
        countdata = {'target_region': key}
        overlaps=[]
        bias_list=[]
        for read, read_positions in read_dir.items():
            if read.split("_")[-1] == key.split("_")[1]:  # If barcodes match
                read_start, read_end = read_positions

                # Check if there is overlap first
                if max(start, read_start) < min(end, read_end):
                    
                    # Calculate overlap between the insertion and the read
                    overlap = min(end, read_end) - max(start, read_start)
                    if overlap >= min_overlap:

                        #calculate bias correction factor
                        read_length = read_end - read_start
                        target_length = end - start   
                        bias = calc_bias(read_length, target_length, min_overlap)
                        # Full-length insertion: check if the read fully covers the insertion
                        if read_start <= start and end <= read_end:
                                if np.random.rand() <= scaling:
                                    bias_list.append(bias)
                                    full_length_count += 1
                                    overlaps.append(overlap)


                        # Partial insertion: check if the read partially overlaps with the insertion
                        #---[-- ]
                        #   [ - ]
                        #   [ --]---
                        #elif (start >= read_start and end > read_end) or \
                        #        ((start <= read_start and end >= read_end) and (start != read_start and end != read_end)) or \
                        #        (start < read_start and end <= read_end):
                        else:
                            # Ensure the overlap is at least 'x'
                            if np.random.rand() <= scaling:
                                bias_list.append(bias)
                                partial_count += 1
                                overlaps.append(overlap)


        countdata['full_matches'] = full_length_count
        countdata['partial_matches'] = partial_count
        countdata['overlap'] = overlaps
        countdata["on_target_bases"] = sum(eval(x) if isinstance(x, str) else x for x in countdata["overlap"])
        countdata['bias'] = bias_list
        data.append(countdata)

    track_usage("count_matches")
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
