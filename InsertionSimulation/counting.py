import logging
import random
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
            start = values[0]
            end = values[1]
        else:
            start = values["start"]
            end = values["end"]
        countdata = {'Insertion': key}
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
                                if random.random() <= scaling:
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
                            if random.random() <= scaling:
                                bias_list.append(bias)
                                partial_count += 1
                                overlaps.append(overlap)


        countdata['full_matches'] = full_length_count
        countdata['partial_matches'] = partial_count
        countdata['overlap'] = overlaps
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

def normalize_ROI_by_length(roi_input_bed, roi_counted_insertions, scaling_factor=1):
    '''
    Normalizes found ROI matches by their length and normalizes by number of reads
    '''
    results_full = []
    results_partial = []
    results_combined=[]
    lengths=[]
    for _, row_counted in roi_counted_insertions.iterrows():
        # Initialize length to None
        length = None
        
        # Iterate through each row in df1 to find partial match
        for _, row_initial in roi_input_bed.iterrows():
            if row_initial.id in row_counted.Insertion:
                # Calculate the length
                length = row_initial.end - row_initial.start
                lengths.append(length)
                break
        
        if length is not None:
            # Calculate the metric
            metric_full = ((row_counted['full_matches'] / length) * scaling_factor) / row_counted.Total_Reads #Reads per base (sf=1, kb sf=1000, or mb)
            results_full.append(metric_full)

            metric_partial = ((row_counted['partial_matches'] / length) * scaling_factor) / row_counted.Total_Reads #Reads per base (sf=1, kb sf=1000, or mb)
            results_partial.append(metric_partial)

            #combined ratio
            metric_combined = (((row_counted['partial_matches'] + row_counted['full_matches']) / length) * scaling_factor) / row_counted.Total_Reads #Reads per base (sf=1, kb sf=1000, or mb)
            results_combined.append(metric_combined)

        else:
            # If no partial match found, append None
            results_full.append(0)
            results_partial.append(0)
    
    # Add the results as a new column to DF2
    roi_counted_insertions['Length'] = lengths
    roi_counted_insertions['Full_Ratio'] = results_full
    roi_counted_insertions['Partial_Ratio'] = results_partial
    roi_counted_insertions['Combined_Ratio'] = results_combined
    
    return roi_counted_insertions