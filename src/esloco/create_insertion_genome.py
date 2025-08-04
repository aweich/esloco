# create_insertion_genome.py
import logging
import numpy as np

#custom
from esloco.utils import get_chromosome

def add_insertions_to_genome_sequence_with_bed(reference_sequence,
                                               insertion_length,
                                               num_insertions,
                                               chromosome_dir,
                                               insertion_number_distribution=None,
                                               bed_df=None,
                                               ):
    '''
    Randomly add insertion sequence into the reference genome or within specified regions.

    For each insertion, a random position on the reference sequence is chosen and its coordinates are stored. 
    For all follwoing insertions, the coordinates of their location are updated if needed (i.e. if they happend to insert at an earlier position)
    In case of a bed-guided insertion, each region in the file is assigned a probability according to its length.
    '''
    position = {}
    if insertion_number_distribution == 'poisson':
        num_insertions = np.random.poisson(num_insertions)
        logging.info(f"Number of insertions drawn from Poisson distribution: {num_insertions}")
    else:
        logging.info(f"Using exactly {num_insertions}.")
    num_insertions = int(num_insertions)
    if bed_df is not None:
        logging.info("BED guided insertion pattern...")
        # Step 1: Calculate probabilities based on region lengths
        logging.info("Calculating insertion probabilities (region length / sum of all regions lengths)...")
        region_lengths = bed_df['end'] - bed_df['start']
        region_probabilities = region_lengths / region_lengths.sum()
        for i in range(num_insertions):
            # Step 2: Randomly select insertion regions #so that each region is selected once!
            if len(bed_df.index) == num_insertions:
                logging.info("Fixed insertion locations selected as provided by the BED.")
                selected_region = bed_df.iloc[i]
            else:
                selected_region_index = np.random.choice(bed_df.index, p=region_probabilities)
                selected_region = bed_df.iloc[selected_region_index]
            # Step 3: Perform insertions within selected regions
            chromosome = selected_region['chrom']
            chromosome_range = chromosome_dir[chromosome]
            insert_position = np.random.randint(selected_region['start'], selected_region['end'])
            # Adjust insertion position to the global genomic coordinates
            global_insert_position = chromosome_range[0] + insert_position
            # Update position table
            for key, value in position.items():
            # If the insertion is after the current position, update the position
                if value[0] >= global_insert_position:
                    position[key][0] += insertion_length
                    position[key][1] += insertion_length
            # Add the new insertion position and (barcoded) name
            insertion_name = chromosome.split('chr')[0] + f"insertion_{i}"
            position[insertion_name] = [global_insert_position, global_insert_position + insertion_length]
    #if no bed is provided
    else:
        for i in range(num_insertions):
            # Choose a random position to insert the smaller sequence
            insert_position = np.random.randint(0, reference_sequence)
            #check in which chr it landed
            chromosome = get_chromosome(insert_position, chromosome_dir)
            # Update position table
            for key, value in position.items():
            # If the insertion is after the current position, update the position
                if value[0] >= insert_position:
                    position[key][0] += insertion_length
                    position[key][1] += insertion_length
            # Add the new insertion position and (barcoded) name
            insertion_name = chromosome.split('chr')[0] + f"insertion_{i}"
            position[insertion_name]= [insert_position, insert_position + insertion_length]
    #new genome length is original length + length of all insertions combined
    #length_updated_reference_sequence = reference_sequence + num_insertions * insertion_length
    #length_updated_reference_sequence,
    return position
