#!/usr/bin/env python3

from Bio import SeqIO
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import time
import gc
import itertools
#import multiprocessing
from functools import partial
from joblib import Parallel, delayed

#will at some point be changes into argparse or sth similar
reference_genome_path = "/home/weichan/permanent/Projects/VIS/dev/VIS_Magdeburg_withBasecalling/hg38.fa" #reads will be created based on this reference
vector_sequence_path = "/home/weichan/permanent/Projects/VIS/dev/VIS_Magdeburg_withBasecalling/pSLCAR-CD19-28z.fasta"#vector #currently 8866 - 42000 (not observed in data): 5kb long should be enough!
sequenced_data_path = "/home/weichan/permanent/Projects/VIS/VIS_integration_site/Results/FullRunAfterModulaization_BUFFERMODE100_CD19_cd247_Vector_integration_site/FASTA/Full_MK025_GFP+.fa"
output_path = "./out/Test_100xIntrons_SummaryTable.csv"
bedpath = "/home/weichan/permanent/Projects/VIS/dev/UCSC/UCSC_GENCODE_V44_Introns_04_24"
coverage=5
insertion_numbers=5
#mean_read_length=5000 #for artificial reads
#Part 1: Create an insertion-infiltrated chromosome and check where the insertions happen.
def timeit(func):
    '''
    Timer decorater to time all the functions. 
    '''
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        elapsed = end - start
        print(f'Time taken: {elapsed:.6f} seconds')
        return result
    return wrapper

def readbed(bedpath, list_of_chromosomes_in_reference):
    bed = pd.read_csv(bedpath, sep='\t', header=None, usecols=[0,1,2], names=['chrom', 'start', 'end'])
    print('Only keeping the following chromosomes...')
    print(list_of_chromosomes_in_reference)
    bed = bed[bed["chrom"].isin(list_of_chromosomes_in_reference)]
    return bed

@timeit
def collapse_fasta(path_to_fasta):
    with open(path_to_fasta, 'r') as fasta_file:
        # extracting multiple data in single fasta file using biopython
        seqList=[]
        for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
            seqList.append(str(record.seq))
    return ''.join(seqList)

@timeit
def pseudo_fasta_coordinates(path_to_fasta):
    '''
    Takes in FASTA with separate entries and returns a table with the pseudo entry borders of each entry after the collapsing of the entries.
    '''
    with open(path_to_fasta, 'r') as fasta_file:
        entries={}
        seqList=[]
        updated_length=0
        for record in SeqIO.parse(fasta_file, 'fasta'):
        	if '_' not in record.id and 'M' not in record.id:# and '2' in record.id: #the 2 modification will be removed later but is good for assessing the performance!
        	    entries[record.id] = [updated_length, updated_length + len(record.seq)]
        	    updated_length = updated_length + len(record.seq)
        	    seqList.append(str(record.seq))
    return ''.join(seqList), entries

@timeit #this step takes long and is dependent on the number of insertions since for each insertion we have to overwrite it in memory! O(n) for chr12,2,20,21,22: 20 insertions take 13 seconds
def add_insertions_to_genome_sequence(reference_sequence, insertion_sequence, num_insertions):
    """Randomly add insertion sequence into the reference sequence."""
    position = {}

    for i in range(num_insertions):
        # Choose a random position to insert the smaller sequence
        insert_position = random.randint(0, len(reference_sequence))
        # Insert the insertion sequence at the chosen position
        reference_sequence = (
            reference_sequence[:insert_position] +
            insertion_sequence +
            reference_sequence[insert_position:]
        )
        # Update position table
        for key, value in position.items():
        # If the insertion is after the current position, update the position
            if value[0] >= global_insert_position:
                position[key][0] += len(insertion_sequence)
                position[key][1] += len(insertion_sequence)

        # Update position table
        position["insertion_%s" % i]= [insert_position, insert_position + len(insertion_sequence)]
    return reference_sequence, position


@timeit
def insertions_per_chromosome(chromosome_dir, insertion_dir):
	"""
	Creates a table that lists the number of insertions per chromosome.
	"""
	# Convert chromosome_dir and insertion_dir to sets for faster lookups
	chromosome_set = {chromosome: tuple(coordinates) for chromosome, coordinates in chromosome_dir.items()}
	insertion_set = {insertion: tuple(coordinates) for insertion, coordinates in insertion_dir.items()}
	chromosome_lengths = {}
	# Iterate through chromosome_dir to get chromosome lengths
	for chromosome, coordinates in chromosome_set.items():
	    chromosome_lengths[chromosome] = coordinates[1] - coordinates[0]

	# Initialize a dictionary to count insertions per chromosome
	insertions_per_chromosome = {chromosome: 0 for chromosome in chromosome_dir}

	# Iterate through insertions
	for insertion, insertion_coordinates in insertion_set.items():
	    # Check if the insertion falls within any chromosome
	    for chromosome, chromosome_coordinates in chromosome_set.items():
	        if chromosome_coordinates[0] <= insertion_coordinates[0] < chromosome_coordinates[1] or \
	           chromosome_coordinates[0] < insertion_coordinates[1] <= chromosome_coordinates[1]:
	            insertions_per_chromosome[chromosome] += 1

	# Convert the dictionary to a DataFrame
	insertions_df = pd.DataFrame.from_dict(insertions_per_chromosome, orient='index', columns=['Num_Insertions'])
	insertions_df.index.name = 'Chromosome'
	# Add a column for chromosome lengths
	insertions_df['Chromosome_Length'] = insertions_df.index.map(chromosome_lengths)

	return insertions_df

#new code for bed consideration option

#This is the basic functionality that needs to be implemented to only have bed file specific insertions based on the weighted coordinates
def add_insertions_to_genome_sequence_with_bed(reference_sequence, insertion_sequence, num_insertions, chromosome_dir=None, bed_df=None):
    """Randomly add insertion sequence into the reference genome within specified regions."""
    position = {}
    
    if bed_df is not None:
        # Step 1: Calculate probabilities based on region lengths
        print("Calculating insertion probabilities (region length / sum of all regions lengths)...")
        region_lengths = bed_df['end'] - bed_df['start']
        region_probabilities = region_lengths / region_lengths.sum()
        updated_reference_sequence = reference_sequence

        for i in range(num_insertions):
            # Step 2: Randomly select insertion regions
            selected_region_index = np.random.choice(bed_df.index, p=region_probabilities)
            selected_region = bed_df.iloc[selected_region_index]

            # Step 3: Perform insertions within selected regions
            chromosome = selected_region['chrom']
            chromosome_range = chromosome_dir[chromosome]

            insert_position = random.randint(selected_region['start'], selected_region['end'])
            # Adjust insertion position to the global genomic coordinates
            global_insert_position = chromosome_range[0] + insert_position

            # Insert the insertion sequence into the reference sequence
            updated_reference_sequence = (
                updated_reference_sequence[:global_insert_position] +
                insertion_sequence +
                updated_reference_sequence[global_insert_position:]
            )

            # Update position table
            for key, value in position.items():
            # If the insertion is after the current position, update the position
                if value[0] >= global_insert_position:
                    position[key][0] += len(insertion_sequence)
                    position[key][1] += len(insertion_sequence)

            # Add the new insertion position
            position["insertion_%s" % i] = [global_insert_position, global_insert_position + len(insertion_sequence)]


        return updated_reference_sequence, position
    #if no bed is provided
    updated_reference_sequence, position = add_insertions_to_genome_sequence(reference_sequence=reference_sequence, insertion_sequence=insertion_sequence, num_insertions=num_insertions)
    return updated_reference_sequence, position


#1 Creates the vector in the right format
#insertion_fasta = collapse_fasta(vector_sequence_path)
#insertion_fasta = 'X' * 5000
#2 Creates a single-string reference genome, while storing the chromosome borders
#fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
#print(chromosome_dir)
#bed_df = readbed(bedpath, chromosome_dir.keys())
#3 Randomly inserts insertio sequence into the single-string reference x times
#mod_fasta, insertion_dir = add_insertions_to_genome_sequence_with_bed(fasta, insertion_fasta, insertion_numbers, chromosome_dir, bed_df)
#4 Creates an overview table of the insertion distribution
#insertions_df = insertions_per_chromosome(chromosome_dir, insertion_dir)
#print(insertions_df)

#sys.exit()
#Part 2: Create artifical reads
@timeit
def get_read_length_distribution_from_real_data(path_to_sequenced_fasta):
	lengths=[]
	for record in SeqIO.parse(path_to_sequenced_fasta, 'fasta'):
	        lengths.append(len(record.seq))
	return lengths

@timeit
def generate_reads(fasta, read_length_distribution, num_reads):
    """Generate reads based on custom read length distribution and pulled from fasta ref."""
    reads = []
    for read_length in read_length_distribution:
        start_position = random.randint(0, len(fasta) - read_length)
        read = fasta[start_position:start_position + read_length]
        reads.append(read)
    return reads

@timeit
def generate_reads_based_on_coverage(fasta, read_length_distribution, coverage, PRECOMPUTE_RANDOM):
    '''
    Randomly pulls a read of size X derived from the read length distribution from the fasta until the fasta is N times covered (coverage). 
    '''
    print("Coverage: " + str(coverage))
    print("Pulling reads...")
    total_length = len(fasta)
    covered_length = 0
    #reads = [] #only a goodf idea if we are not testing high coverages, otherwise memory is floated
    read_coordinates = {}
    while covered_length < coverage * total_length: #6 seconds per 3x10^9 string length for mean read length of 1k
        #read_length = random.choice(read_length_distribution)
        read_length = PRECOMPUTE_RANDOM.pop()
        start_position = random.randint(0, total_length - read_length)
        #read = fasta[start_position:start_position + read_length]
        #reads.append(read)
        covered_length += read_length
        read_coordinates["Read_%s" % len(read_coordinates.keys())] =[start_position, start_position + read_length]
    return read_coordinates #reads optional
'''
def generate_reads_based_on_coverage_parallel(fasta, read_length_distribution, coverage, num_processes):
    """
    Parallel version of generate_reads_based_on_coverage using multiprocessing.
    """
    print("Pulling reads...")

    #total_length = len(fasta)
    #covered_length = 0
    #read_coordinates = {}

    # Create a multiprocessing pool
    pool = multiprocessing.Pool(processes=num_processes)

    # Generate reads in parallel
    results = pool.apply_async(generate_reads_based_on_coverage, args=(fasta, read_length_distribution, coverage))
    
    # Close the pool
    pool.close()
    pool.join()

    return results.get()
'''
#@timeit
def generate_read_length_distribution(num_reads, mean_read_length, distribution='lognormal', **kwargs):
    """
    Generate a list of read lengths with a customizable mean read length and distribution.
    """
    if distribution == 'normal':
        # Generate read lengths from a normal distribution
        read_lengths = np.random.normal(mean_read_length, mean_read_length / 10, num_reads)
    elif distribution == 'lognormal':
        # Generate read lengths from a log-normal distribution
        read_lengths = np.random.lognormal(mean=np.log(mean_read_length), sigma=1.0, size=num_reads)
    else:
        raise ValueError("Unsupported distribution. Supported options: 'normal', 'lognormal'.")
    
    # Ensure that all read lengths are positive integers
    read_lengths = np.round(np.abs(read_lengths)).astype(int)
    # Filter out zero-length reads
    read_lengths = read_lengths[read_lengths > 0]
    
    return read_lengths #not returning list

#Part 3: Simulate how many insertions can be found using our sequencing approach with different parameter settings
#@timeit
def count_insertions(insertion_dir, read_dir):
    """
    Count the number of full-length and partial insertions for each insertion.
    """
    data = []

    for insertion, insertion_positions in insertion_dir.items():
        full_length_count = 0
        partial_count = 0
        insertion_start, insertion_end = insertion_positions

        for read, read_positions in read_dir.items():
            read_start, read_end = read_positions

            if read_start <= insertion_start and insertion_end <= read_end:  # Full-length insertion
                full_length_count += 1
            elif (insertion_start < read_start and insertion_end > read_start) or \
                 (insertion_start < read_end and insertion_end > read_end):  # Partial insertion
                partial_count += 1
        
        data.append({'Insertion': insertion, 'Full_Length_Matches': full_length_count, 'Partial_Matches': partial_count})

    df = pd.DataFrame(data)
    full=df['Full_Length_Matches'].sum()
    partial=df['Partial_Matches'].sum()
    return df, full, partial

'''
### Calling everything in the right order
#1 Creates the vector in the right format
insertion_fasta = collapse_fasta(vector_sequence_path)
#2 Creates a single-string reference genome, while storing the chromosome borders
fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
#3 Randomly inserts insertio sequence into the single-string reference x times
mod_fasta, insertion_dir = add_insertions_to_genome_sequence(fasta, insertion_fasta, insertion_numbers)
#4 Creates an overview table of the insertion distribution
insertions_df = insertions_per_chromosome(chromosome_dir, insertion_dir)
print(insertions_df)
#5 creates histogram from actual sequenced distribution and the randomly pulled reads
sequence_read_length_distribution = get_read_length_distribution_from_real_data(sequenced_data_path)
plt.hist(sequence_read_length_distribution, bins=400)
reads = generate_reads(mod_fasta, sequence_read_length_distribution, num_reads=len(sequence_read_length_distribution))
plt.hist([len(i) for i in reads], bins=400)

cov_reads, cov_coordinates = generate_reads_based_on_coverage(mod_fasta, sequence_read_length_distribution, coverage=coverage)
plt.hist([len(i) for i in cov_reads], bins=400)
plt.savefig('sequence_based_read_length_distribution.pdf')
plt.close()

#6 creates histogram for artificially created distribution and the randomly pulled reads
custom_read_length_distribution = generate_read_length_distribution(500000, mean_read_length = mean_read_length, distribution = 'lognormal') #currently only normal and lognormal work
plt.hist(custom_read_length_distribution, bins=400)

custom_cov_reads, custom_cov_coordinates = generate_reads_based_on_coverage(mod_fasta, custom_read_length_distribution, coverage=coverage)
plt.hist([len(i) for i in custom_cov_reads], bins=400)
plt.savefig('calculated_read_length_distribution.pdf')
#7 creates final table with detected insertions
detected, full, partial = count_insertions(insertion_dir, custom_cov_coordinates)
print("Total Reads: " + str(len(custom_cov_reads)))
print(full)
print(partial)
'''
t0 = time.time()
coverages = [1, 5, 10, 15, 20] #* 100 #ten replicates?
mean_read_lengths = [1000, 2000, 3000, 4000, 5000, 6000,7000,8000,9000,10000,15000, 20000, 25000, 30000]
#mean_read_lengths = [12000]
#coverages=[1] * 10 #,5,10] #* 10 #coverage with some influence on the runtime
#mean_read_lengths=[5000, 6000]
combinations = itertools.product(mean_read_lengths, coverages)
#coverages=[1,2]
#mean_read_lengths=[5000, 6000]

#run once
#1 Creates the vector in the right format
#insertion_fasta = collapse_fasta(vector_sequence_path)
insertion_fasta = 'X' * 5000 #artificial insertion
#2 Creates a single-string reference genome, while storing the chromosome borders
fasta, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path)
#3 import bed with possible insertion positions 
bed_df = readbed(bedpath, chromosome_dir.keys())
#3 Randomly inserts insertio sequence into the single-string reference x times
mod_fasta, insertion_dir = add_insertions_to_genome_sequence_with_bed(fasta, insertion_fasta, insertion_numbers, chromosome_dir, bed_df)

del fasta
del insertion_fasta
del bed_df


# Define a function to process each combination
def process_combination(mean_read_length, coverage, insertion_numbers, mod_fasta, insertion_dir):
    
    #Creates a read length distribution based mean read length and draws artificial reads from it based on coverage settings. Then it compares the coordinates of the previously created
    #insertions in the reference fasta and checks whether they are partially or fully contained within the artificial reads.
    
    #custom_read_length_distribution = get_read_length_distribution_from_real_data(sequenced_data_path) #for experimental data
    print('Calculating for:' + str(mean_read_length))
    custom_read_length_distribution = generate_read_length_distribution(1000000, mean_read_length=mean_read_length, distribution='lognormal')
    PRECOMPUTE_RANDOM = [random.choice(custom_read_length_distribution) for _ in range(100000000)]
    #custom_cov_coordinates = generate_reads_based_on_coverage_parallel(mod_fasta, custom_read_length_distribution, coverage=coverage, num_processes=10)
    custom_cov_coordinates = generate_reads_based_on_coverage(mod_fasta, custom_read_length_distribution, coverage, PRECOMPUTE_RANDOM)
    detected, full, partial = count_insertions(insertion_dir, custom_cov_coordinates)
    return [mean_read_length, coverage, insertion_numbers, full, partial]
    #save memory from overload
    del custom_read_length_distribution
    del custom_cov_coordinates
    gc.collect()

'''
results = [] #1270 without parallelization for 20
for mean_read_length, coverage in combinations:
    out = process_combination(mean_read_length, coverage, insertion_numbers, mod_fasta, insertion_dir)
    results.append(out)

'''
def parallel_duplicates(n_interations): #x100 down to ~5.5h
    results=[]
    print(n_interations)
    for mean_read_length, coverage in combinations:
        out = process_combination(mean_read_length, coverage, insertion_numbers, mod_fasta, insertion_dir)
        results.append(out)
    return results

parallel_results = Parallel(n_jobs=20)(delayed(parallel_duplicates)(i) for i in range(100))
finaloutput = list(itertools.chain(*parallel_results))

'''
def process_combination_wrapper(args):
    return process_combination(*args)
# Create a multiprocessing pool
pool = multiprocessing.Pool(processes=20) #x100: 148530.24656414986 s with 25 cores

args_list = [(mean_read_length, coverage, insertion_numbers, mod_fasta, insertion_dir) for mean_read_length, coverage in combinations]

print(args_list)
# Parallel execution of process_combination function for each combination
#results = []
#for _ in range(10):
#    results.extend(pool.map(process_combination_wrapper, args_list))
# Retrieve results
#finaloutput = []
#for result in results:
#    finaloutput.append(result.get())
#    print(result.get())
# Close the pool and wait for all processes to finish
pool.close()
pool.join()
'''
#finaloutput = results
# Convert finaloutput to DataFrame
finaloutput_df = pd.DataFrame(finaloutput, columns=['mean_read_length', 'coverage', 'insertion_number', 'full_matches', 'partial_matches'])

print(finaloutput_df)
t1 = time.time()
total = t1-t0
print(total)
# Write DataFrame to CSV
finaloutput_df.to_csv(output_path, sep='\t', header=True, index=False)
