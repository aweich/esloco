from create_insertion_genome import add_insertions_to_genome_sequence_with_bed
from utils import profile, check_barcoding, roi_barcoding, barcode_genome
from fasta_operations import pseudo_fasta_coordinates
from bed_operations import readbed, chromosome_to_global_coordinates

@profile
def create_barcoded_insertion_genome(reference_genome_path, bedpath, blocked_regions_bedpath, restriction, insertion_length, insertion_numbers, insertion_number_distribution, n_barcodes):
    '''
    Pre-processment step of the insertion mode.

    1.) The reference genome cooridnates are transformed into a string-like format (One single FASTA string) and the chromsome borders are stored
    2.) Based on user input, masked regions are defined and weighted
    3.) For each barcode, the chromosome borders get a prefix and inserttion positions are randomly chosen. If an insertion-bed is defined, the insertions are only placed within these regions accoridng to their length-based weights.
    4.) The length of the reference genome (for coverage estimation), the insertion cooridnates, the masked regions, and the chromosome border dict are returned
    ''' 
    print("Create Barcoded Genome...")
    collected_insertion_dict={}
    #1
    ref_genome_size, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path, restriction)
    
    #2 #create blocked regions file
    if blocked_regions_bedpath is not None:
        blocked_bed = readbed(blocked_regions_bedpath, chromosome_dir.keys())
        masked_regions = chromosome_to_global_coordinates(blocked_bed, chromosome_dir)
    else:
        masked_regions = None
    
    print(masked_regions)

    #3
    for i in range(n_barcodes):
        barcoded_chromosome_dir = barcode_genome(chromosome_dir, i)
        #optional bed-guided insertion
        bed_df = readbed(bedpath, barcoded_chromosome_dir.keys(), barcoding=check_barcoding(n_barcodes))
        if not bed_df:
            print("Insertions will be placed randomly...")
        genome_size, insertion_dict = add_insertions_to_genome_sequence_with_bed(ref_genome_size, insertion_length, insertion_numbers, barcoded_chromosome_dir, insertion_number_distribution, bed_df)
        collected_insertion_dict.update(insertion_dict)
        del barcoded_chromosome_dir, bed_df, insertion_dict
    #4
    return genome_size, collected_insertion_dict, masked_regions, chromosome_dir

def create_barcoded_roi_genome(reference_genome_path, restriction, roi_bedpath, n_barcodes, blocked_regions_bedpath):
    '''
    Pre-processment step of the roi mode.

    1.) The reference genome cooridnates are transformed into a string-like format (One single FASTA string) and the chromsome borders are stored
    2.) Based on user input, masked regions are defined and weighted
    3.) The pre-defined ROIs are barcoded.
    4.) The length of the reference genome (for coverage estimation), the insertion cooridnates, the masked regions, and the chromosome border dict are returned
    '''
    #create global cooridnates from bed based on provided genome ref to adjust ROIs to string-like genome
    genome_size, chromosome_dir = pseudo_fasta_coordinates(reference_genome_path, restriction)
    bed = readbed(roi_bedpath, chromosome_dir.keys())
    roi_dict = chromosome_to_global_coordinates(bed, chromosome_dir)
    roi_dict = roi_barcoding(roi_dict, n_barcodes)
    
    #create blocked regions file
    if blocked_regions_bedpath is not None:
        blocked_bed = readbed(blocked_regions_bedpath, chromosome_dir.keys()) #barcoding should be default false
        masked_regions = chromosome_to_global_coordinates(blocked_bed, chromosome_dir)
    
    return roi_dict, bed, masked_regions, genome_size, chromosome_dir