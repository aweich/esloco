import math
import logging
import pandas as pd

def bedcheck(bed):
    """
    Check if the first column of a df contains the string "chr" and if the second and thirs column contain an integer.  
    """
    if bed.iloc[:, 0].str.startswith("chr").all() and \
        bed.iloc[:, 1].apply(lambda x: isinstance(x, int)).all() and \
        bed.iloc[:, 2].apply(lambda x: isinstance(x, int)).all():
        return True
    logging.info("The file does not seem to be in BED format. Make sure the data looks like: chrN integer integer")
    return False

def extract_bed_columns(df):
    """
    Extracts the first three columns of a DataFrame and assigns default column names if missing.
    Assumes the first three columns correspond to 'chrom', 'start', and 'end'.
    """
    # Extract the first three columns
    df = df.iloc[:, :3].copy()
    # Assign default column names if they don't exist
    df.columns = ["chrom", "start", "end"]
    return df

def readbed(bedpath, list_of_chromosomes_in_reference, barcoding=False):
    """
    Reads bed files in a flexible manner. Adds None fr columns that are non-existent. 
    """
    bed = None
    try:
        # Read only the available columns initially
        bed = pd.read_csv(bedpath, sep="\t", header=None, dtype=str)
        # Check if first row is a header (contains "chr" and column 1 is "start")
        if "chr" in bed.iloc[0, 0].lower() and bed.iloc[0, 1].lower() == "start":
            bed = pd.read_csv(bedpath, sep="\t", header=0, dtype=str)
        else:
            logging.info(f"No header in {bedpath}. Extracting columns 1-3 and assigning default column names.")
            bed = extract_bed_columns(bed)
        # Convert start and end columns to integer
        bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
        bed["end"] = pd.to_numeric(bed["end"], errors="coerce")
        if not bedcheck(bed):
            return None
        bed.columns = bed.columns.str.lower()
        all_cols = ["chrom", "start", "end", "id", "barcode", "weight"]
        missing_cols = list(set(all_cols) - set(bed.columns))
        logging.info(f"Your BED is missing the following columns:{missing_cols}. Trying to fill in default values where possible.")
        for missing_col in missing_cols:
            if missing_col == "weight":
                bed[missing_col] = 1
            elif missing_col == "id":
                bed[missing_col] = bed["chrom"].astype(str) + "_" + bed["start"].astype(str) + "_" + bed["end"].astype(str) + "_" + bed["weight"] .astype(str)
            else:
                bed[missing_col] = [[]] * len(bed)
        # Fill empty or NaN values in the 'weight' column with the default value of 1
        bed['weight'] = bed['weight'].apply(lambda x: 1 if x in ["", None] or (isinstance(x, float) and math.isnan(x)) else x)
        # Fill NaN values in the 'Barcode' column with an empty list
        bed['barcode'] = bed['barcode'].apply(lambda x: [] if x in ["", None] or (isinstance(x, float) and math.isnan(x)) else x)
        # Apply barcoding transformation if required
        if barcoding:
            logging.info('Barcoding selected: Transforming the chromosome names in the bed...')
            barcode = '_'.join(list(list_of_chromosomes_in_reference)[0].split("_")[:-1])
            bed['chrom'] = barcode + '_' + bed['chrom'].astype(str)
        # Filter out chromosomes not in the reference list
        logging.info(f'Only keeping the following chromosomes: {", ".join(sorted(list_of_chromosomes_in_reference))}')
        bed = bed[bed["chrom"].isin(list_of_chromosomes_in_reference)]
    except Exception as e:
        logging.info(f"Not defined or error reading the BED file: {e}")
        return None
    return bed

def chromosome_to_global_coordinates(beddf, input_chromosome_dict):
    '''
    Uses a df (bed-like) and a dict of chromosome broders in global format and returns the global cooridnates of the bed entries.
    Checks if start and stop are 0s and if so, sets the coordinates to the full chromosome length for this row. 
    '''
    try:
        updated_coordinates = {}
        #barcode_exists = 'Barcode' in beddf.columns
        #weight_exists = 'weight' in beddf.columns
        for _, row in beddf.iterrows():
            ID=row["id"]
            chrom = row['chrom']
            if row["start"] == row["end"] == 0: #full chromosome mode
                logging.info(f"Start and End coordinates in are both 0. Full {chrom} used for {ID}.")
                start = input_chromosome_dict[chrom][0]
                end = input_chromosome_dict[chrom][1]
            else:
                start = row['start'] + input_chromosome_dict[chrom][0]
                end = row['end'] + input_chromosome_dict[chrom][0]
            # Prepare the entry for updated_coordinates based on conditions
            entry = {'start': start, 'end': end, 'weight': row['weight'], 'barcode': row['barcode']}
            updated_coordinates[ID] = entry
        return updated_coordinates
    except:
        logging.warning("BED coordinates could not be updated.")
        return None

def global_to_chromosome_coordinates(chromosome_dict, items_dict):
    '''
    Converts items dictionary to DataFrame in BED format using chromosome dictionary.
    '''
    bed_data = []
    for item, coordinates in items_dict.items():
        global_start, global_end = coordinates
        for chrom, dimensions in chromosome_dict.items():
            chr_start, chr_end = dimensions
            if chr_start < global_start < chr_end:
                bed_data.append([chrom, global_start - chr_start, global_end - chr_start, item])
    beddf = pd.DataFrame(bed_data, columns=['chrom', 'start', 'end', 'item'])
    return beddf
