import os
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d

from utils import track_usage

def bin_coverage(coverage, bin_size):
    """
    Aggregate coverage data into bins of a specified size.
    """
    num_bins = int(np.ceil(len(coverage) / bin_size))
    binned_coverage = np.zeros(num_bins)
    
    for i in range(num_bins):
        start = i * bin_size
        end = start + bin_size
        binned_coverage[i] = np.mean(coverage[start:end])
    
    return binned_coverage

def plot_reads_coverage(ref_length,bin_size, reads_dict, mean_read_length, current_coverage, insertion_dict, outputpath):
    """
    Plots a coverage-like plot using the reference genome and reads information.
    """

    if not outputpath or outputpath.lower() == "none":
        logging.info("Skipping plot due to missing output path.")
        return  # Stop execution
    
    output_file = f"{outputpath}/{mean_read_length}_{current_coverage}_coverage.png"
    
    if os.path.exists(output_file):
        logging.info(f"Skipping plot generation: Output file '{output_file}' already exists.")
        return  # Stop execution
    else:
        os.makedirs(outputpath, exist_ok=True)
        logging.info(f"Generating plot: {output_file}")

    #start plotting

    smooth_sigma = 3
    # Initialize coverage array
    coverage = np.zeros(ref_length)
    
    # Extract unique suffixes and assign colors
    unique_suffixes = set(read_id.split('_')[-1] for read_id in reads_dict.keys())
    colors = list(mcolors.TABLEAU_COLORS.keys())
    suffix_color_map = {suffix: colors[i % len(colors)] for i, suffix in enumerate(unique_suffixes)}
    
    # Populate coverage array based on reads
    for read_id, (start, stop) in reads_dict.items():
        suffix = read_id.split('_')[-1]
        coverage[start:stop] += 1
    
    # Bin the coverage
    binned_coverage = bin_coverage(coverage, bin_size)
    bin_positions = np.arange(0, ref_length, bin_size)
    
    # Smooth the binned coverage using a Gaussian filter
    smoothed_binned_coverage = gaussian_filter1d(binned_coverage, sigma=smooth_sigma)

    # Plot the binned coverage
    plt.figure(figsize=(20, 6))
    plt.plot(bin_positions[:len(smoothed_binned_coverage)], smoothed_binned_coverage, drawstyle='steps-pre', color='gray', alpha=0.5)
    #plt.plot(bin_positions[:len(binned_coverage)], binned_coverage, color='r', marker='o')
    
    # Plot individual reads with different colors based on suffix
    for suffix in unique_suffixes:
        coverage_suffix = np.zeros(ref_length)
        for read_id, (start, stop) in reads_dict.items():
            if read_id.endswith(suffix):
                coverage_suffix[start:stop] += 1
        binned_coverage_suffix = bin_coverage(coverage_suffix, bin_size)
        smoothed_binned_coverage_suffix = gaussian_filter1d(binned_coverage_suffix, sigma=smooth_sigma)
        plt.plot(bin_positions[:len(smoothed_binned_coverage_suffix)], smoothed_binned_coverage_suffix, drawstyle='steps-pre', color=suffix_color_map[suffix], label=suffix, alpha=0.5)

    # Add vertical lines at specified positions
    for key, positions in insertion_dict.items():
        suffix = key.split('_')[1]

        if isinstance(positions, dict):
            positions = list(positions.values())  # Convert dict_values to a list
            positions = positions[0:2]
        for pos in positions:
            max_height = max(smoothed_binned_coverage) * 1.1  # Slightly above the highest line
            plt.plot(pos, max_height, marker='v', color=suffix_color_map[suffix], markersize=4, linestyle='None', alpha=0.7)
        

    # Create custom legend
    legend_handles = [Line2D([0], [0], color=suffix_color_map[suffix], lw=2, label=suffix) for suffix in unique_suffixes]
    legend_handles.append(Line2D([0], [0], color='gray', lw=2, label="Combined")) #manual adding of the combined line
    plt.legend(handles=legend_handles, title='Barcode', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.xlabel('Position on "one-string" Reference Genome (1e6 binned)')
    plt.ylabel('Read Coverage')
    plt.title('Read Coverage Plot')
    #plt.show()

    # Save the plot
    track_usage("plot_reads_coverage")
    plt.savefig(output_file)
    plt.close()
    
    logging.info(f"Plot saved as {output_file}")
