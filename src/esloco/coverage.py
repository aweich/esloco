import os
import logging
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

from esloco.utils import track_usage

def get_barcode_color_mapping(barcodes):
    if isinstance(barcodes, (list, pd.Series, np.ndarray)):
        barcodes = pd.Series(barcodes)
    barcodes = sorted(barcodes)
    colors = px.colors.qualitative.Plotly
    barcode_colors = (colors * (len(barcodes) // len(colors) + 1))[:len(barcodes)]
    return dict(zip(barcodes, barcode_colors))

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


def plot_reads_coverage(ref_length,
                        bin_size,
                        reads_dict,
                        mean_read_length,
                        current_coverage,
                        insertion_dict,
                        outputpath
                        ):
    """
    Plots a coverage-like plot using the reference genome and reads information.
    """
    output_file_svg = f"{outputpath}/{mean_read_length}_{current_coverage}_coverage.svg"
    output_file_html = f"{outputpath}/{mean_read_length}_{current_coverage}_coverage.html"
    if os.path.exists(output_file_svg) and os.path.exists(output_file_html):
        logging.info(f"Skipping plot generation: Output files {output_file_svg} and {output_file_html} already exist.")
        return  # Stop execution
    # Start plotting
    smooth_sigma = 3
    # Initialize coverage array
    coverage = np.zeros(ref_length)
    # Extract unique suffixes and assign colors
    unique_suffixes = set(read_id.split('_')[-1] for read_id in reads_dict.keys())
    suffix_color_map = get_barcode_color_mapping(unique_suffixes)
    # Populate coverage array based on reads
    for read_id, (start, stop) in reads_dict.items():
        suffix = read_id.split('_')[-1]
        coverage[start:stop] += 1
    # Bin the coverage
    binned_coverage = bin_coverage(coverage, bin_size)
    bin_positions = np.arange(0, ref_length, bin_size)
    # Smooth the binned coverage using a Gaussian filter
    smoothed_binned_coverage = gaussian_filter1d(binned_coverage, sigma=smooth_sigma)
    # Create the plotly figure
    fig = go.Figure()
    # Plot the binned coverage
    fig.add_trace(go.Scatter(x=bin_positions[:len(smoothed_binned_coverage)],
                             y=smoothed_binned_coverage,
                             mode='lines',
                             name='Combined',
                             line=dict(color='gray',
                                       width=2)))
    # Plot individual reads with different colors based on suffix
    for suffix in unique_suffixes:
        coverage_suffix = np.zeros(ref_length)
        for read_id, (start, stop) in reads_dict.items():
            if read_id.endswith(suffix):
                coverage_suffix[start:stop] += 1
        binned_coverage_suffix = bin_coverage(coverage_suffix, bin_size)
        smoothed_binned_coverage_suffix = gaussian_filter1d(binned_coverage_suffix, sigma=smooth_sigma)
        fig.add_trace(go.Scatter(x=bin_positions[:len(smoothed_binned_coverage_suffix)],
                                 y=smoothed_binned_coverage_suffix,
                                 mode='lines',
                                 name=suffix,
                                 line=dict(color=suffix_color_map[suffix],
                                           width=2)))
    # Add triangles at specified positions
    for key, positions in insertion_dict.items():
        split_key = key.rsplit("_", 1)  # Split the key into the main part and suffix, handling cases like id_id_id_0, id_0, and id_id_0
        suffix = split_key[1]
        if isinstance(positions, dict):
            positions = list(positions.values())  # Convert dict_values to a list
            positions = positions[0:2]
        # Adjust y value based on suffix
        y_adjustment = 0.01 * int(suffix)  # Adjust dynamically based on suffix
        max_height = max(smoothed_binned_coverage) * (1.1 + y_adjustment)  # Slightly above the highest line
        try:
            fig.add_trace(go.Scatter(x=positions,
                                     y=[max_height] * len(positions),
                                     mode='markers',
                                     name=key,
                                     marker=dict(symbol='triangle-down',
                                                 size=10,
                                                 color=suffix_color_map[suffix])))
        except Exception:
            fig.add_trace(go.Scatter(x=positions,
                                     y=[max_height] * len(positions),
                                     mode='markers',
                                     name=key,
                                     marker=dict(symbol='triangle-down',
                                                 size=10,
                                                 color="black")))
    # Update layout
    fig.update_layout(
        title=f'Read Coverage Plot for {mean_read_length} bp Reads and {current_coverage}x Coverage',
        xaxis_title='Position on "one-string" Reference Genome (1e6 binned)',
        yaxis_title='Read Coverage',
        legend_title='Barcode',
        template='plotly_white',
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,  # Adjust this value to ensure the legend is below the x-axis title
            xanchor="center",
            x=0.5
        )
    )
    # Save the plot as HTML and SVG
    track_usage("plot_reads_coverage")
    pio.write_image(fig, output_file_svg)
    pio.write_html(fig, output_file_html)
    logging.info(f"Plots saved as {output_file_svg} and {output_file_html}")
