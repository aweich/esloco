#!/usr/bin/env python3
import sys
import os
from esloco.config_handler import parse_config
from esloco.plot_functions import read_data, barplot_absolute_matches, barplot_absolute_matches_barcodes, plot_barcode_distribution, plot_lineplot, plot_isolated_lineplot, plot_log_data, generate_html_report

def print_plot_help():
    print("")
    print("Usage: ")
    print("         plot_esloco --config <config_file>")
    print("")
    print("Generate plots from the results of a local coverage estimation simulation and create an HTML report.")
    print("")
    print("Options:")
    print("     --config  <config_file>    Path to the configuration file.")
    print("     --help                     Show this help message.")
    print("")
    sys.exit(0)

def main():
    """ Main function to plot all results at once. """

    # Load configuration
    if len(sys.argv) != 3 or sys.argv[1] == "--help" or sys.argv[1] != "--config":
        print_plot_help()
        sys.exit(1)

    if sys.argv[1] == "--config":
        config_file = sys.argv[2]

    if not os.path.exists(config_file):
        print(f"Configuration file '{config_file}' does not exist.")
        sys.exit(1)

    try:
        param_dictionary = parse_config(config_file)
    except Exception as e:
        print(f"Error parsing config: {e}")
        sys.exit(1)

    # immediate CLI feedback after starting the plotting
    print("Checking config...")

    # output location
    output_path = param_dictionary.get('output_path')
    output_path_plots = param_dictionary.get('output_path_plots')
    experiment_name = param_dictionary.get('experiment_name')
    combinations = param_dictionary.get('combinations')

    if not os.path.exists(output_path_plots):
        #just in case the config has changed between running the simulation and plotting the results
        os.makedirs(output_path_plots)

    basic_data = read_data(f"{output_path}/{experiment_name}_barcode_distribution_table.csv")
    matches_data = read_data(f"{output_path}/{experiment_name}_matches_table.csv")
    log = f"{output_path}/{experiment_name}_log.log"
    coverage_plots = [f"{output_path_plots}/{combination[0]}_{combination[1]}_coverage.html" for combination in combinations]

    #specific plots
    print("Starting plot generation...")
    barplot_absolute= barplot_absolute_matches(experiment_name, matches_data, output_path=output_path_plots)
    barplot_barcodes = barplot_absolute_matches_barcodes(experiment_name, matches_data, output_path=output_path_plots)
    total_reads = plot_barcode_distribution(experiment_name, basic_data, output_path=output_path_plots)
    lineplot_full, lineplot_partial, lineplot_otb = plot_lineplot(experiment_name, matches_data, output_path=output_path_plots)
    lineplot_panel_full, lineplot_panel_partial, lineplot_panel_otb = plot_isolated_lineplot(experiment_name, matches_data, output_path=output_path_plots, filterplots=12) #make command line args at some point
    ressources = plot_log_data(experiment_name, log, output_path=output_path_plots)

    all_plots = [barplot_barcodes,
                 total_reads,
                 lineplot_full, lineplot_partial, lineplot_otb,
                 lineplot_panel_full, lineplot_panel_partial, lineplot_panel_otb,
                 ressources] + coverage_plots

    # HTML report
    # The html report only points to the plots it displays.
    # As the html report is supposed to be placed in the plot output folder, it is necessary to change the relative paths.
    # This whole report generation is not very robust, but it is a quick solution to get an overview of the results.
    all_plots_relative = [os.path.relpath(plot, output_path) for plot in all_plots]
    generate_html_report(all_plots_relative, config=config_file, output_html=f"{output_path}{experiment_name}_report.html")

if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        print("Configuration error:", e)
        sys.exit(1)
