import sys
import os

from config_handler import parse_config
from plot_functions import read_data, plot_barcode_distribution, plot_lineplot, plot_isolated_lineplot, plot_log_data, generate_html_report

def main():
    """ Main function to plot all results at once. """

    # Load configuration
    if len(sys.argv) != 2:
        print("Usage: python simplot.py <config_file>")
        sys.exit(1)
    
    config_file = sys.argv[1]
    
    try:
        param_dictionary = parse_config(config_file)
    except Exception as e:
        print(f"Error parsing config: {e}")
        sys.exit(1)
    
    # output location
    output_path = param_dictionary.get('output_path')
    output_path_plots = param_dictionary.get('output_path_plots')
    experiment_name = param_dictionary.get('experiment_name')
    
    if not os.path.exists(output_path_plots):
        #just in case the config has changed between running the simulation and plotting the results
        os.makedirs(output_path_plots)

    basic_data = read_data(f"{output_path}/{experiment_name}_barcode_distribution_table.csv")
    matches_data = read_data(f"{output_path}/{experiment_name}_matches_table.csv")
    log = f"{output_path}/{experiment_name}_log.log"
    
    #specific plots
    print("Starting plot generation...")
    plot1, plot2 = plot_barcode_distribution(basic_data, output_path=output_path_plots)
    plot3, plot4 = plot_lineplot(matches_data, output_path=output_path_plots)
    plot5, plot6 = plot_isolated_lineplot(matches_data, output_path=output_path_plots, filter=20) #make command line args possible? 
    plot7 = plot_log_data(log, output_path=output_path_plots)
    all_plots = [plot1, plot2, plot3, plot4, plot5, plot6, plot7]
    
    # HTML report
    generate_html_report(all_plots, config=config_file, output_html=f"{experiment_name}_report.html")

if __name__ == "__main__":
	try:
		main()
	except ValueError as e:
		print("Configuration error:", e)
		sys.exit(1)