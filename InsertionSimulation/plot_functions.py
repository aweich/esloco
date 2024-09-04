#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import sys
import scipy
from matplotlib.patches import Patch
import re
#coverage plot
import pysam
import ast


out_dir="./out/Noise_Simulation/plots/"
inputdata="./out/Noise_Simulation/Barcodes_100_Introns_with_50kb_blocked_AND_200_Overlap_AND_Poisson5_matches_table.csv"
prefix="Barcodes_100_Introns_with_50kb_blocked_AND_200_Overlap_AND_Poisson5" #'Weight_1_I_DominanceSimulation' #"Combined_" #'Weight_4_I_DominanceSimulation' #sample name for output plot
mode="I"

def plot_matches(data, out_dir):
	"""
	Plot partial and full matches over the mean read lengths with group ID being the coverage and save the plot as a JPG file.
	"""
	# Create a seaborn catplot for full matches
	plt.figure()
	g = sns.catplot(data=data, x='mean_read_length', y='full_matches', hue='coverage', kind='point', errorbar='sd', legend_out=False)
	g.tick_params(axis='x', labelrotation=90)
	g.despine(left=True)
	g.set_axis_labels('Mean Read Length', 'Full Matches')
	g.legend.set_title('Coverage')
	#g.set(yticks=(0,5,25,50,75,100)) #custom ticks
	#g.set_yticklabels(['0','5','25','50','75','100']) 
	# Save the plot as a JPG file
	plotname = prefix + 'full_matches_plot.jpg'
	full_matches_plot_path = os.path.join(out_dir, plotname)
	g.savefig(full_matches_plot_path)

	# Create a seaborn catplot for partial matches
	plt.figure()
	g = sns.catplot(data=data, x='mean_read_length', y='partial_matches', hue='coverage', kind='point', legend_out=False)
	g.tick_params(axis='x', labelrotation=90)
	g.despine(left=True)
	g.set_axis_labels('Mean Read Length', 'Partial Matches')
	g.legend.set_title('Coverage')

	# Save the plot as a JPG file
	plotname = prefix + 'partial_matches_plot.jpg'
	partial_matches_plot_path = os.path.join(out_dir, plotname)
	g.savefig(partial_matches_plot_path)
	plt.close()
 

def lineplot_matches(data, value_column, x_axis="mean_read_length", hue="coverage", out_dir=out_dir):
	"""
	Plot partial and full matches over the mean read lengths with group ID being the coverage and save the plot as JPG files.
	"""
		# Create a line plot for the original data
	plt.figure(figsize=(16, 9))
	plt.title(prefix, y=1.1)
	
	# Plot the original data with a line plot
	g = sns.lineplot(data=data, x=x_axis, y=value_column, hue=hue, legend='full', errorbar='sd', linewidth = 3, palette=sns.color_palette("viridis"), alpha=0.5)
	g.set_xticks(data[x_axis].unique())
	g.tick_params(axis='x', labelrotation=90)
	g.set(xlabel=x_axis, ylabel=value_column)

	# Add a regression line for the original data
	try:
		sns.regplot(data=data, x=x_axis, y=value_column, scatter=False, color='black', line_kws={'linewidth':3, 'linestyle':'--', "alpha":0.8})
	except:
		print("no regression possible")

	# Add the additional points from your table
	new_points_x = [0.51, 6.34, 0.54, 0.79, 0.36, 0.6]
	new_points_y = [4, 44, 5, 3, 3, 4]
	
	# Plot the additional points
	plt.scatter(new_points_x, new_points_y, color='red', label='Approximated from sequencing data.', zorder=5)

	# Add a regression line for the new points
	try:
		sns.regplot(x=new_points_x, y=new_points_y, scatter=False, color='red', ci=None, line_kws={'linewidth':3, 'linestyle':'--', "alpha":0.8})
	except:
		print("no regression possible for new points")

	# Add the legend
	# Add the legend with all labels (original and added regression lines)
	handles, labels = g.get_legend_handles_labels()
	plt.legend(handles=handles + [plt.Line2D([0], [0],linestyle='--', color="black", lw=3, label="Regression based on simulation"), 
	                              plt.Line2D([0], [0],linestyle='--', color="red", lw=3, label="Regression based on sequencing data")],
	           bbox_to_anchor=(1.05, 1), loc='upper left')

	# Save the plot as a JPG file
	plotname = prefix + "_" + value_column + '_lineplot_with_additional_regression.jpg'
	plot_path = os.path.join(out_dir, plotname)
	plt.savefig(plot_path, bbox_inches='tight')
	plt.close()


def plot_points(data, option, color):
	mean_data = data.groupby('coverage')[['full_matches', 'partial_matches']].mean().reset_index()
	if option == "partial_matches":
		#experimental value
		plt.plot(5200, 22, marker='o', color='red', mec='black', mew=1)
		#plt.text(5200 + 500, 22, 'Observed value', verticalalignment='bottom', horizontalalignment='left') #+10 for space
		#for n,i in enumerate(mean_data['coverage']):
		#   plt.plot(data['mean_read_length'].unique(), mean_data['partial_matches'][n], marker='o', color=color[n], mec='black', mew=1)
	elif option == "full_matches":
		plt.plot(5200, 5, marker='o', color='red', mec='black', mew=1)
		plt.text(5200 + 500, 5, 'Observed value', verticalalignment='bottom', horizontalalignment='left') #+10 for space
		#for n,i in enumerate(mean_data['coverage']):
		#   plt.plot(data['mean_read_length'].unique(), mean_data['full_matches'][n],  marker='o', color=color[n], mec='black', mew=1)
	else:
		print("no valid option")

def plot_barcode_barplot_panel(df, match_type, out_dir, prefix):
	"""
	Same plot as below but with panels
	"""
	# Calculate the total full matches for each coverage and mean read length
	df["Barcode"] = df["Insertion"].str.split("_insertion").str[0]

	# Group the DataFrame by coverage
	grouped_df = df.groupby("coverage")

	# Calculate the number of subplots needed based on the number of unique coverage values
	num_subplots = len(grouped_df)
	num_cols = 2  # Number of columns in the subplot grid
	num_rows = (num_subplots + 1) // num_cols  # Number of rows in the subplot grid

	# Set up the figure and axis
	fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
	plt.subplots_adjust(hspace=0.5)  # Adjust vertical space between subplots

	# Flatten the axes array if there's only one row
	if num_rows == 1:
		axes = axes.reshape(1, -1)

	legend_handles = []  # List to store legend handles

	for i, (coverage, group) in enumerate(grouped_df):
		# Calculate the total matches (full or partial) for each coverage, mean read length, and Barcode
		pivot_df = group.pivot_table(index="mean_read_length", columns="Barcode", values=match_type, aggfunc="mean", fill_value=0)

		# Plot stacked bar chart
		ax = axes[i // num_cols, i % num_cols]
		pivot_df.plot(kind='bar', stacked=True, ax=ax, legend=False)  # Remove legend for individual subplot

		# Set labels and title
		ax.set_xlabel("Mean Read Length")
		ax.set_ylabel(match_type.capitalize() + " Matches")
		ax.set_title(f"Barcode Contribution to {match_type.capitalize()} Matches for Coverage {coverage}")


		# Get legend handles for the first subplot
		if i == 0:
			legend_handles, _ = ax.get_legend_handles_labels()

	# Create a single legend outside the subplots
	fig.legend(legend_handles, df["Barcode"].unique(), loc='upper right')


	# Save the plot
	plotname = f"{prefix}_{match_type.capitalize()}_stacked_bar_panel.jpg"
	plot_path = os.path.join(out_dir, plotname)
	plt.savefig(plot_path, bbox_inches='tight')
	plt.close()


def plot_barcode_barplot(df, match_type, out_dir, prefix):
	# Calculate the total full matches for each coverage and mean read length
	#df["Barcode"] = df["Insertion"].str.split("_insertion").str[0]

	# Calculate the total matches (full or partial) for each coverage, mean read length, and Barcode
	pivot_df = df.pivot_table(index=["mean_read_length", "coverage"], columns="Barcode", values=match_type, aggfunc="sum", fill_value=0)

	# Plot stacked bar chart
	ax = pivot_df.plot(kind='bar', stacked=True, figsize=(15, 8))

	# Set labels and title
	ax.set_xlabel("Mean Read Length")
	ax.set_ylabel(match_type.capitalize() + " Matches")
	ax.set_title(f"Barcode Contribution to {match_type.capitalize()} Matches Grouped by Coverage")

	# Create legend
	ax.legend(title="Barcode")

	# Save the plot
	plotname = f"{prefix}_{match_type.capitalize()}_stacked_bar.jpg"
	plot_path = os.path.join(out_dir, plotname)
	plt.title(prefix, y=1.1)
	plt.savefig(plot_path, bbox_inches='tight')
	plt.close()

def lineplot_matches_barcode(data, unique_column, value_column, x_axis="mean_read_length", hue="coverage", palette=None, out_dir=out_dir, prefix=prefix):
	sns.set_style("ticks")
	"""
	Plot partial and full matches over the mean read lengths and coverages. 
	This splits by barcodes, so we essentially get each barcode's individual numbers (as mean per x iterations) as a single line. 
	This shows that there are no preferred barcodes.
	"""
	# Create unique plots for each unique value in the unique_column
	unique_values = data[unique_column].unique()
	unique_values = unique_values[np.argsort(unique_values)]
	num_unique = len(unique_values)
	cols = 7  # Number of columns in the subplot grid
	rows = (num_unique + cols - 1) // cols  # Number of rows in the subplot grid

	fig, axes = plt.subplots(rows, cols, figsize=(25, rows * 5), sharey=True)
	axes = axes.flatten()  # Flatten the 2D array of axes
	handles=[]
	labels=[]
	ordered=data[hue].unique()[np.argsort(data[hue].unique())]
	print(ordered)


	for n, parameter in enumerate(unique_values):
		ax = axes[n]
		# Capture handles and labels for the legend from the first plot
		parameter_data = data[data[unique_column] == parameter]
		if n == 0:
			sns.lineplot(data=parameter_data, x=x_axis, y=value_column, hue=hue, hue_order=ordered,  palette=palette, legend=True, ax=ax, errorbar=None)
			handles, labels = ax.get_legend_handles_labels()
			ax.get_legend().remove()
		else:
			sns.lineplot(data=parameter_data, x=x_axis, y=value_column, hue=hue,hue_order=ordered,  palette=palette, legend=False, ax=ax, errorbar=None)
		ax.set_xticks(data[x_axis].unique())
		ax.tick_params(axis='x', labelrotation=90)
		ax.set(xlabel=x_axis, ylabel=value_column)
		ax.set_title(f"{unique_column}: {parameter}")

	# Remove any unused subplots
	for i in range(n + 1, len(axes)):
		fig.delaxes(axes[i])

	# Create a single legend for the hue
	fig.legend(handles, labels, bbox_to_anchor=(0.9, 1), loc='upper left', title=hue)

	# Adjust layout
	plt.tight_layout(rect=[0, 0, 0.85, 1])

	# Save the figure
	plotname = f"{prefix}_{unique_column}_{value_column}_barcode_lineplot.jpg"
	plot_path = os.path.join(out_dir, plotname)
	plt.savefig(plot_path, bbox_inches='tight')
	plt.close()

def plot_insertions(data, insertion_name, value_column_1, value_column_2, out_dir, prefix=""):
	"""
	Plot two value columns for each insertion side by side.

	"""
	# Filter data for the specified insertion
	insertion_data = data[data['ROI'] == insertion_name]

	# Plot first value column
	plt.figure(figsize=(15, 6))
	plt.subplot(1, 2, 1)
	sns.lineplot(x='mean_read_length', y=value_column_1, hue='coverage', data=insertion_data)
	plt.title(f'{value_column_1} for {insertion_name}')
	plt.xlabel('Mean Read Length')
	plt.ylabel(value_column_1)
	plt.legend(title='Coverage')

	# Plot second value column
	plt.subplot(1, 2, 2)
	sns.lineplot(x='mean_read_length', y=value_column_2, hue='coverage', data=insertion_data)
	plt.title(f'{value_column_2} for {insertion_name}')
	plt.xlabel('Mean Read Length')
	plt.ylabel(value_column_2)
	plt.legend(title='Coverage')

	# Save plot
	plotname = prefix + f"ROI_{insertion_name}_lineplot.jpg"
	plotpath = os.path.join(out_dir, plotname)
	plt.savefig(plotpath, bbox_inches='tight')
	plt.close()

def extract_numbers_from_filename(filename):
	numbers = re.findall(r'\d+', filename)
	if len(numbers) > 1:
		print(numbers)
		print(int(numbers[0]))
		print("more than one int in filename. Using the first one...")
		return int(numbers[0])
		#sys.exit()
	else:
		return int(numbers[0])

def combine_files_with_id(input_files):
	combined_df = pd.DataFrame()
	
	for idx, file_path in enumerate(input_files):
		print(file_path)
		df = pd.read_csv(file_path, sep='\t')
		numbers = extract_numbers_from_filename(file_path)
		print(numbers)
		df['id'] = numbers
		combined_df = pd.concat([combined_df, df], ignore_index=True)
	
	return combined_df

def combine_value_columns(df, value_column_1, value_column_2):
	df["combined_values"] = df[value_column_1] + df[value_column_2]
	return df

def roi_plot_clustermap(df):
	# Step 1: Process the sample names to extract the part before the first '_'
	print(df['Insertion'])
	df['Sample'] = df['Insertion'].str.split("_").str[0]

	# Step 2: Sum full_matches and partial_matches
	df['Total_Matches'] = df['full_matches'] + df['partial_matches']
	df["Normalized_Matches"] = df["Total_Matches"] / df["Length"]
	df_grouped = df.groupby(['Sample', 'coverage']).agg({'Normalized_Matches': 'mean'}).reset_index()

	# Step 4: Pivot the table to have samples as rows and coverage as columns
	df_pivot = df_grouped.pivot(index='Sample', columns='coverage', values='Normalized_Matches')
	print(df_grouped.head())
	print(df_pivot.head())
	# Step 5: Plot the clustermap
	g = sns.clustermap(df_pivot, cmap='viridis', annot=True, dendrogram_ratio=(0.05, 0.05),cbar_kws={'label': 'Mean value across all iterations (sum of matches / length of panel region)'}, row_cluster=False,)
	g.ax_heatmap.set_xlabel('Genome-wide coverage n*6e9')
	g.ax_heatmap.set_ylabel('Panel regions')
	g.ax_row_dendrogram.set_visible(False) #suppress row dendrogram
	g.cax.set_position([1.05, .2, .03, .45])
	# Save plot
	plotname = prefix + "_roi" + '_heatmap.jpg'
	plot_path = os.path.join(out_dir, plotname)
	plt.savefig(plot_path, bbox_inches='tight', dpi=300)
	plt.close()

def overlap_plot(df):
	"""
	Plots a histogram of the overlaps between read and I/ROI 
	"""
	
	#transform into lists of integers and not lists of chracters
	df['overlap'] = df['overlap'].apply(lambda x: ast.literal_eval(x))
	overlap_values = [item for sublist in df['overlap'] for item in sublist]
	sns.histplot(overlap_values, bins=50, kde=False)
	plt.axvline(200, color='black', linestyle="--")
	plt.text(200, plt.ylim()[1]*0.8, f'Min. overlap: 200', color='black', fontsize='xx-small')
	plt.axvline(5000, color='black', linestyle="--")
	plt.text(5000, plt.ylim()[1]*0.8, f'Max. overlap (I/ROI size): 5000', color='black', fontsize='xx-small')
	plt.xlabel('Overlap Values')
	plt.ylabel('Frequency')
	plt.title('Histogram of Overlap Values')

	# Save plot
	plotname = prefix + "_overlap_hist.jpg"
	plot_path = os.path.join(out_dir, plotname)
	plt.savefig(plot_path, bbox_inches='tight', dpi=300)
	plt.close()



# Example usage:
if mode == "I":
	inputdata_df=pd.read_csv(inputdata, sep='\t')
	#MK025_data = pd.read_csv("./out/MK025_experimental_data.csv", sep='\t')

	inputdata_df["Barcode"] = inputdata_df["Insertion"].str.split("_insertion").str[0]
	inputdata_df["Iteration"] = inputdata_df["Insertion"].str.split("_").str[4]


	#1 first combine the detections of an iteration across all barcodes (SUM)
	inputdata_df = combine_value_columns(inputdata_df, "full_matches", "partial_matches")
	print(inputdata_df.head())
	barcodes_summed_up = inputdata_df.groupby(['coverage', 'mean_read_length', 'Iteration']).agg({'full_matches': 'sum', 'partial_matches': 'sum', 'combined_values': 'sum'}).reset_index() #barcodes summed
	print(barcodes_summed_up.head())
	#2 second calculate the mean
	mean_per_iteration = barcodes_summed_up.groupby(['coverage', 'mean_read_length']).agg({'full_matches': 'mean', 'partial_matches': 'mean', 'combined_values': 'mean'}).reset_index()
	print("mean per iteration")
	print(mean_per_iteration)
	lineplot_matches(mean_per_iteration, "combined_values", x_axis='coverage', hue='mean_read_length')

	#3 deal with barcodes
	#combined_mean_df = inputdata_df.groupby(['coverage', 'mean_read_length', 'Barcode']).agg({'combined_values': 'mean'}).reset_index() #mean over iterations
	#print(combined_mean_df.head())
	lineplot_matches(inputdata_df, "partial_matches", x_axis='coverage', hue='Barcode')

	### plot overlap
	print(inputdata_df.head())
	overlap_plot(inputdata_df)
	#old
	sys.exit()
	#combined full and partial matches
	combined_df = combine_value_columns(inputdata_df, "full_matches", "partial_matches")
	print(combined_df)
	combined_mean_df = combined_df.groupby(['coverage', 'mean_read_length', 'Barcode']).agg({'combined_values': 'mean'}).reset_index() #mean over iterations
	combined_mean_df = combined_mean_df.groupby(['coverage', 'mean_read_length']).agg({'combined_values': 'sum'}).reset_index() #sum over barcode
	print(combined_mean_df)
	combined_mean_df['Coverage_ReadLength'] = combined_mean_df['coverage'].astype(str) + "_" + combined_mean_df['mean_read_length'].astype(str)
	print(combined_mean_df.head())
	#not bad but also not good heatmaps
	#lineplot_matches(combined_mean_df, "combined_values", x_axis='Coverage_ReadLength', hue='Barcode')
	lineplot_matches(combined_mean_df, "combined_values", x_axis='coverage', hue='mean_read_length')
	#sums up the full matches and partial matches across the barcodes but keeps the iterations for the statistics int he plot! = This means that the plot shows all detected insertions across the barcodes, irrespective of the barcode
	summed_df = inputdata_df.groupby(['coverage', 'mean_read_length', 'Iteration']).agg({'full_matches': 'sum', 'partial_matches': 'sum'}).reset_index()
	print(summed_df)
	#lineplot_matches(summed_df, "full_matches", x_axis='mean_read_length', hue='coverage')
	#lineplot_matches(summed_df, "partial_matches", x_axis='mean_read_length', hue='coverage')

	#calculates the mean value for each barcode across the iterations
	mean_df = summed_df.groupby(['coverage', 'mean_read_length']).agg({'full_matches': 'mean', 'partial_matches': 'mean'}).reset_index()
	print(mean_df)
	sys.exit()
	#plot_barcode_barplot(mean_df, "full_matches", out_dir, prefix)
	#plot_barcode_barplot(mean_df, "partial_matches", out_dir, prefix)

	#mean
	#mean of the full matches and partial matches across the the iterations for the statistics int he plot
	mean_df = inputdata_df.groupby(['coverage', 'mean_read_length','Barcode']).agg({'full_matches': 'mean', 'partial_matches': 'mean'}).reset_index()
	#summed_df = summed_df.sort_values(by=['coverage', 'mean_read_length'])
	mean_df['Coverage_ReadLength'] = mean_df['coverage'].astype(str) + "_" + mean_df['mean_read_length'].astype(str)
	print(mean_df)
	#not bad but also not good heatmaps
	lineplot_matches(mean_df, "full_matches", x_axis='Coverage_ReadLength', hue='Barcode')
	lineplot_matches(mean_df.sort_values(by="mean_read_length"),"partial_matches", x_axis='Coverage_ReadLength', hue='Barcode')

elif mode=="ROI":
	inputdata_df=pd.read_csv(inputdata, sep='\t')
	print(inputdata_df.head())
	roi_plot_clustermap(inputdata_df) #qucik chatgpt function
	#overlap plot
	overlap_plot(inputdata_df)


#if combined weighted plot
else:
	#input_files = ["./out/DominanceSimulation/RandomInsertions_Weight_%s_I_DominanceSimulation_matches_table.csv" %i for i in range(1,6)]#
	input_files = ["./out/Noise_Simulation/Barcodes_%i_Introns_with_50kb_blocked_matches_table.csv" %i for i in [1,10,100,1000]]
	print(input_files)
	combined = combine_files_with_id(input_files)
	combined["Barcode"] = combined["Insertion"].str.split("_insertion").str[0]
	combined["Iteration"] = combined["Insertion"].str.split("_").str[4]

	#combine for total matches and then sum up n barcodes per iteration and cov
	combined_df = combine_value_columns(combined, "full_matches", "partial_matches")
	print(combined_df.tail())
	barcodes_summed_up = combined_df.groupby(['coverage', 'mean_read_length', 'Iteration', 'id']).agg({'full_matches': 'sum', 'partial_matches': 'sum', 'combined_values': 'sum'}).reset_index() #barcodes summed
	print(barcodes_summed_up.tail())

	#2 second calculate the mean over iterations
	mean_per_iteration = barcodes_summed_up.groupby(['coverage', 'mean_read_length', 'id']).agg({'full_matches': 'mean', 'partial_matches': 'mean', 'combined_values': 'mean'}).reset_index()
	print("mean per iteration")
	print(mean_per_iteration)
	
	lineplot_matches(mean_per_iteration, "combined_values", x_axis='coverage', hue='id')
	sys.exit()
	lineplot_matches(mean_df, "full_matches", x_axis='Coverage_ReadLength', hue='Barcode') #all matches per weight: "full_matches", x_axis='Coverage_ReadLength', hue='0_weight'
	lineplot_matches(mean_df.sort_values(by="mean_read_length"), "partial_matches", x_axis='Coverage_ReadLength', hue='Barcode')
	#
	lineplot_matches_barcode(mean_df, "0_weight", "full_matches", x_axis='Coverage_ReadLength', hue='Barcode',palette=sns.color_palette("deep"))
	lineplot_matches_barcode(mean_df.sort_values(by="mean_read_length"), "0_weight", "partial_matches", x_axis='Coverage_ReadLength', hue='Barcode', palette=sns.color_palette("deep"))