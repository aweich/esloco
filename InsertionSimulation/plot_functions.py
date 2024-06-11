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

out_dir="./out/DominanceSimulation/plots/"
inputdata="./out/DominanceSimulation/Homogeneous_I_DominanceSimulation_matches_table.csv"
prefix="Homogeneous_Fixed_" #'Weight_1_I_DominanceSimulation' #"Combined_" #'Weight_4_I_DominanceSimulation' #sample name for output plot
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
 

def lineplot_matches(data, value_column, x_axis="mean_read_length", hue="coverage", out_dir=out_dir):
	"""
	Plot partial and full matches over the mean read lengths with group ID being the coverage and save the plot as JPG files.
	"""
	# Create a line plot for full matches
	plt.figure()
	plt.title(prefix, y=1.1)
	g = sns.lineplot(data=data, x=x_axis, y=value_column, hue=hue, legend='full', errorbar='sd')
	g.set_xticks(data[x_axis].unique())
	g.tick_params(axis='x', labelrotation=90)
	g.set(xlabel=x_axis, ylabel=value_column)
	g.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=hue)

	# Save the plot as a JPG file
	#add experimental line
	#plt.axvline(x=5200, color='black', linestyle=':')
	#plt.text(5200, plt.ylim()[1], 'MK025', verticalalignment='bottom', horizontalalignment='center')
	# Add point at y=5 on the vertical line
	#plot_points(MK025_data, option=value_column, color=sns.cubehelix_palette())

	plotname = prefix + "_" + value_column + '_lineplot.jpg'
	plot_path = os.path.join(out_dir, plotname)
	plt.savefig(plot_path, bbox_inches='tight')

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
		print("more than one int in filename")
		sys.exit()
	else:
		return int(numbers[0])

def combine_files_with_id(input_files):
	combined_df = pd.DataFrame()
	
	for idx, file_path in enumerate(input_files):
		print(file_path)
		df = pd.read_csv(file_path, sep='\t')
		numbers = extract_numbers_from_filename(file_path)
		print(numbers)
		df['0_weight'] = numbers
		combined_df = pd.concat([combined_df, df], ignore_index=True)
	
	return combined_df

def create_heatmap(df,id, value_column):

	# Create a unique identifier for rows
	df['identifier'] = df.apply(lambda x: f"{x['Barcode']}_{x[id]}", axis=1)
	
	# Create pivot tables
	matches_pivot = df.pivot(index='identifier', columns=['coverage', 'mean_read_length'], values=value_column)
	matches_pivot = matches_pivot.fillna(0)
	# Create a DataFrame for row colors
	row_colors_df = df[['identifier', 'Barcode', id]].drop_duplicates().set_index('identifier')

	# Map colors to unique values of 'Barcode' and '0_weight'
	unique_barcodes = row_colors_df['Barcode'].unique()
	barcode_colors = sns.color_palette('hsv', len(unique_barcodes))
	barcode_color_dict = dict(zip(unique_barcodes, barcode_colors))

	unique_id = row_colors_df[id].unique()
	id_colors = sns.color_palette('coolwarm', len(unique_id))
	id_color_dict = dict(zip(unique_id, id_colors))

	# Add color annotations
	row_colors_df['Barcode_color'] = row_colors_df['Barcode'].map(barcode_color_dict)
	row_colors_df['id_color'] = row_colors_df[id].map(id_color_dict)
	row_colors = row_colors_df[['Barcode_color', 'id_color']]

	# Create custom legend handles
	barcode_handles = [Patch(facecolor=barcode_color_dict[bc], edgecolor='black', label=bc) for bc in unique_barcodes]
	id_handles = [Patch(facecolor=id_color_dict[wt], edgecolor='black', label=f'{wt}') for wt in unique_id]

	# Plot the heatmaps with annotations
	sns.set(font_scale=0.7)
	g = sns.clustermap(matches_pivot, row_colors=row_colors, cmap="viridis", figsize=(20, 10), dendrogram_ratio=0.02)
	g.ax_heatmap.set_title(value_column)
	g.ax_heatmap.set_xlabel('Coverage, Mean Read Length')
	g.ax_heatmap.set_ylabel('Barcode, Weight')
	
	g.ax_cbar.set_position((-0.2,0.75,0.05,0.2))
	g.ax_heatmap.legend(handles=barcode_handles + id_handles, title='Legends', bbox_to_anchor=(-0.2, 1), loc='upper left', ncol=1)
	# Save plot
	plotname = prefix + "_" + str(value_column) + '_heatmap.jpg'
	plot_path = os.path.join(out_dir, plotname)
	plt.savefig(plot_path, bbox_inches='tight')



# Example usage:
if mode == "I":
	inputdata_df=pd.read_csv(inputdata, sep='\t')
	MK025_data = pd.read_csv("./out/MK025_experimental_data.csv", sep='\t')

	inputdata_df["Barcode"] = inputdata_df["Insertion"].str.split("_insertion").str[0]
	inputdata_df["Iteration"] = inputdata_df["Insertion"].str.split("_").str[4]

	#sums up the full matches and partial matches across the barcodes but keeps the iterations for the statistics int he plot! = This means that the plot shows all detected insertions across the barcodes, irrespective of the barcode
	summed_df = inputdata_df.groupby(['coverage', 'mean_read_length', 'Iteration']).agg({'full_matches': 'sum', 'partial_matches': 'sum'}).reset_index()
	print(summed_df)
	#plot_matches(inputdata_df, out_dir)
	#plot_combined_matches(finputdata_df, out_dir)
	lineplot_matches(summed_df, "full_matches", x_axis='mean_read_length', hue='coverage')
	lineplot_matches(summed_df, "partial_matches", x_axis='mean_read_length', hue='coverage')
	#lineplot_matches_barcode(summed_df, "full_matches", out_dir)
	#lineplot_matches_barcode(summed_df, "partial_matches", out_dir)
	#lineplot_matches(inputdata_df, "full_matches", out_dir)
	#lineplot_matches(inputdata_df, "partial_matches", out_dir)

	#calculates the mean value for each barcode across the iterations
	print(inputdata_df)
	mean_df = inputdata_df.groupby(['coverage', 'mean_read_length', 'Barcode']).agg({'full_matches': 'sum', 'partial_matches': 'sum'}).reset_index()
	print(mean_df)
	#plot_barcode_barplot(mean_df, "full_matches", out_dir, prefix)
	#plot_barcode_barplot(mean_df, "partial_matches", out_dir, prefix)
	create_heatmap(mean_df, "coverage", "partial_matches")

	#mean
	#mean of the full matches and partial matches across the the iterations for the statistics int he plot
	mean_df = inputdata_df.groupby(['coverage', 'mean_read_length','Barcode']).agg({'full_matches': 'mean', 'partial_matches': 'mean'}).reset_index()
	#summed_df = summed_df.sort_values(by=['coverage', 'mean_read_length'])
	mean_df['Coverage_ReadLength'] = mean_df['coverage'].astype(str) + "_" + mean_df['mean_read_length'].astype(str)
	print(mean_df)
	#not bad but also not good heatmaps
	lineplot_matches(mean_df, "full_matches", x_axis='Coverage_ReadLength', hue='Barcode')
	lineplot_matches(mean_df.sort_values(by="mean_read_length"), "partial_matches", x_axis='Coverage_ReadLength', hue='Barcode')

elif mode=="ROI":
	inputdata_df=pd.read_csv(inputdata, sep='\t')
	print(inputdata_df.head())
	#inputdata_df["Barcode"] = inputdata_df["Insertion"].str.rsplit(pat="_", n=2).str[1]
	inputdata_df["ROI"] = inputdata_df["Insertion"].str.rsplit(pat="_", n=2).str[0]
	inputdata_df["Barcode"] = inputdata_df["Insertion"].str.rsplit(pat="_", n=2).str[1]
	print(inputdata_df)
	#aggregate over iterations (=replications)
	grouped = inputdata_df.groupby(['ROI','Barcode', 'coverage', 'mean_read_length']).agg({
		'full_matches': 'sum',
		'partial_matches': 'sum'
	}).reset_index()
	print(grouped)
	#plot_barcode_barplot(grouped, "full_matches", out_dir, prefix)
	#plot_barcode_barplot(grouped, "partial_matches", out_dir, prefix)
	#does not work and needs a better plan!
	plot_insertions(grouped,"BAP1", "partial_matches", "full_matches", out_dir)

	#sums up the full matches and partial matches across the barcodes but keeps the iterations for the statistics int he plot! = This means that the plot shows all detected insertions across the barcodes, irrespective of the barcode
	#inputdata_df["Iteration"] = inputdata_df["Insertion"].str.rsplit(pat="_", n=1).str[1]
	#summed_df = inputdata_df.groupby(['coverage', 'mean_read_length', 'Iteration']).agg({'full_matches': 'sum', 'partial_matches': 'sum'}).reset_index()
	#print(summed_df)
	#plot_matches(inputdata_df, out_dir)
	#plot_combined_matches(finputdata_df, out_dir)
	#lineplot_matches(summed_df, "full_matches", out_dir)
	#lineplot_matches(summed_df, "partial_matches", out_dir)

else:
	#input_files = ["./out/DominanceSimulation/RandomInsertions_Weight_%s_I_DominanceSimulation_matches_table.csv" %i for i in range(1,6)]#
	input_files = ["./out/DominanceSimulation/Weight_%s_I_DominanceSimulation_matches_table.csv" %i for i in [1,2,3,4,5,10,20]]
	print(input_files)
	combined = combine_files_with_id(input_files)
	combined["Barcode"] = combined["Insertion"].str.split("_insertion").str[0]
	combined["Iteration"] = combined["Insertion"].str.split("_").str[4]
	print(combined)

	#mean of the full matches and partial matches across the the iterations for the statistics int he plot
	mean_df = combined.groupby(['coverage', 'mean_read_length','Barcode', '0_weight']).agg({'full_matches': 'mean', 'partial_matches': 'mean'}).reset_index()
	#summed_df = summed_df.sort_values(by=['coverage', 'mean_read_length'])
	mean_df['Coverage_ReadLength'] = mean_df['coverage'].astype(str) + "_" + mean_df['mean_read_length'].astype(str)
	print(mean_df)
	#not bad but also not good heatmaps
	create_heatmap(mean_df, "0_weight", "full_matches")
	create_heatmap(mean_df, "0_weight", "partial_matches")
	lineplot_matches(mean_df, "full_matches", x_axis='Coverage_ReadLength', hue='Barcode') #all matches per weight: "full_matches", x_axis='Coverage_ReadLength', hue='0_weight'
	lineplot_matches(mean_df.sort_values(by="mean_read_length"), "partial_matches", x_axis='Coverage_ReadLength', hue='Barcode')
	#
	lineplot_matches_barcode(mean_df, "0_weight", "full_matches", x_axis='Coverage_ReadLength', hue='Barcode',palette=sns.color_palette("deep"))
	lineplot_matches_barcode(mean_df.sort_values(by="mean_read_length"), "0_weight", "partial_matches", x_axis='Coverage_ReadLength', hue='Barcode', palette=sns.color_palette("deep"))