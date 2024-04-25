#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

out_dir="./out/"
inputdata="./out/Homogenous_NonWeigthed_20_Iterations_SummaryTable.csv"
prefix='Homogenous_NonWeigthed_20_Iterations' #sample name for output plot

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
 

def lineplot_matches(data, value_column, out_dir):
	"""
	Plot partial and full matches over the mean read lengths with group ID being the coverage and save the plot as JPG files.
	"""
	# Create a line plot for full matches
	plt.figure()
	g = sns.lineplot(data=data, x='mean_read_length', y=value_column, hue='coverage', legend='full', errorbar='sd')
	g.set_xticks(data['mean_read_length'].unique())
	g.tick_params(axis='x', labelrotation=90)
	g.set(xlabel='Mean Read Length', ylabel=value_column)
	g.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Coverage')

	# Save the plot as a JPG file
	#add experimental line
	plt.axvline(x=5200, color='black', linestyle=':')
	plt.text(5200, plt.ylim()[1], 'MK025', verticalalignment='bottom', horizontalalignment='center')
	# Add point at y=5 on the vertical line
	plot_points(MK025_data, option=value_column, color=sns.cubehelix_palette())

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
		#	plt.plot(data['mean_read_length'].unique(), mean_data['partial_matches'][n], marker='o', color=color[n], mec='black', mew=1)
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
	plt.savefig(plot_path, bbox_inches='tight')

def lineplot_matches_barcode(data, value_column, out_dir):
	"""
	Plot partial and full matches over the mean read lengths and coverages. This splits by barcodes, so we essentially get the each barcodes individual numbers (as mean per x iterations) as a single line. This shows that there are no preferred barcodes.
	"""
	# Create a line plot for full matches
	plt.figure()
	#data["Barcode"] = data["Insertion"].str.split("_insertion").str[0]
	unique_barcodes = data["Iteration"].unique()
	for n, barcode in enumerate(unique_barcodes):
		barcode_data = data[data["Iteration"] == barcode]
		if n == 0:
			g = sns.lineplot(data=barcode_data, x='mean_read_length', y=value_column, hue='coverage', legend=True,errorbar=None) #collapse multiple iterations
		g = sns.lineplot(data=barcode_data, x='mean_read_length', y=value_column, hue='coverage', legend=False, errorbar=None) #collapsed 10
	
	g.set_xticks(data['mean_read_length'].unique())
	g.tick_params(axis='x', labelrotation=90)
	g.set(xlabel='Mean Read Length', ylabel=value_column)
	g.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Coverage')

	# add experimental line
	#max value of data column
	max_value = max(data[value_column].tolist())
	print(max_value)
	print(data.sort_values(by=[value_column]))
	plt.axvline(x=5200, color='black', linestyle=':')
	plt.text(5200,max_value, 'MK025', verticalalignment='bottom', horizontalalignment='center')
	# Add point at y=5 on the vertical line
	plot_points(MK025_data, option=value_column, color=sns.cubehelix_palette())

	plotname = prefix + str(value_column) + '_barcode_lineplot.jpg'
	full_matches_plot_path = os.path.join(out_dir, plotname)
	plt.savefig(full_matches_plot_path, bbox_inches='tight')

# Example usage:
inputdata_df=pd.read_csv(inputdata, sep='\t')
MK025_data = pd.read_csv("./out/MK025_10xSummaryTable.csv", sep='\t')

inputdata_df["Barcode"] = inputdata_df["Insertion"].str.split("_insertion").str[0]
inputdata_df["Iteration"] = inputdata_df["Insertion"].str.split("_").str[4]

#sums up the full matches and partial matches across the barcodes but keeps the iterations for the statistics int he plot!
summed_df = inputdata_df.groupby(['coverage', 'mean_read_length', 'Iteration']).agg({'full_matches': 'sum', 'partial_matches': 'sum'}).reset_index()
print(summed_df)
#plot_matches(inputdata_df, out_dir)
#plot_combined_matches(finputdata_df, out_dir)
lineplot_matches(summed_df, "full_matches", out_dir)
lineplot_matches(summed_df, "partial_matches", out_dir)
lineplot_matches_barcode(summed_df, "full_matches", out_dir)
lineplot_matches_barcode(summed_df, "partial_matches", out_dir)
#lineplot_matches(inputdata_df, "full_matches", out_dir)
#lineplot_matches(inputdata_df, "partial_matches", out_dir)

#calculates the mean value for each barcode across the iterations
print(inputdata_df)
mean_df = inputdata_df.groupby(['coverage', 'mean_read_length', 'Barcode']).agg({'full_matches': 'sum', 'partial_matches': 'sum'}).reset_index()
print(mean_df)
plot_barcode_barplot(mean_df, "full_matches", out_dir, prefix)
plot_barcode_barplot(mean_df, "partial_matches", out_dir, prefix)