#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os

out_dir="./out/"
finaloutput="./out/10xSummaryTable.csv"
prefix='Added_MK025_10x'

def plot_matches(data, out_dir):
    """
    Plot partial and full matches over the mean read lengths with group ID being the coverage and save the plot as a JPG file.

    Parameters:
        data (pandas.DataFrame): DataFrame containing the data with columns 'mean_read_length', 'coverage', 'insertion_number', 'full_matches', and 'partial_matches'.
        out_dir (str): Output directory where the JPG file will be saved.
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
 

def lineplot_matches(data, out_dir):
    """
    Plot partial and full matches over the mean read lengths with group ID being the coverage and save the plot as JPG files.
    """
    # Create a line plot for full matches
    plt.figure()
    g = sns.lineplot(data=data, x='mean_read_length', y='full_matches', hue='coverage', legend='full', errorbar='sd')
    g.set_xticks(data['mean_read_length'].unique())
    g.set_yticks((0,5,25,50,75,100)) #custom ticks
    g.set_yticklabels(['0','5','25','50','75','100']) 
    g.tick_params(axis='x', labelrotation=90)
    g.set(xlabel='Mean Read Length', ylabel='Full Matches')
    g.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Coverage')

    # Save the plot as a JPG file
    #add experimental line
    plt.axvline(x=5200, color='black', linestyle=':')
    plt.text(5200, plt.ylim()[1], 'MK025', verticalalignment='bottom', horizontalalignment='center')
    # Add point at y=5 on the vertical line
    plot_points(MK025_data, option="full", color=sns.cubehelix_palette())

    plotname = prefix + 'line_full_matches_plot.jpg'
    full_matches_plot_path = os.path.join(out_dir, plotname)
    plt.savefig(full_matches_plot_path, bbox_inches='tight')

    # Create a line plot for partial matches
    plt.figure()
    g = sns.lineplot(data=data, x='mean_read_length', y='partial_matches', hue='coverage', legend='full', errorbar='sd')
    g.set_xticks(data['mean_read_length'].unique())
    g.tick_params(axis='x', labelrotation=90)
    g.set(xlabel='Mean Read Length', ylabel='Partial Matches')
    g.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Coverage')

    # Save the plot as a JPG file
    plt.axvline(x=5200, color='black', linestyle=':')
    plt.text(5200, plt.ylim()[1], 'MK025', verticalalignment='bottom', horizontalalignment='center')
    # Add point at y=5 on the vertical line
    
    plot_points(MK025_data, option="partial", color=sns.cubehelix_palette())
    plotname = prefix + 'line_partial_matches_plot.jpg'
    partial_matches_plot_path = os.path.join(out_dir, plotname)
    plt.savefig(partial_matches_plot_path,bbox_inches='tight')

def plot_points(data, option, color):
    
    mean_data = data.groupby('coverage')[['full_matches', 'partial_matches']].mean().reset_index()
    print(mean_data.head())
    if option == "partial":
        #experimental value
        plt.plot(5200, 22, marker='o', color='red', mec='black', mew=1)
        #plt.text(5200 + 500, 22, 'Observed value', verticalalignment='bottom', horizontalalignment='left') #+10 for space
        for n,i in enumerate(mean_data['coverage']):
            plt.plot(data['mean_read_length'].unique(), mean_data['partial_matches'][n], marker='o', color=color[n], mec='black', mew=1)
    elif option == "full":
        plt.plot(5200, 5, marker='o', color='red', mec='black', mew=1)
        plt.text(5200 + 500, 5, 'Observed value', verticalalignment='bottom', horizontalalignment='left') #+10 for space
        for n,i in enumerate(mean_data['coverage']):
            plt.plot(data['mean_read_length'].unique(), mean_data['full_matches'][n],  marker='o', color=color[n], mec='black', mew=1)
    else:
        print("no valid option")

# Example usage:
# Assuming 'finaloutput_df' is your DataFrame containing the data
# Assuming 'out_dir' is the output directory

# Example usage:
finaloutput_df=pd.read_csv(finaloutput, sep='\t')
MK025_data = pd.read_csv("./out/MK025_10xSummaryTable.csv", sep='\t')
#plot_matches(finaloutput_df, out_dir)
#plot_combined_matches(finaloutput_df, out_dir)
lineplot_matches(finaloutput_df, out_dir)