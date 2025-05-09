#!/bin/python3

#%%
import pandas as pd
from scipy.stats import pearsonr, sem, ttest_1samp
import sys
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

import re
from datetime import datetime
pio.kaleido.scope.mathjax = None
pio.templates.default = "simple_white"

from src.esloco.config_handler import seq_read_data

#%%
seq_read_lengths = seq_read_data("/home/weichan/permanent/Projects/VIS/Data/VIS_Magdeburg/20240205_1448_MN35428_FAX66700_46d2ede9/MK025_GFPpos_sup_dorado_ref_simulation_1kb.fasta.gz", distribution=True, min_read_length=0)
print(len(seq_read_lengths))
#%%
'''
if len(seq_read_lengths) > 1000000:
    sampled_seq_read_lengths = np.random.choice(seq_read_lengths, size=1000000, replace=False)
else:
    sampled_seq_read_lengths = seq_read_lengths

# Subsample 1,000,000 reads if the total number of reads exceeds this limit

transformed = np.log10(sampled_seq_read_lengths)
hist, bin_edges = np.histogram(transformed, bins=100)

# Plot the histogram
plt.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), edgecolor="black", align="edge")
plt.xlabel("Read Length")
plt.ylabel("Frequency")
plt.title("Histogram of Sequence Read Lengths")
plt.show()

histograms = {"seq": transformed}
#sns.histplot(seq_read_lengths, bins=10000, kde=True)
# %%
all_histograms = []
all_bin_edges = []
all_means = []
all_ids = []
dir = "/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_bigtest/plots/"
# Process each .npy file
for file in os.listdir(dir):
    if file.endswith(".npy"):
        file_path = os.path.join(dir, file)
        data = np.load(file_path)
        print(f"{file}: Mean = {np.mean(data)}")
        print(f"{file}: n_reads = {len(data)}")
        
        # Compute histogram and bin edges
        transformed_data = np.log10(data)
        histograms[file] = transformed_data

        #hist, bin_edges = np.histogram(transformed_data, bins=100)
        #all_histograms.append(hist)
        #all_bin_edges.append(bin_edges)
        #all_means.append(np.log10(np.mean(data)))
        #all_ids.append(file.split("_")[0] + "-" + file.split("_")[1])

# Add the original histogram
all_histograms.append(hist)
all_bin_edges.append(bin_edges)
all_means.append(np.log10(np.mean(sampled_seq_read_lengths)))
all_ids.append("Sequence Read Lengths")
'''
#%%
# Plot all histograms
'''
fig, axes = plt.subplots(int(np.ceil(len(all_histograms)/2)), 2, figsize=(9, 9))
axes = axes.flatten()  # Flatten the axes array to iterate over it
for i, (hist, bin_edges, mean, id, ax) in enumerate(zip(all_histograms, all_bin_edges, all_means, all_ids, axes[:len(all_histograms)])):
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    density = hist / np.sum(hist) / np.diff(bin_edges)  # Normalize to create density
    sns.lineplot(x=bin_centers, y=density, ax=ax, label=id, alpha=0.3, linewidth=5, color=sns.color_palette("husl", len(all_histograms))[i])
    ax.axvline(mean, color=sns.color_palette("husl", len(all_histograms))[i], linestyle="--", linewidth=1)
    ax.set_title(f"Histogram for {id}")
    ax.set_ylabel("Density")
    ax.legend()

plt.tight_layout()
plt.xlabel("")
plt.ylabel("")
fig.supxlabel('Log10(Read Length)')
#plt.title("Combined Histograms of Sequence Read Lengths")
plt.legend()
plt.show()

fig, ax = plt.subplots(figsize=(9, 6))

for i, (hist, bin_edges, mean, id) in enumerate(zip(all_histograms, all_bin_edges, all_means, all_ids)):
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    density = hist / np.sum(hist) / np.diff(bin_edges)  # Normalize to create density
    sns.lineplot(x=bin_centers, y=density, ax=ax, label=id, alpha=0.5, linewidth=5, linestyle=":",
                 color=sns.color_palette("husl", len(all_histograms))[i])
    ax.axvline(mean, color=sns.color_palette("husl", len(all_histograms))[i], linestyle="--", linewidth=1, alpha=0.7)

ax.set_title("Combined Histograms of Sequence Read Lengths")
ax.set_xlabel("Log10(Read Length)")
ax.set_ylabel("Density")
ax.legend(title="Dataset")
plt.tight_layout()
'''
#%%
'''
for key in histograms.keys():
    print(f"{key}: {np.mean(histograms[key])}")
    sns.histplot(histograms[key], bins=100, kde=True, label=key)
'''

# %% plot with true positive lines
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
pio.kaleido.scope.mathjax = None
pio.templates.default = "plotly_white"

def barplot_absolute_matches(data, tp_dict, noplot=False):

    # Output paths
    #output_svg = os.path.join(output_path, f"{experiment_name}_Barplot_absolute_numbers.svg")
    #output_html = os.path.join(output_path, f"{experiment_name}_Barplot_absolute_numbers.html")

    # Extract iteration number from the 'Insertion' column
    #data['Iteration'] = data['target_region'].str.extract(r'_(\d+)$').astype(int)

    # Group by mean_read_length, coverage, and iteration # sum up all barcodes
    summary = data.groupby(['mean_read_length', 'coverage', 'iteration']).agg(
        full_matches_total=('full_matches', 'sum'),
        partial_matches_total=('partial_matches', 'sum'),
        bases_on_target_total=('on_target_bases', 'sum')
    ).reset_index()

    print(summary.head())
    # Aggregate across iterations to compute mean and standard error
    final_summary = summary.groupby(['mean_read_length', 'coverage']).agg(
        full_matches_mean=('full_matches_total', 'mean'),
        full_matches_se=('full_matches_total', lambda x: x.std() / np.sqrt(len(x))),
        partial_matches_mean=('partial_matches_total', 'mean'),
        partial_matches_se=('partial_matches_total', lambda x: x.std() / np.sqrt(len(x))),
        on_target_bases_mean=('bases_on_target_total', 'median'),
        on_target_bases_se=('bases_on_target_total', lambda x: x.std() / np.sqrt(len(x)))
    ).reset_index()

    print(final_summary)

    if noplot:
        return final_summary
    
    x_labels = [f"{row['mean_read_length']}, {row['coverage']}" for _, row in final_summary.iterrows()]
    x = np.arange(len(x_labels))  # Numerical indices for the bars

    fig = make_subplots(rows=1, cols=3, subplot_titles=["Full Matches", "Partial Matches", "OTBs"], 
                       x_title="Mean Read Length, Coverage", y_title="Matches (mean)")

    fig.add_trace(go.Bar(
        x=x,
        y=final_summary['full_matches_mean'],
        error_y=dict(type='data', array=final_summary['full_matches_se']),
        name='Full Matches',
        marker_color='black'
    ), row=1, col=1)

    # Add a horizontal line at y = dict["full"]
    fig.add_hline(
        y=tp_dict["full"],
        line_dash="dash",
        line_color="red",
        annotation_text=f"Full: {tp_dict['full']}",
        annotation_position="top left",
        row=1,
        col=1
    )

    fig.add_trace(go.Bar(
        x=x,
        y=final_summary['partial_matches_mean'],
        error_y=dict(type='data', array=final_summary['partial_matches_se']),
        name='Partial Matches',
        marker_color='grey'
    ), row=1, col=2)

    # Add a horizontal line at y = dict["full"]
    fig.add_hline(
        y=tp_dict["partial"],
        line_dash="dash",
        line_color="red",
        annotation_text=f"Partial: {tp_dict['partial']}",
        annotation_position="top left",
        row=1,
        col=2
    )

    fig.add_trace(go.Bar(
        x=x,
        y=final_summary['on_target_bases_mean'],
        error_y=dict(type='data', array=final_summary['on_target_bases_se']),
        name='OTBs',
        marker_color='lightgrey'
    ), row=1, col=3)
    # Add a horizontal line at y = dict["full"]
    fig.add_hline(
        y=tp_dict["otb"],
        line_dash="dash",
        line_color="red",
        annotation_text=f"OTB: {tp_dict['otb']}",
        annotation_position="top left",
        row=1,
        col=3
    )

    # Update x-axis labels for all subplots
    fig.update_layout(
        xaxis=dict(
            tickmode='array',
            tickvals=x,
            ticktext=x_labels
        ),
        xaxis2=dict(
            tickmode='array',
            tickvals=x,
            ticktext=x_labels
        ),
        xaxis3=dict(
            tickmode='array',
            tickvals=x,
            ticktext=x_labels
        )
    )
    fig.show()

    #fig.write_html(output_html)
    #fig.write_image(output_svg, width=600, height=400)
    #print(f"Barplot absolute numbers saved as {output_html}")
    #return final_summary
#%%
def read_data(filepath, header=0):
    data = pd.read_csv(filepath, sep='\t', header=header)
    return data

#%%
data = read_data('/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_5_VCN/1_5_10k_matches_table.csv')
print(data.head())
#%% seq data + VIS

vis = read_data('/home/weichan/temporary/Projects/VIS_out/Simulation/final/localization/ExactInsertions_one_run_28z.bed', header=None)
vis["length"] = vis[2] - vis[1]
vis["length"] = vis["length"].astype(int)
print(vis.head())

# Create the dictionary with the required keys and values
vis_dict = {
    "full": len(vis[vis["length"] > 5000]),
    "partial": len(vis[vis["length"] <= 5000]),
    "otb": vis["length"].sum()
}

print(vis_dict)
print(vis_dict["full"])
#%%

barplot_absolute_matches(data, vis_dict)

# %%
paths7 =  [(1,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/1_7_matches_table.csv'),
         (10,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/10_7_matches_table.csv'),
         (100,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/100_7_matches_table.csv'),
         (1000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/1000_7_matches_table.csv'),
         (10000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/10000_7_matches_table.csv')]

paths5 = [(1,5, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/1_5_matches_table.csv'),
         (10,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/10_5_matches_table.csv'),
         (100,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/100_5_matches_table.csv'),
         (1000,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/1000_5_matches_table.csv')]

paths10 =  [(1,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/1_10_matches_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/10_10_matches_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/100_10_matches_table.csv'),
         (1000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/1000_10_matches_table.csv'),
         (10000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/10000_10_matches_table.csv')]

paths12 = [(1000,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/12/1000_12_matches_table.csv'),
           (10000,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/12/10000_12_matches_table.csv')]

paths15 =  [(1,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/1_15_matches_table.csv'),
         (10,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/10_15_matches_table.csv'),
         (100,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/100_15_matches_table.csv'),
         (1000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/1000_15_matches_table.csv'),
         (10000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/10000_15_matches_table.csv')]

paths = [list(paths7), list(paths5), list(paths10),list(paths12), list(paths15)]
print(paths)
#%%
summaries = []
for path_group in paths:
    for n, vcn, i in path_group:
        data = read_data(i)
        summary = barplot_absolute_matches(data, vis_dict, noplot=True)
        summary["VCN"] = vcn
        summary["n"] = n
        summary["log_n"] = np.log10(n)
        summaries.append(summary)

#%%
# Combine all summaries into a single DataFrame
# Assign colors to the respective VCNs
vcn_colors = {
    5: "#0000FF",  # Blue
    7: "#008000",  # Green
    10: "#FFA500",  # Orange
    12: "#A9A500",  
    15: "#800080"  # Purple
}

combined_summary = pd.concat(summaries, ignore_index=True)

# Plot using seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(data=combined_summary, x="log_n", y="on_target_bases_mean", hue="VCN", 
             palette=vcn_colors, marker="o", linewidth=5)

# Add error bars
plt.errorbar(combined_summary["log_n"], combined_summary["on_target_bases_mean"], 
             yerr=combined_summary["on_target_bases_se"], fmt='none', c='black', capsize=3, label="SE")

# Find the intersection point
intersection_x = np.log10(15000)
intersection_y = vis_dict["otb"]

# Plot horizontal and vertical lines ending at the intersection point
plt.plot([combined_summary["log_n"].min(), intersection_x], [intersection_y, intersection_y],
          label="True OTB", linewidth=2, alpha=0.5, color='red', linestyle="--")
plt.plot([intersection_x, intersection_x], [combined_summary["on_target_bases_mean"].min(), intersection_y],
          label="Max. n", linewidth=2, alpha=0.5, color='red', linestyle="--")

# Add a label at the intersection point
plt.text(intersection_x, intersection_y, f"({intersection_x:.2f}, {intersection_y:.2f})", 
         color='red', fontsize=10, ha='left', va='bottom')

plt.title("On Target Bases Mean vs n")
plt.xlabel("log10(n)")
plt.ylabel("On Target Bases Median")
plt.legend(title="VCN")
plt.tight_layout()
plt.show()

#%%
selection = 10
clonality = []
for n, vcn, i in paths10:
    data = read_data(i)
    summary = data.groupby(['barcode','mean_read_length', 'coverage']).agg(
        bases_on_target_total=('on_target_bases', 'mean')
    ).reset_index()
    summary["VCN"] = vcn
    summary["n"] = n
    summary["log_n"] = np.log10(n)
    print(summary.head())
    clonality.append(summary)

combined_clonality = pd.concat(clonality, ignore_index=True)

#rgba colors corresponding to the hexcodes above, might need to be changed later
vcn_colors = {
    5: "rgba(0, 0, 255",  # Blue with 80% opacity
    7: "rgba(0, 128, 0",  # Green with 80% opacity
    10: "rgba(255, 165, 0",  # Orange with 80% opacity
    12: "rgba(155, 105, 0",  # Orange with 80% opacity
    15: "rgba(128, 0, 128"  # Purple with 80% opacity
}

print(combined_clonality['n'].unique())
# Assign a single color based on the VCN
vcn_color = vcn_colors[selection]
combined_clonality_sorted = combined_clonality.sort_values(by='n')
fig = px.pie(combined_clonality_sorted, values='bases_on_target_total', names='n',
             color_discrete_sequence=[vcn_color + "," + str((i+0.2)*0.2)+ ")" for i in range(len(combined_clonality['n'].unique()))])
# Update layout to remove legend and add labels directly on the sectors
fig.update_layout(title="OTBs n / Total OTBs in %", width=500, height=500, margin=dict(l=50, r=50, t=50, b=50), showlegend=False)
fig.update_traces(textinfo='label+percent', textfont_size=30, marker=dict(line=dict(color='black', width=5)))
fig.update_traces(sort=False, selector=dict(type='pie'))
fig.show()

#%%
'''
print(combined_clonality)
for n in combined_clonality["n"].unique():
    print(f"n: {n}")
    filtered = combined_clonality[combined_clonality["n"] == n]
    fig = px.pie(filtered, values='bases_on_target_total', names='barcode', labels=None)
    fig.update_layout(width=500, height=500)
    fig.update_traces(textfont_size=20)
    fig.show()
'''
#%%    
combined_clonality['percentage'] = combined_clonality.groupby('n')['bases_on_target_total'].transform(lambda x: x / x.sum() * 100)
combined_clonality['barcode'] = combined_clonality['barcode'].astype(str)  # Convert barcode to string for better labeling

# Define a custom color mapping
def custom_color_mapping(barcode):
    if barcode == '0':
        return 'red'
    else:
        return f'rgba(128, 128, 128, {0.5 + 0.9 * (int(barcode) % 10) / 10})'  # Grayscale with varying transparency

# Apply the custom color mapping
combined_clonality['color'] = combined_clonality['barcode'].apply(custom_color_mapping)

# Create a stacked bar chart with custom colors
fig_bar = px.bar(combined_clonality, x='n', y='percentage', color='barcode',
                 title="",
                 labels={'percentage': 'Percentage (%)', 'barcode': 'Barcode'},
                 color_discrete_map={row['barcode']: row['color'] for _, row in combined_clonality.iterrows()})
fig_bar.update_layout(barmode='stack', width=600, height=600, xaxis=dict(type='category'), font=dict(size=20))
fig_bar.update_xaxes(showline=True, linewidth=2, linecolor='black')
fig_bar.update_yaxes(showline=True, linewidth=2, linecolor='black')
fig_bar.show()

# %%
# plot the barcode distribution table for each configuration (-> supposed to show that reads are equally distirbuted across barcodes)

bpaths7 =  [(1,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/1_7_barcode_distribution_table.csv'),
         (10,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/10_7_barcode_distribution_table.csv'),
         (100,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/100_7_barcode_distribution_table.csv'),
          (1000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/1000_7_barcode_distribution_table.csv'),
          (10000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/7/10000_7_barcode_distribution_table.csv')]

bpaths5 = [(1,5, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/1_5_barcode_distribution_table.csv'),
         (10,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/10_5_barcode_distribution_table.csv'),
         (100,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/100_5_barcode_distribution_table.csv'),
         (1000,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/1000_5_barcode_distribution_table.csv'),
         (10000,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/5/10000_5_barcode_distribution_table.csv')
         ]

bpaths10 =  [(1,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/1_10_barcode_distribution_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/10_10_barcode_distribution_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/100_10_barcode_distribution_table.csv'),
         (1000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/1000_10_barcode_distribution_table.csv'),
         (10000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10/10000_10_barcode_distribution_table.csv')]

bpaths15 =  [(1,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/1_15_barcode_distribution_table.csv'),
         (10,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/10_15_barcode_distribution_table.csv'),
         (100,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/100_15_barcode_distribution_table.csv'),
         (1000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/1000_15_barcode_distribution_table.csv'),
         (10000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/15/10000_15_barcode_distribution_table.csv')]


selection = 10
barcode_distribution = []
for n, vcn, i in bpaths10:
    data = read_data(i)
     # Identify barcode columns based on their position before the 'coverage' column
    coverage_index = data.columns.get_loc('coverage')
    barcode_columns = data.columns[:coverage_index]
    summary = data[barcode_columns].mean().reset_index()
    summary.columns = ['barcode', 'mean']
    summary["VCN"] = vcn
    summary["n"] = n
    summary["log_n"] = np.log10(n)
    barcode_distribution.append(summary)

combined_barcode_distribution = pd.concat(barcode_distribution, ignore_index=True)

#rgba colors corresponding to the hexcodes above, might need to be changed later
vcn_colors = {
    5: "rgba(0, 0, 255",  # Blue with 80% opacity
    7: "rgba(0, 128, 0",  # Green with 80% opacity
    10: "rgba(255, 165, 0",  # Orange with 80% opacity
    15: "rgba(128, 0, 128"  # Purple with 80% opacity
}

print(combined_barcode_distribution['n'].unique())
# Assign a single color based on the VCN
vcn_color = vcn_colors[selection]
# Sort the data by 'n' to ensure names are in ascending order
combined_barcode_distribution_sorted = combined_barcode_distribution.sort_values(by='n')
fig = px.pie(combined_barcode_distribution_sorted, values='mean', names='n', 
             color_discrete_sequence=[vcn_color + "," + str((i+0.2)*0.2)+ ")" for i in range(len(combined_barcode_distribution['n'].unique()))])
# Update layout to remove legend and add labels directly on the sectors
fig.update_traces(sort=False, selector=dict(type='pie'))
fig.update_layout(
    title="Reads n / Total reads in %",
    width=500, height=500, margin=dict(l=50, r=50, t=50, b=50), showlegend=False
)
fig.update_traces(textinfo='label+percent', textfont_size=30, marker=dict(line=dict(color='black', width=5)))

fig.show()

#%%
#variability of individual combined_barcode_distribution

print(combined_barcode_distribution)
print(combined_clonality)
combined_clonality['variance_per_barcode']=combined_clonality.groupby('n')['bases_on_target_total'].transform(lambda x: np.var(x)) #population scale sd
#combined_clonality['mean_per_barcode']=combined_clonality.groupby('n')['bases_on_target_total'].transform(lambda x: np.mean(x)) #population scale sd
combined_clonality['sd_per_barcode']=combined_clonality.groupby('n')['bases_on_target_total'].transform(lambda x: (np.std(x) * len(x))/ np.mean(x)) #population scale sd
combined_clonality['sem_per_barcode']=combined_clonality.groupby('n')['bases_on_target_total'].transform(lambda x: sem(x/len(x))) #population scale sd
print(combined_clonality)


#%%

# Calculate summary statistics
summary = combined_clonality.groupby("n")['bases_on_target_total'].agg(
    mean="mean", 
    std=(lambda x: np.std(x, ddof=1)),  # Population scale standard deviation
    count="count"
).reset_index()
summary["CV"] = summary["std"] / summary["mean"]

# Plot the data using Plotly
fig = px.bar(summary, x="n", y="CV", text="CV",
             color_discrete_sequence=px.colors.sequential.Greys_r[2:7], 
             title="Coefficient of Variation per n")
fig.update_traces(texttemplate='%{text:.2f}', textposition='auto')
fig.update_layout(
    title=dict(text="Coefficient of Variation per Clonality", font=dict(size=20)),
    width=500, height=500, font=dict(size=20)
)
fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black', title="n", showticklabels=True, ticks="")
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', title="CV", showticklabels=True, ticks="")
fig.update_yaxes(tickvals=[0, 1, 2, 3, 4])

fig.show()

#%%
print(combined_clonality.head())
fig = px.box(combined_clonality, x='n', y='bases_on_target_total', points="all", facet_col="n",
              color="n", color_discrete_sequence=px.colors.sequential.Greys_r[0:5])
fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black', tickmode='array', tickvals=combined_clonality['n'].unique())
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', showticklabels=True, matches=None)  # Set dtick=1 to keep only full number y-axis values
fig.update_layout(width=1200, height=300, font=dict(size=16), showlegend=False, 
                  title=dict(text="Individual barcode contribution to total OTBs", font=dict(size=20)))
fig.show()

#%%
# contirbutors
combined_clonality['percentage'] = combined_clonality.groupby(['n'])['bases_on_target_total'].transform(lambda x: x / x.sum() * 100)
combined_clonality['zero_contributor'] = combined_clonality['bases_on_target_total'] == 0

# Aggregate data to count zero contributors and total barcodes per 'n'
zero_contributor_summary = combined_clonality.groupby(['n']).agg(
    zero_contributor_count=('zero_contributor', 'sum'),
    total_barcodes=('barcode', 'count')
).reset_index()

# Calculate the percentage of zero contributors
zero_contributor_summary['zero_contributor_percentage'] = (zero_contributor_summary['zero_contributor_count'] / zero_contributor_summary['total_barcodes']) * 100
zero_contributor_summary['non_zero_contributor_percentage'] = 100 - zero_contributor_summary['zero_contributor_percentage']

print(zero_contributor_summary)

# Create a pie chart for each combination of n and Dominance

for _, row in zero_contributor_summary.iterrows():
    n = row['n']
    zero_percentage = row['zero_contributor_percentage']
    non_zero_percentage = row['non_zero_contributor_percentage']
    
    fig = px.pie(
        values=[zero_percentage, non_zero_percentage],
        names=["No Contributor", "Contributor"],
        title=f"n={n}",
        color_discrete_sequence=['red', 'green']
    )
    fig.update_traces(sort=False, selector=dict(type='pie'))
    fig.update_traces(textinfo='percent', textfont_size=16, marker=dict(line=dict(color='black', width=4)))
    fig.update_layout(width=400, height=400, showlegend=True)
    fig.show()

# %%
##############
############
##########
########
#######

#Dominance dilution experiment
dd10path =  [#(5,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/1_5_dd10_matches_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10_10_dd10_matches_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/100_10_dd10_matches_table.csv')]

dd1path =  [#(10,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10_10_dd1_matches_table.csv'),
            (100,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/100_10_dd1_matches_table.csv'),
            (1000,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/1000_10_dd1_matches_table.csv')]

dd01path =  [#(10,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10_10_dd1_matches_table.csv'),
            (1000,0.1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/1000_10_dd01_matches_table.csv'),
            (10000,0.1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10000_10_dd01_matches_table.csv')]


paths = [list(dd10path), list(dd1path), list(dd01path)]
print(paths)
#%%
summaries = []
for path_group in paths:
    for n, dominance, i in path_group:
        data = read_data(i)
        summary = barplot_absolute_matches(data, vis_dict, noplot=True)
        summary = data.groupby(['barcode','mean_read_length', 'coverage']).agg(
        bases_on_target_total=('on_target_bases', 'mean')
    ).reset_index()
        summary["Dominance"] = dominance
        summary["n"] = n
        summary["log_n"] = np.log10(n)
        summaries.append(summary)

#%%
# Combine all summaries into a single DataFrame
# Assign colors to the respective VCNs
vcn_colors = {
    5: "#0000FF",  # Blue
    7: "#008000",  # Green
    10: "#FFA500",  # Orange
    12: "#A9A500",  
    15: "#800080"  # Purple
}

combined_dominance = pd.concat(summaries, ignore_index=True)
print(combined_dominance.head())
# %%
combined_dominance['percentage'] = combined_dominance.groupby(['n', 'Dominance'])['bases_on_target_total'].transform(lambda x: x / x.sum() * 100)
combined_dominance['barcode'] = combined_dominance['barcode'].astype(str)  # Convert barcode to string for better labeling

# Define a custom color mapping
def custom_color_mapping(barcode):
    if barcode == '0':
        return 'orange'
    else:
        return f'rgba(128, 128, 128, {0.5 + 0.9 * (int(barcode) % 10) / 10})'  # Grayscale with varying transparency

# Apply the custom color mapping
combined_dominance['color'] = combined_dominance['barcode'].apply(custom_color_mapping)

# Create a stacked bar chart with custom colors and facet on VCN
fig_bar = px.bar(combined_dominance, x='n', y='percentage', color='barcode', facet_col='Dominance',
                 title="",
                 labels={'percentage': 'Percentage (%)', 'barcode': 'Barcode'},
                 color_discrete_map={row['barcode']: row['color'] for _, row in combined_dominance.iterrows()})
fig_bar.update_layout(barmode='stack', width=1200, height=500, font=dict(size=20), showlegend=False)
fig_bar.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
fig_bar.update_yaxes(showline=True, linewidth=2, linecolor='black')
fig_bar.update_xaxes(showline=True, linewidth=2, linecolor='black')
fig_bar.update_yaxes(showline=True, linewidth=2, linecolor='black')
fig_bar.show()

#%% CV of dominance 
summary = combined_dominance.groupby(["n", "Dominance"])['bases_on_target_total'].agg(
    mean="mean", 
    std=(lambda x: np.std(x, ddof=1)),  # Population scale standard deviation
    count="count"
).reset_index()
summary["CV"] = summary["std"] / summary["mean"]
print(summary)
# Plot the data using Plotly
fig = px.bar(summary, x="n", y="CV", text="CV", facet_col="Dominance",
             color_discrete_sequence=px.colors.sequential.Greys_r[2:7])
fig.update_traces(texttemplate='%{text:.2f}', textposition='auto')
fig.update_layout(
    title=dict(text="Coefficient of Variation per Clonality", font=dict(size=20)),
    width=1200, height=500, font=dict(size=20)
)

fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black', title="n", showticklabels=True, ticks="")
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', showticklabels=True, ticks="", matches='y')
fig.update_yaxes(tickvals=[0, 1, 2, 3, 4])
fig.show()
#%%
# separate boxplots for eahc dominance
dd10 = combined_dominance[combined_dominance["Dominance"]==10]
dd1 = combined_dominance[combined_dominance["Dominance"]==1]
dd01 = combined_dominance[combined_dominance["Dominance"]==0.1]

print(combined_dominance[combined_dominance["barcode"]=="0"]['n'])

fig = px.box(dd10, x='n', y='percentage', color="n", color_discrete_sequence=px.colors.sequential.Greys_r[2:7],
             labels={'percentage': 'Percentage (%)'})
fig.add_trace(go.Scatter(
    y=dd10[dd10["barcode"]=="0"]['percentage'],
    x=dd10[dd10["barcode"]=="0"]['n'],
    mode='markers',
    marker=dict(color="orange", size=15, line=dict(color='black', width=2)),
    showlegend=False
))
fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(showline=True, linewidth=2, linecolor='black') 
fig.update_layout(width=500, height=500, font=dict(size=20), showlegend=False)
fig.show()


fig = px.box(dd1, x='n', y='percentage', color="n", color_discrete_sequence=px.colors.sequential.Greys_r[2:7],
             labels={'percentage': 'Percentage (%)'})
fig.add_trace(go.Scatter(
    y=dd1[dd1["barcode"]=="0"]['percentage'],
    x=dd1[dd1["barcode"]=="0"]['n'],
    mode='markers',
    marker=dict(color="orange", size=15, line=dict(color='black', width=2)),
    showlegend=False
))
fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(showline=True, linewidth=2, linecolor='black') 
fig.update_layout(width=500, height=500, font=dict(size=20), showlegend=False)
fig.show()

fig = px.box(dd01, x='n', y='percentage', color="n", color_discrete_sequence=px.colors.sequential.Greys_r[2:7],
             labels={'percentage': 'Percentage (%)'})
fig.add_trace(go.Scatter(
    y=dd01[dd01["barcode"]=="0"]['percentage'],
    x=dd01[dd01["barcode"]=="0"]['n'],
    mode='markers',
    marker=dict(color="orange", size=15, line=dict(color='black', width=2)),
    showlegend=False
))
fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(showline=True, linewidth=2, linecolor='black') 
fig.update_layout(width=500, height=500, font=dict(size=20), showlegend=False)
fig.show()
# %% # facet box plots

fig = px.box(combined_dominance, x='n', y='percentage', color="n", color_discrete_sequence=px.colors.sequential.Greys_r[2:7],
             labels={'percentage': 'Percentage (%)'}, facet_col="Dominance")

# Add orange dots for each facet
for dominance_value in combined_dominance["Dominance"].unique():
    filtered_data = combined_dominance[(combined_dominance["barcode"] == "0") & (combined_dominance["Dominance"] == dominance_value)]
    fig.add_trace(go.Scatter(
        y=filtered_data['percentage'],
        x=filtered_data['n'],
        mode='markers',
        marker=dict(color="orange", size=10, line=dict(color='black', width=2)),
        showlegend=False
    ), row=1, col=list(combined_dominance["Dominance"].unique()).index(dominance_value) + 1)

fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(showline=True, linewidth=2, linecolor='black') 
fig.update_layout(width=1200, height=500, font=dict(size=20), showlegend=False,
                  title=dict(text="OTBs per n", font=dict(size=20)))
fig.show()

# %%
#pie chart with number of n contributing 0
combined_dominance['percentage'] = combined_dominance.groupby(['n', 'Dominance'])['bases_on_target_total'].transform(lambda x: x / x.sum() * 100)
combined_dominance['zero_contributor'] = combined_dominance['bases_on_target_total'] == 0

# Aggregate data to count zero contributors and total barcodes per 'n' and 'Dominance'
zero_contributor_summary = combined_dominance.groupby(['n', 'Dominance']).agg(
    zero_contributor_count=('zero_contributor', 'sum'),
    total_barcodes=('barcode', 'count')
).reset_index()

# Calculate the percentage of zero contributors
zero_contributor_summary['zero_contributor_percentage'] = (zero_contributor_summary['zero_contributor_count'] / zero_contributor_summary['total_barcodes']) * 100
zero_contributor_summary['non_zero_contributor_percentage'] = 100 - zero_contributor_summary['zero_contributor_percentage']

print(zero_contributor_summary)
#%%
# Create a pie chart for each combination of n and Dominance
for _, row in zero_contributor_summary.iterrows():
    n = row['n']
    dominance = row['Dominance']
    zero_percentage = row['zero_contributor_percentage']
    non_zero_percentage = row['non_zero_contributor_percentage']
    
    fig = px.pie(
        values=[zero_percentage, non_zero_percentage],
        names=["No Contributor", "Contributor"],
        title=f"n={n}, Dominance={dominance}",
        color_discrete_sequence=['red', 'green']
    )
    fig.update_traces(sort=False, selector=dict(type='pie'))
    fig.update_traces(textinfo='percent', textfont_size=16, marker=dict(line=dict(color='black', width=4)))
    fig.update_layout(width=400, height=400, showlegend=True)
    fig.show()
#%%
# %%


#####
#####
#####
sample=dd10
#Dominance significance test
# dominant clone higher than background
sample_test = {}
for n in sample["n"].unique():
    print(f"n: {n}")
    filtered = sample[sample["n"] == n]
    background=filtered[filtered["barcode"]!="0"]['bases_on_target_total']
    dominant =filtered[filtered["barcode"]=="0"]['bases_on_target_total']
    print(f"background: {background.mean()}, dominant: {dominant.mean()}")
    stat, p_value = ttest_1samp(background, dominant.mean(), alternative='less')
    print(f"t-statistic: {stat}, p-value: {p_value}")
    # Rare event detection via non-parametric empirical p values
    p_empirical = np.mean(background >= dominant.mean()) #how many extreme values are larger the dominant value
    print(f"empirical p-value: {p_empirical}")
    sample_test[n]=[p_value, p_empirical]
# %%

index = ["One-sided t-test", "Empirical test"]
df = pd.DataFrame(sample_test, index=index)

# Add significance labels
def add_significance_label(value):
    if value > 0.05:
        return "n.s."
    elif value <= 0.001:
        return "***"
    elif value <= 0.01:
        return "**"
    elif value <= 0.05:
        return "*"

# Apply the function to create annotated labels
annotated_labels = df.map(lambda x: f"{x:.3f} \n {add_significance_label(x)}").values
print(annotated_labels)
fig = px.imshow(
    df,
    color_continuous_scale=["green", "white"],
    labels=dict(x="Clonality", y="Test Type"),
    text_auto=True
)
fig.update_traces(text=annotated_labels, texttemplate="%{text}")
fig.update_coloraxes(cmin=0, cmax=0.05)  # Set the color scale range to threshold at 0.05

fig.update_layout(
    title="",
    xaxis_title="Clonality (n)",
    yaxis_title="p-value",
    font=dict(size=20),
)
fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
fig.show()

# %%
