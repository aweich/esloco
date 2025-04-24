#!/bin/python3

#%%
import pandas as pd
from scipy.stats import pearsonr
import sys
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

from config_handler import seq_read_data

#%%
seq_read_lengths = seq_read_data("/home/weichan/permanent/Projects/VIS/Data/VIS_Magdeburg/20240205_1448_MN35428_FAX66700_46d2ede9/MK025_GFPpos_sup_dorado_ref_simulation_1kb.fasta.gz", distribution=True, min_read_length=0)
print(len(seq_read_lengths))
#%%

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

#%%
# Plot all histograms
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

#%%
for key in histograms.keys():
    print(f"{key}: {np.mean(histograms[key])}")
    sns.histplot(histograms[key], bins=100, kde=True, label=key)


# %% plot with true positive lines
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
pio.kaleido.scope.mathjax = None
pio.templates.default = "plotly_white"

def barplot_absolute_matches(data, tp_dict):

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
    return final_summary
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
paths = [(1,5, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/1_genome_matches_table.csv'),
         (10,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10_genome_matches_table.csv'),
         (100,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10_genome_matches_table.csv'),
         (1,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/1_10_genome_matches_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10_10_genome_matches_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/100_10_genome_matches_table.csv'),
         (1,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/1_7_genome_matches_table.csv'),
         (10,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10_7_genome_matches_table.csv'),
         (100,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/100_7_genome_matches_table.csv')]#,
         #(10000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/10000_7_genome_matches_table.csv')]
summaries = []
for n, vcn, i in paths:
    data = read_data(i)
    summary = barplot_absolute_matches(data, vis_dict)
    summary["VCN"] = vcn
    summary["n"] = n
    summaries.append(summary)

#%%
# Combine all summaries into a single DataFrame
combined_summary = pd.concat(summaries, ignore_index=True)

# Plot using seaborn
plt.figure(figsize=(8, 6))
sns.scatterplot(data=combined_summary, x="n", y="on_target_bases_mean", hue="VCN", palette="viridis", s=100)
plt.errorbar(combined_summary["n"], combined_summary["on_target_bases_mean"], 
             yerr=combined_summary["on_target_bases_se"], fmt='none', c='black', capsize=3, label="SE")
plt.axhline(y=vis_dict["otb"], color='red', linestyle='-', label="True OTB")
plt.title("On Target Bases Mean vs n")
plt.xlabel("n")
plt.ylabel("On Target Bases Median")
plt.legend(title="VCN")
plt.tight_layout()
plt.show()

# %%
