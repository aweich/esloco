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
        on_target_bases_mean=('bases_on_target_total', 'mean'),
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
paths7 =  [(1,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/1_7_matches_table.csv'),
         (10,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/10_7_matches_table.csv'),
         (100,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/100_7_matches_table.csv'),
         (1000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/1000_7_matches_table.csv'),
         (10000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/10000_7_matches_table.csv')]

paths5 = [(1,5, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/1_5_matches_table.csv'),
         (10,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/10_5_matches_table.csv'),
         (100,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/100_5_matches_table.csv'),
         (1000,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/1000_5_matches_table.csv'),
         (10000,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/10000_5_matches_table.csv')]

paths10 =  [(1,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/1_10_matches_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/10_10_matches_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/100_10_matches_table.csv'),
         (1000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/1000_10_matches_table.csv'),
         (10000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/10000_10_matches_table.csv')]

paths12 = [(1,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/1_12_matches_table.csv'),
    (10,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/10_12_matches_table.csv'),
    (100,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/100_12_matches_table.csv'),
    (1000,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/1000_12_matches_table.csv'),
           (10000,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/10000_12_matches_table.csv')]

paths15 =  [(1,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/1_15_matches_table.csv'),
         (10,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/10_15_matches_table.csv'),
         (100,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/100_15_matches_table.csv'),
         (1000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/1000_15_matches_table.csv'),
         (10000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/10000_15_matches_table.csv')]

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
    5: "#8BE6E9",  # Blue
    7: "#aeeb9b",  # Green
    10: "#FFA500",  # Orange
    12: "#DB5050",  
    15: "#e697df"  # Purple
}

combined_summary = pd.concat(summaries, ignore_index=True)

# Plot using seaborn
plt.figure(figsize=(6, 7))
sns.lineplot(data=combined_summary, x="log_n", y="on_target_bases_mean", hue="VCN", 
             palette=vcn_colors, marker="o", linewidth=5)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

# Add error bars
plt.errorbar(combined_summary["log_n"], combined_summary["on_target_bases_mean"], 
             yerr=combined_summary["on_target_bases_se"], fmt='none', c='black', capsize=3, label="SE")

# Find the intersection point
#intersection_x = np.log10(15000)
#intersection_y = vis_dict["otb"]

# Plot horizontal and vertical lines ending at the intersection point
#plt.plot([combined_summary["log_n"].min(), intersection_x], [intersection_y, intersection_y],
#          label="True OTB", linewidth=2, alpha=0.5, color='red', linestyle="--")
#plt.plot([intersection_x, intersection_x], [combined_summary["on_target_bases_mean"].min(), intersection_y],
#          label="Max. n", linewidth=2, alpha=0.5, color='red', linestyle="--")

# Add a label at the intersection point
#plt.text(intersection_x, intersection_y, f"({intersection_x:.2f}, {intersection_y:.2f})", 
 #        color='red', fontsize=10, ha='left', va='bottom')
plt.xticks(ticks=combined_summary["log_n"].unique())#, labels=[f"{int(10**x)}" for x in combined_summary["log_n"].unique()])
plt.title("")
plt.xlabel("log10(n)")
plt.ylabel("Mean  OTBs")
plt.legend(title="VCN", loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=len(vcn_colors)+1, fontsize=10, title_fontsize=11)
#plt.savefig(f"../ProgressReport/simulation_description/plots/lineplot_fig3.svg", format='svg', bbox_inches='tight')
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
        return "#FFA500"
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
#fig_bar.write_image(f"../ProgressReport/simulation_description/plots/barplot_fig3.png")
fig_bar.show()

# %%
# plot the barcode distribution table for each configuration (-> supposed to show that reads are equally distirbuted across barcodes)

bpaths7 =  [(1,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/1_7_barcode_distribution_table.csv'),
         (10,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/10_7_barcode_distribution_table.csv'),
         (100,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/100_7_barcode_distribution_table.csv'),
          (1000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/1000_7_barcode_distribution_table.csv'),
          (10000,7, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/7/10000_7_barcode_distribution_table.csv')
          ]

bpaths5 = [(1,5, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/1_5_barcode_distribution_table.csv'),
         (10,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/10_5_barcode_distribution_table.csv'),
         (100,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/100_5_barcode_distribution_table.csv'),
         (1000,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/1000_5_barcode_distribution_table.csv'),
         (10000,5,'/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/5/10000_5_barcode_distribution_table.csv')
         ]

bpaths10 =  [(1,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/1_10_barcode_distribution_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/10_10_barcode_distribution_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/100_10_barcode_distribution_table.csv'),
         (1000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/1000_10_barcode_distribution_table.csv'),
         (10000,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/10/10000_10_barcode_distribution_table.csv')
         ]

bpaths12 =  [(1,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/1_12_barcode_distribution_table.csv'),
         (10,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/10_12_barcode_distribution_table.csv'),
         (100,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/100_12_barcode_distribution_table.csv'),
         (1000,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/1000_12_barcode_distribution_table.csv'),
         (10000,12, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/12/10000_12_barcode_distribution_table.csv')
         ]

bpaths15 =  [(1,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/1_15_barcode_distribution_table.csv'),
         (10,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/10_15_barcode_distribution_table.csv'),
         (100,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/100_15_barcode_distribution_table.csv'),
         (1000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/1000_15_barcode_distribution_table.csv'),
         (10000,15, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example1/15/10000_15_barcode_distribution_table.csv')
         ]


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
print(combined_barcode_distribution.head())
combined_barcode_distribution['log_mean'] = np.log10(combined_barcode_distribution['mean'])

fig = px.box(combined_barcode_distribution, x='n', y='mean', points="all", facet_col="n",
             color="n", color_discrete_sequence=["black"])
fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black', tickmode='array', tickvals=combined_barcode_distribution['n'].unique())
fig.update_yaxes(showline=True, showgrid=False,  linewidth=2, linecolor='black', showticklabels=True, matches=None, tickformat=".0f")  # Limit y-axis values to two decimal places
fig.update_layout(width=1200, height=300, font=dict(size=16), showlegend=False, 
                  title=dict(text="Mean Reads from Individual Barcodes", font=dict(size=20)))
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
             color_discrete_sequence=px.colors.sequential.Greys_r[2:10], 
             title="Coefficient of Variation per n")
fig.update_traces(texttemplate='%{text:.2f}', textposition='auto')
fig.update_layout(
    title=dict(text="", font=dict(size=20)),
    width=500, height=500, font=dict(size=20)
)
fig.update_xaxes(type='category', showline=True, linewidth=3, linecolor='black', title="n", showticklabels=True, ticks="")
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', title="Coefficient of Variation", showticklabels=True, ticks="")
fig.update_yaxes(tickvals=[0, 1, 2, 3, 4])
#fig.write_image(f"../ProgressReport/simulation_description/plots/cv.svg")
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
combined_clonality['zero_contributor'] = combined_clonality['bases_on_target_total'] != 0

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
# Facet all pies in one row using plotly.subplots

fig = make_subplots(
    rows=1, cols=len(zero_contributor_summary),
    specs=[[{'type': 'domain'}] * len(zero_contributor_summary)],
    subplot_titles=[f"n={row['n']}" for _, row in zero_contributor_summary.iterrows()]
)

for idx, (_, row) in enumerate(zero_contributor_summary.iterrows()):
    zero_percentage = row['zero_contributor_percentage']
    non_zero_percentage = row['non_zero_contributor_percentage']
    fig.add_trace(
        go.Pie(
            values=[zero_percentage, non_zero_percentage],
            labels=["Barcode without OTBs", "Barcode with OTBs"],
            marker=dict(colors=['dimgrey', 'lightgrey'], line=dict(color='black', width=2.5)),
            textinfo='percent',
            sort=False,
            showlegend=True
        ),
        row=1, col=idx+1
    )

fig.update_layout(
    width=300 * len(zero_contributor_summary),
    height=300,
    title_text="",
    font=dict(size=14)
)
#fig.write_image(f"../ProgressReport/simulation_description/plots/piecharts.svg")
fig.show()


# %%
# Calculate the percentage of zero and non-zero contributors for each 'n'
zero_contributor_summary = combined_clonality.groupby(['n', 'zero_contributor']).size().reset_index(name='count')
zero_contributor_summary['percentage'] = zero_contributor_summary.groupby('n')['count'].transform(lambda x: x / x.sum() * 100)

print(zero_contributor_summary)
#zero_contributor_summary["zero_contributor"] = zero_contributor_summary["zero_contributor"].map({True: "Contributing Barcode", False: "Not contributing Barcode"})
# Create a stacked bar plot
fig_bar = px.bar(zero_contributor_summary, x='n', y='percentage', color='zero_contributor',
                 title="",
                 labels={'percentage': 'Percentage (%)', 'zero_contributor': 'Contribution to total OTBs'},
                 color_discrete_map={True: 'lightgrey', False: 'dimgrey'})

fig_bar.update_layout(barmode='stack', width=800, height=300, xaxis=dict(type='category'), font=dict(size=18))
fig_bar.update_xaxes(showline=True, linewidth=2, linecolor='black', title="n")
fig_bar.update_yaxes(showline=True, linewidth=2, linecolor='black', title="Percentage (%)", range=[0, 100])
#fig_bar.write_image(f"../ProgressReport/simulation_description/plots/barplot_NonContributors_fig3.svg")
fig_bar.show()
# %%