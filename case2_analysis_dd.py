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

def read_data(filepath, header=0):
    data = pd.read_csv(filepath, sep='\t', header=header)
    return data

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

#%%
#Dominance dilution experiment
dd10path =  [#(5,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/1_5_dd10_matches_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/10_10_dd10_matches_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/100_10_dd10_matches_table.csv')]

dd1path =  [#(10,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10_10_dd1_matches_table.csv'),
            (100,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/100_10_dd1_matches_table.csv'),
            (1000,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/1000_10_dd1_matches_table.csv')]

dd01path =  [#(10,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10_10_dd1_matches_table.csv'),
            (1000,0.1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/1000_10_dd01_matches_table.csv'),
            (10000,0.1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/10000_10_dd01_matches_table.csv')]


paths = [list(dd10path), list(dd1path), list(dd01path)]
print(paths)
#%%
summaries = []
for path_group in paths:
    for n, dominance, i in path_group:
        data = read_data(i)
        vis_dict={}
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

# %%
# %%
## create the barplot and add labels to the orange marked dot; the labels should be: 
# 1.) one sided t test p value, 
# 2.) empirical p value,
# 3.) rank of the dominant clone
def significance_test(sample, value_column='bases_on_target_total'):
    sample_test = {}
    for n in sample["n"].unique():
        print(f"n: {n}")
        filtered = sample[sample["n"] == n]
        background = filtered[filtered["barcode"] != "0"][value_column]
        dominant = filtered[filtered["barcode"] == "0"][value_column]
        print(f"background: {background.mean()}, dominant: {dominant.mean()}")
        stat, p_value = ttest_1samp(background, dominant.mean(), alternative='less')
        print(f"t-statistic: {stat}, p-value: {p_value}")
        # Rare event detection via non-parametric empirical p values
        p_empirical = np.mean(background >= dominant.mean())  # how many extreme values are larger the dominant value
        print(f"empirical p-value: {p_empirical}")
        # Test if dominant mean is higher than all background values and get ranking
        is_highest = dominant.mean() > background.max()
        rank = int((background > dominant.mean()).sum() + 1)  # 1-based rank (1 = highest)
        print(f"Dominant mean higher than all background: {is_highest}")
        total = len(background) + 1  # total number of samples including the dominant clone
        print(f"Dominant mean rank among all (1=highest): {rank} of {total}")
        sample_test[n] = [p_value, p_empirical, rank, total]
    return sample_test


# Prepare test results for annotation
def get_annotation_text(n, sample_test):
    if n in sample_test:
        pval, emp, rank, total = sample_test[n]
        return f"p-value t-test: {pval:.2g}<br>empirical p-value: {emp:.2g}<br>rank: {rank} of {total}"
    else:
        return ""

dd10 = combined_dominance[combined_dominance["Dominance"]==10]
dd1 = combined_dominance[combined_dominance["Dominance"]==1]
dd01 = combined_dominance[combined_dominance["Dominance"]==0.1]

print(combined_dominance[combined_dominance["barcode"]=="0"]['n'])

dilutions = [dd10, dd1, dd01]
import random
for dilution in dilutions:
    fig = px.box(dilution,  x='n', y='percentage', color="n", color_discrete_sequence=px.colors.sequential.Greys_r[2:7],
                labels={'percentage': 'Percentage (%)'}, title="Contribution to total OTBs (per Barcode)")
    # Calculate significance test results for this dilution
    sample_test = significance_test(dilution)
    # Add orange marker for dominant clone with annotation
    for n in dilution["n"].unique():
        dominant = dilution[(dilution["barcode"] == "0") & (dilution["n"] == n)]
        if not dominant.empty:
            annotation = get_annotation_text(n, sample_test)
            # Add the orange marker
            fig.add_trace(go.Scatter(
                y=dominant['percentage'],
                x=dominant['n'],
                mode='markers',
                marker=dict(color="orange", size=15, line=dict(color='black', width=2)),
                showlegend=False
            ))
            # Add annotation next to the plot with an arrow

            n_index = list(sorted(dilution["n"].unique())).index(n)
            y_offset = np.mean(dominant['percentage'].values) - n_index * np.mean(dominant['percentage'].values)

            fig.add_annotation(
                x=dominant['n'],
                y=dominant['percentage'].values[0] + y_offset,
                text=annotation,
                ax=50,  # Adjust the x offset for the annotation
                ay=0,
                font=dict(size=14, color="black"),
                bgcolor="white",
                bordercolor="black",
                borderwidth=1,
                showarrow=True
            )
    fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black') 
    fig.update_layout(width=500, height=500, font=dict(size=20), showlegend=False)
    #fig.show()
    fig.write_image(f"../ProgressReport/simulation_description/plots/boxplot_OTB_dilution_{dilution['Dominance'].iloc[0]}.svg")
# %%
# plots of reads per barcode for each dilution vs OTB per barcode

bdd10path =  [#(5,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/1_5_dd10_matches_table.csv'),
         (10,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/10_10_dd10_barcode_distribution_table.csv'),
         (100,10, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/100_10_dd10_barcode_distribution_table.csv')]

bdd1path =  [#(10,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10_10_dd1_matches_table.csv'),
            (100,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/100_10_dd1_barcode_distribution_table.csv'),
            (1000,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/1000_10_dd1_barcode_distribution_table.csv')]

bdd01path =  [#(10,1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_calc_VCN/dd/10/10_10_dd1_matches_table.csv'),
            (1000,0.1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/1000_10_dd01_barcode_distribution_table.csv'),
            (10000,0.1, '/home/weichan/temporary/Data/Simulation/I_CAR_test/Case2_example2/dd/10/10000_10_dd01_barcode_distribution_table.csv')]


barcode_distribution = []
for path_group in [bdd10path, bdd1path, bdd01path]:
    for n, dominance, i in path_group: 
        data = read_data(i)
        # Identify barcode columns based on their position before the 'coverage' column
        coverage_index = data.columns.get_loc('coverage')
        barcode_columns = data.columns[:coverage_index]
        summary = data[barcode_columns].mean().reset_index()
        summary.columns = ['barcode', 'mean']
        summary["Dominance"] = dominance
        summary["n"] = n
        summary["log_n"] = np.log10(n)
        barcode_distribution.append(summary)

combined_barcode_distribution = pd.concat(barcode_distribution, ignore_index=True)

combined_barcode_distribution['percentage'] = combined_barcode_distribution.groupby(['n'])['mean'].transform(lambda x: x / x.sum() * 100)
print(combined_barcode_distribution.head())


#%%
# Define a custom color mapping
def custom_color_mapping(barcode):
    if barcode == '0':
        return 'red'
    else:
        return f'rgba(128, 128, 128, {0.5 + 0.9 * (int(barcode) % 10) / 10})'  # Grayscale with varying transparency

# Apply the custom color mapping
combined_barcode_distribution['color'] = combined_barcode_distribution['barcode'].apply(custom_color_mapping)

# Create a stacked bar chart with custom colors
fig_bar = px.bar(combined_barcode_distribution, x='n', y='percentage', color='barcode',
                 title="Mean Reads per Barcode: 0.1% Dilution",
                 labels={'percentage': 'Percentage (%)', 'barcode': 'Barcode'},
                 color_discrete_map={row['barcode']: row['color'] for _, row in combined_barcode_distribution.iterrows()})
fig_bar.update_layout(barmode='stack', width=600, height=600, xaxis=dict(type='category'), font=dict(size=20))
fig_bar.update_xaxes(showline=True, linewidth=2, linecolor='black')
fig_bar.update_yaxes(showline=True, linewidth=2, linecolor='black')
fig_bar.show()
# %%
# repeat the boxplot with the labels for the reads (also with the significane tests!)
print(combined_barcode_distribution.head())
dd10 = combined_barcode_distribution[combined_barcode_distribution["Dominance"]==10]
dd1 = combined_barcode_distribution[combined_barcode_distribution["Dominance"]==1]
dd01 = combined_barcode_distribution[combined_barcode_distribution["Dominance"]==0.1]

print(combined_barcode_distribution[combined_barcode_distribution["barcode"]=="0"]['n'])

dilutions = [dd10, dd1, dd01]
import random
from plotly.subplots import make_subplots
for dilution in dilutions:
    fig = px.box(dilution,  x='n', y='percentage', color="n", color_discrete_sequence=px.colors.sequential.Greys_r[2:7],
                labels={'percentage': 'Percentage (%)'}, title=f"Prevalence of Clone 0: {dilution['Dominance'].iloc[0]}%")
    # Calculate significance test results for this dilution
    sample_test = significance_test(dilution, value_column='mean')
    # Add orange marker for dominant clone with annotation
    for n in dilution["n"].unique():
        dominant = dilution[(dilution["barcode"] == "0") & (dilution["n"] == n)]
        if not dominant.empty:
            annotation = get_annotation_text(n, sample_test)
            # Add the orange marker
            fig.add_trace(go.Scatter(
                y=dominant['percentage'],
                x=dominant['n'],
                mode='markers',
                marker=dict(color="orange", size=15, line=dict(color='black', width=2)),
                showlegend=False
            ))
            # Add annotation next to the plot with an arrow

            n_index = list(sorted(dilution["n"].unique())).index(n)
            y_offset = np.mean(dominant['percentage'].values) - n_index * np.mean(dominant['percentage'].values)

            fig.add_annotation(
                x=dominant['n'],
                y=dominant['percentage'].values[0] + y_offset,
                text=annotation,
                ax=50,  # Adjust the x offset for the annotation
                ay=0,
                font=dict(size=14, color="black"),
                bgcolor="white",
                bordercolor="black",
                borderwidth=1,
                showarrow=True
            )
    fig.update_xaxes(type='category', showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black') 
    fig.update_layout(width=500, height=500, font=dict(size=20), showlegend=False)
    fig.write_image(f"../ProgressReport/simulation_description/plots/boxplot_reads_dilution_{dilution['Dominance'].iloc[0]}.svg")
# %%
# insertion level OTB plot
ils = []
for path_group in paths:
    for n, dominance, i in path_group:
        data = read_data(i)
        summary = data.groupby(['mean_read_length', 'coverage', 'iteration', 'target']).agg(
        full_matches_total=('full_matches', 'sum'),
        partial_matches_total=('partial_matches', 'sum'),
        bases_on_target_total=('on_target_bases', 'sum')
        ).reset_index()
        final_summary = summary.groupby(['target','mean_read_length', 'coverage']).agg(
        full_matches_mean=('full_matches_total', 'mean'),
        full_matches_se=('full_matches_total', lambda x: x.std() / np.sqrt(len(x))),
        partial_matches_mean=('partial_matches_total', 'mean'),
        partial_matches_se=('partial_matches_total', lambda x: x.std() / np.sqrt(len(x))),
        on_target_bases_mean=('bases_on_target_total', 'median'),
        on_target_bases_se=('bases_on_target_total', lambda x: x.std() / np.sqrt(len(x)))
        ).reset_index()
        summary["Dominance"] = dominance
        summary["n"] = n
        summary["log_n"] = np.log10(n)
        ils.append(summary)

combined_ils = pd.concat(ils, ignore_index=True)

print(combined_ils.head())

#%%

ils_otb = combined_ils.groupby(['target', 'n', "Dominance"]).agg(bases_on_target_total=('bases_on_target_total', 'mean')).reset_index()
print(ils_otb.head())

#%%
# Get unique (n, Dominance) pairs, sorted by n increasing, dominance decreasing
unique_pairs = ils_otb[['n', "Dominance"]].drop_duplicates().sort_values(['n', 'Dominance'], ascending=[True, False]).values

# Create subplot grid: one column per (n, Dominance) pair
fig = make_subplots(
    rows=1,
    cols=len(unique_pairs),
    subplot_titles=[f"n={n}\nDominance={dominance}" for n, dominance in unique_pairs],
    shared_yaxes=False,  # Set to False for independent y-axes
    horizontal_spacing=0.05,
    x_title='Mean OTBs',
    y_title='Top 20 Target regions by OTBs'
)

for idx, (n, dominance) in enumerate(unique_pairs, start=1):
    target_data = ils_otb[(ils_otb['n'] == n) & (ils_otb["Dominance"] == dominance)]
    top_20_targets = target_data.nlargest(20, 'bases_on_target_total')
    # Assign colors: orange for targets containing "Barcode_0", grey for others
    top_20_targets['color'] = top_20_targets['target'].apply(
        lambda x: 'orange' if 'Barcode_0' in str(x) else 'grey'
    )
    fig.add_trace(
        go.Bar(
            y=top_20_targets['target'],
            x=top_20_targets['bases_on_target_total'],
            marker_color=top_20_targets['color'],
            orientation='h',
            width=0.5,
            showlegend=False
        ),
        row=1,
        col=idx
    )
    fig.update_xaxes(
        showline=True,
        linewidth=2,
        linecolor='black',
        tickangle=90,
        #title_text='OTBs',
        row=1,
        col=idx
    )
    fig.update_yaxes(
        showline=True,
        linewidth=2,
        linecolor='black',
        tickfont=dict(size=8, color='rgba(0,0,0,0)'),  # Hide tick labels by making them transparent
        automargin=True,
        tickmode='array',
        tickvals=[],  # Remove tick marks
        showticklabels=False,  # Hide tick labels
        #title_text="Top 20 Target regions by OTBs",
        row=1,
        col=idx
    )

fig.update_layout(
    width=200 * len(unique_pairs),
    height=400,
    font=dict(size=14),
    barmode='relative',
    showlegend=True
)
#fig.show()
fig.write_image(f"../ProgressReport/simulation_description/plots/ils_barplots.svg")
# %%
