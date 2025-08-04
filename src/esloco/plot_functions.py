#%%
import os
import sys
import math
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from datetime import datetime
pio.templates.default = "plotly_white"
from esloco.coverage import get_barcode_color_mapping

def read_data(filepath):
    data = pd.read_csv(filepath, sep='\t')
    return data

#%%
def barplot_absolute_matches(experiment_name, data, output_path):
    # Output paths
    output_svg = os.path.join(output_path, f"{experiment_name}_Barplot_absolute_numbers.svg")
    output_html = os.path.join(output_path, f"{experiment_name}_Barplot_absolute_numbers.html")
    # Extract iteration number from the 'Insertion' column
    data['Iteration'] = data['target_region'].str.extract(r'_(\d+)$').astype(int)

    # Group by mean_read_length, coverage, and iteration # sum up all barcodes
    summary = data.groupby(['mean_read_length', 'coverage', 'Iteration']).agg(
        full_matches_total=('full_matches', 'sum'),
        partial_matches_total=('partial_matches', 'sum'),
        bases_on_target_total=('on_target_bases', 'sum')
    ).reset_index()
    # Aggregate across iterations to compute mean and standard error
    final_summary = summary.groupby(['mean_read_length', 'coverage']).agg(
        full_matches_mean=('full_matches_total', 'mean'),
        full_matches_se=('full_matches_total', lambda x: x.std() / np.sqrt(len(x))),
        partial_matches_mean=('partial_matches_total', 'mean'),
        partial_matches_se=('partial_matches_total', lambda x: x.std() / np.sqrt(len(x))),
        on_target_bases_mean=('bases_on_target_total', 'mean'),
        on_target_bases_se=('bases_on_target_total', lambda x: x.std() / np.sqrt(len(x)))
    ).reset_index()

    #print(final_summary)
    x_labels = [f"{row['mean_read_length']}, {row['coverage']}" for _, row in final_summary.iterrows()]
    x = np.arange(len(x_labels))

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=x,
        y=final_summary['full_matches_mean'],
        error_y=dict(type='data', array=final_summary['full_matches_se']),
        name='Full Matches',
        marker_color='black'
    ))

    fig.add_trace(go.Bar(
        x=x,
        y=final_summary['partial_matches_mean'],
        error_y=dict(type='data', array=final_summary['partial_matches_se']),
        name='Partial Matches',
        marker_color='grey'
    ))

    fig.add_trace(go.Bar(
        x=x,
        y=final_summary['on_target_bases_mean'],
        error_y=dict(type='data', array=final_summary['on_target_bases_se']),
        name='OTBs',
        marker_color='lightgrey'
    ))
    fig.update_layout(
        title='Full and Partial Matches Across Conditions',
        xaxis=dict(
            tickmode='array',
            tickvals=x,
            ticktext=x_labels,
            title='Mean Read Length, Coverage'
        ),
        yaxis=dict(
            title='Mean Count across Iterations'
        ),
        barmode='group'
    )

    #fig.show()
    fig.write_html(output_html)
    fig.write_image(output_svg, width=600, height=400)

    print(f"Barplot absolute numbers saved as {output_html}")
    return str(output_html)
#%%
def barplot_absolute_matches_barcodes(experiment_name, data, output_path):
     # Output paths
    output_svg = os.path.join(output_path, f"{experiment_name}_Barplot_Barcode_absolute_numbers.svg")
    output_html = os.path.join(output_path, f"{experiment_name}_Barplot_Barcode_absolute_numbers.html")

    data['Barcode'] = data['target_region'].str.extract(r'_(\d+)_').astype(int)
    data['Iteration'] = data['target_region'].str.extract(r'_(\d+)$').astype(int)

    barcode_color_map = get_barcode_color_mapping(data['Barcode'].unique())
    # Stacked bar plot colored by barcodes
    # Group by mean_read_length, coverage, and barcode #sum across iterations
    barcode_summary = data.groupby(['mean_read_length', 'coverage', "Barcode"]).agg(
        full_matches_total=('full_matches', 'sum'),
        partial_matches_total=('partial_matches', 'sum'),
        on_target_bases_total=('on_target_bases', 'sum')
    ).reset_index()

    # Divide the values by the number of iterations
    iterations_count = data['Iteration'].nunique()

    barcode_summary['full_matches_total'] /= iterations_count
    barcode_summary['partial_matches_total'] /= iterations_count
    barcode_summary['on_target_bases_total'] /= iterations_count

    # Melt the dataframe for easier plotting with Plotly
    barcode_summary_melted = barcode_summary.melt(id_vars=['mean_read_length', 'coverage', 'Barcode'],
                                                value_vars=['full_matches_total', 'partial_matches_total', 'on_target_bases_total'],
                                                var_name='match_type', value_name='count')
    # Create combined bar plot
    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=['Full Matches', 'Partial Matches', 'Bases on Target'],
        shared_yaxes=False,
        shared_xaxes=True,
        #x_title="Coverage and Mean Read Length",
        y_title="Mean Count across Iterations"
    )
    fig.update_annotations(font_size=12)
    for match_type in barcode_summary_melted['match_type'].unique():
        subset = barcode_summary_melted[barcode_summary_melted['match_type'] == match_type]
        subset['color'] = subset['Barcode'].map(barcode_color_map)
        col = 1 if match_type == 'full_matches_total' else (2 if match_type == 'partial_matches_total' else 3)
        for barcode in subset['Barcode'].unique():
            barcode_data = subset[subset['Barcode'] == barcode]
            fig.add_trace(go.Bar(
                x=barcode_data['mean_read_length'].astype(str) + ', ' + barcode_data['coverage'].astype(str),
                y=barcode_data['count'],
                name=f"{barcode}",
                marker_color=barcode_color_map[barcode],
                showlegend=(col == 1)
            ), row=1, col=col)

    fig.update_xaxes(title_text='Coverage, Mean Read Length', title_font=dict(size=8), title_standoff=5)
    fig.update_layout(title_text='Mean Count (by barcode)', showlegend=True)
    #fig.show()
    fig.write_html(output_html)
    fig.write_image(output_svg, width=1200, height=400)

    print(f"Barplot absolute numbers with barcodes saved as {output_html}")
    return str(output_html)


#%%
def plot_barcode_distribution(experiment_name, data, output_path):
    # Output paths
    output_path_total = os.path.join(output_path, f"{experiment_name}_barplot_total_reads.svg")
    output_html_total = os.path.join(output_path, f"{experiment_name}_barplot_total_reads.html")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure the folder exists

    # Identify barcode columns based on their position before the 'coverage' column
    coverage_index = data.columns.get_loc('coverage')
    barcode_columns = data.columns[:coverage_index]

    # Group by coverage and mean_read_length and calculate mean across iterations
    grouped = data.groupby(['coverage', 'mean_read_length'])[barcode_columns].mean().reset_index()
    grouped['coverage_mean_read_length'] = grouped['coverage'].astype(str) + '_' + grouped['mean_read_length'].astype(str)
    grouped["sum"] = grouped[barcode_columns].sum(axis=1)

    # Melt the dataframe for easier plotting
    grouped_melted = grouped.melt(id_vars=['coverage_mean_read_length', 'sum'],
                                  value_vars=barcode_columns,
                                  var_name='barcode',
                                  value_name='count')
    grouped_melted = grouped_melted.sort_values('barcode', key=pd.to_numeric)

    # Color mapping for barcodes
    barcode_color_map = get_barcode_color_mapping(barcode_columns)

    # Interactive stacked bar plot with Plotly
    fig = px.bar(grouped_melted,
                 x='coverage_mean_read_length',
                 y='count',
                 color='barcode',
                 title='Mean Total Reads by Coverage and Mean Read Length (Stacked by Barcode)',
                 color_discrete_map=barcode_color_map)
    fig.update_xaxes(title_text='Coverage, Mean Read Length', title_font=dict(size=12))
    fig.update_yaxes(title_text='Mean Total Reads', title_font=dict(size=12))
    fig.show()
    fig.write_html(output_html_total)
    # Static plot
    fig.update_layout(legend=dict(font=dict(size=8), orientation="h", yanchor="bottom", y=-0.75, x=0.5, xanchor="center"))
    fig.write_image(output_path_total, scale=10)

    print(f"Mean Total Reads Stacked Barplot saved as {output_html_total}")
    return str(output_html_total)

# %%
def plot_lineplot(experiment_name, data, output_path):
    #output
    output_path_partial = os.path.join(output_path, f"{experiment_name}_lineplot_partial_matches.svg")
    output_path_full = os.path.join(output_path, f"{experiment_name}_lineplot_full_matches.svg")
    output_path_otb = os.path.join(output_path, f"{experiment_name}_lineplot_otb_matches.svg")
    output_html_partial = os.path.join(output_path, f"{experiment_name}_lineplot_partial_matches.html")
    output_html_full = os.path.join(output_path, f"{experiment_name}_lineplot_full_matches.html")
    output_html_otb = os.path.join(output_path, f"{experiment_name}_lineplot_otb_matches.html")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure the folder exists

    #examples
    #I mode: Barcode_0_insertion_0_0
    #ROI mode: TRA1_0_0
    if data["target_region"].str.contains("insertion").any():
        data[["temp1","barcode","ID1","ID2","Iteration"]] = data["target_region"].str.split("_", expand=True)
        data["id"] = data["ID1"] + "_" + data["ID2"]
    else:
        split_data = data["target_region"].str.rsplit("_", n=2)
        data["id"] = split_data.str[0]
        data["barcode"] = split_data.str[1]
        data["Iteration"] = split_data.str[2]

    numeric_cols = ["full_matches",
                    "partial_matches",
                    "on_target_bases", 
                    "mean_read_length", 
                    "coverage",
                    "barcode", 
                    "Iteration"]

    data = data[["id"] + numeric_cols]

    data.loc[:, numeric_cols] = data[numeric_cols].apply(pd.to_numeric, errors="coerce")

    # Partial matches plot
    grouped = data.groupby(['coverage', 'mean_read_length', 'id'])['partial_matches'].mean().reset_index()

    # Static plot with seaborn
    plt.figure(figsize=(40, 10))
    sns.lineplot(data=grouped, x='id', y='partial_matches',
                 hue='coverage', style='mean_read_length',
                 markers=True, dashes=False)
    plt.xlabel('ID')
    plt.ylabel('Mean Partial Matches')
    plt.title('Partial Matches by Coverage and Mean Read Length')
    plt.legend(title='Coverage and Mean Read Length')
    plt.xticks(rotation=90)
    plt.savefig(output_path_partial, dpi=300, bbox_inches="tight")
    plt.close()

    # Interactive plot with Plotly for partial matches
    fig = px.line(grouped, x='id', y='partial_matches', color='coverage', line_dash='mean_read_length', markers=True, title='Target-specific Partial Matches (by Coverage, Mean Read Length)')
    fig.update_layout(title_font=dict(size=14))
    fig.update_xaxes(title_text='ID', title_font=dict(size=12))
    fig.update_yaxes(title_text='Mean Partial Matches', title_font=dict(size=12))
    fig.write_html(output_html_partial)

    # Full matches plot
    grouped = data.groupby(['coverage', 'mean_read_length', 'id'])['full_matches'].mean().reset_index()

    # Static plot with seaborn
    plt.figure(figsize=(40, 10))
    sns.lineplot(data=grouped, x='id', y='full_matches',
                 hue='coverage', style='mean_read_length',
                 markers=True, dashes=False)
    plt.xlabel('ID')
    plt.ylabel('Mean Full Matches')
    plt.title('Full Matches by Coverage and Mean Read Length')
    plt.legend(title='Coverage and Mean Read Length')
    plt.xticks(rotation=90)
    plt.savefig(output_path_full, dpi=300, bbox_inches="tight")
    plt.close()

    # Interactive plot with Plotly
    fig = px.line(grouped, x='id', y='full_matches', color='coverage', line_dash='mean_read_length', markers=True, title='Target-specific Full Matches (by Coverage, Mean Read Length)')
    fig.update_layout(title_font=dict(size=14))
    fig.update_xaxes(title_text='ID', title_font=dict(size=12))
    fig.update_yaxes(title_text='Mean Full Matches', title_font=dict(size=12))
    fig.write_html(output_html_full)

    # OTB plot
    grouped = data.groupby(['coverage', 'mean_read_length', 'id'])['on_target_bases'].mean().reset_index()

    # Static plot with seaborn
    plt.figure(figsize=(40, 10))
    sns.lineplot(data=grouped, x='id', y='on_target_bases',
                 hue='coverage', style='mean_read_length',
                 markers=True, dashes=False)
    plt.xlabel('ID')
    plt.ylabel('Mean OTBs')
    plt.title('OTBs by Coverage and Mean Read Length')
    plt.legend(title='Coverage and Mean Read Length')
    plt.xticks(rotation=90)
    plt.savefig(output_path_otb, dpi=300, bbox_inches="tight")
    plt.close()

    # Interactive plot with Plotly
    fig = px.line(grouped, x='id', y='on_target_bases', color='coverage', line_dash='mean_read_length', markers=True, title='Target-specific OTBs (by Coverage, Mean Read Length)')
    fig.update_layout(title_font=dict(size=14))
    fig.update_xaxes(title_text='ID', title_font=dict(size=12))
    fig.update_yaxes(title_text='Mean OTBs', title_font=dict(size=12))
    fig.write_html(output_html_otb)

    print(f"Partial Lineplot saved as {output_html_partial}")
    print(f"Full Lineplot saved as {output_html_full}")
    print(f"OTBs Lineplot saved as {output_html_otb}")
    return str(output_html_full), str(output_html_partial), str(output_html_otb)

#%%
def plot_isolated_lineplot(experiment_name, data, output_path, filterplots=20, id_list=[]):

    if filterplots > 20 or len(id_list) > 20:
        print("Individual plots are limited to 20 IDs. Please reduce the filter value or select up to 20 specific IDs in `plot_esloco` itself.")
        sys.exit()

    # Output paths
    output_path_partial = os.path.join(output_path, f"{experiment_name}_panel_lineplot_partial_matches.svg")
    output_path_full = os.path.join(output_path, f"{experiment_name}_panel_lineplot_full_matches.svg")
    output_path_otb = os.path.join(output_path, f"{experiment_name}_panel_lineplot_otb_matches.svg")
    output_html_partial = os.path.join(output_path, f"{experiment_name}_panel_lineplot_partial_matches.html")
    output_html_full = os.path.join(output_path, f"{experiment_name}_panel_lineplot_full_matches.html")
    output_html_otb = os.path.join(output_path, f"{experiment_name}_panel_lineplot_otb_matches.html")

    # Use loaded data
    if data["target_region"].str.contains("insertion").any():
        data[["temp1","barcode","ID1","ID2","Iteration"]] = data["target_region"].str.split("_", expand=True)
        data["id"] = data["ID1"] + "_" + data["ID2"]
    else:
        split_data = data["target_region"].str.rsplit("_", n=2)
        data["id"] = split_data.str[0]
        data["barcode"] = split_data.str[1]
        data["Iteration"] = split_data.str[2]

    numeric_cols = ["full_matches",
                    "partial_matches",
                    "on_target_bases",
                    "mean_read_length",
                    "coverage",
                    "barcode",
                    "Iteration"]

    data = data[["id"] + numeric_cols]

    data.loc[:, numeric_cols] = data[numeric_cols].apply(pd.to_numeric, errors="coerce")

    grouped_partial = data.groupby(['coverage', 'mean_read_length', 'barcode', 'Iteration', 'id'])['partial_matches'].mean().reset_index()
    grouped_full = data.groupby(['coverage', 'mean_read_length', 'barcode', 'Iteration', 'id'])['full_matches'].mean().reset_index()
    grouped_otb = data.groupby(['coverage', 'mean_read_length', 'barcode', 'Iteration', 'id'])['on_target_bases'].mean().reset_index()

    # Create a new column for the combination of coverage and mean_read_length
    grouped_partial['coverage_mean_read_length'] = grouped_partial['coverage'].astype(str) + '_' + grouped_partial['mean_read_length'].astype(str)
    grouped_full['coverage_mean_read_length'] = grouped_full['coverage'].astype(str) + '_' + grouped_full['mean_read_length'].astype(str)
    grouped_otb['coverage_mean_read_length'] = grouped_otb['coverage'].astype(str) + '_' + grouped_otb['mean_read_length'].astype(str)

    # Filter IDs
    if id_list:
        ids = grouped_partial['id'].unique()
        filtered_ids = [id for id in ids if id in id_list]
        filtered_partial = grouped_partial[grouped_partial['id'].isin(filtered_ids)]
        filtered_full = grouped_full[grouped_full['id'].isin(filtered_ids)]
        filtered_otb = grouped_otb[grouped_otb['id'].isin(filtered_ids)]
    else:
        first_ids = grouped_partial['id'].unique()[:filterplots]
        filtered_partial = grouped_partial[grouped_partial['id'].isin(first_ids)]
        filtered_full = grouped_full[grouped_full['id'].isin(first_ids)]
        filtered_otb = grouped_otb[grouped_otb['id'].isin(first_ids)]

    # Calculate mean and standard deviation for partial matches
    partial_stats = filtered_partial.groupby(['coverage_mean_read_length', 'id', 'barcode']).agg(
        mean_partial_matches=('partial_matches', 'mean'),
        std_partial_matches=('partial_matches', 'std')
    ).reset_index()

    # Calculate mean and standard deviation for full matches
    full_stats = filtered_full.groupby(['coverage_mean_read_length', 'id', 'barcode']).agg(
        mean_full_matches=('full_matches', 'mean'),
        std_full_matches=('full_matches', 'std')
    ).reset_index()

    # Calculate mean and standard deviation for full matches
    otb_stats = filtered_otb.groupby(['coverage_mean_read_length', 'id', 'barcode']).agg(
        mean_otb_matches=('on_target_bases', 'mean'),
        std_otb_matches=('on_target_bases', 'std')
    ).reset_index()

    #dimensions
    num_plots = len(filtered_partial['id'].unique())

    if num_plots < 5:
        rows = 1
        cols = num_plots
    else:
        rows = math.ceil(math.sqrt(num_plots))  # At least a square root in rows
        cols = math.ceil(num_plots / rows)  # Distribute plots evenly

    # Interactive plot for partial matches
    #color
    barcode_color_map = get_barcode_color_mapping(filtered_partial["barcode"].unique())
    fig = make_subplots(rows=rows, cols=cols, shared_yaxes=True, subplot_titles=filtered_partial['id'].unique())
    for i, unique_id in enumerate(filtered_partial['id'].unique(), start=1):
        subset = partial_stats[partial_stats['id'] == unique_id]
        row = (i - 1) // cols + 1
        col = (i - 1) % cols + 1
        for barcode in subset['barcode'].unique():
            barcode_data = subset[subset['barcode'] == barcode]
            fig.add_trace(go.Scatter(x=barcode_data['coverage_mean_read_length'],
                                     y=barcode_data['mean_partial_matches'],
                                     mode='lines+markers',
                                     name=str(barcode),
                                     legendgroup=str(barcode),
                                     showlegend=(i == 1),
                                     line=dict(color=barcode_color_map[barcode])),
                                     row=row,
                                     col=col)
    fig.update_xaxes(title_text='Coverage, Mean Read Length', title_font=dict(size=8), title_standoff=5)
    fig.update_yaxes(title_text='Mean Partial Matches', title_font=dict(size=8), title_standoff=5)
    fig.update_layout(title_text='Target-specific Partial Matches (by barcode)', showlegend=True)
    fig.write_html(output_html_partial)
    fig.write_image(output_path_partial, scale=3, width=1200, height=1200)

    # Interactive plot for full matches
    fig = make_subplots(rows=rows, cols=cols, shared_yaxes=True, subplot_titles=filtered_full['id'].unique())
    for i, unique_id in enumerate(filtered_full['id'].unique(), start=1):
        subset = full_stats[full_stats['id'] == unique_id]
        row = (i - 1) // cols + 1
        col = (i - 1) % cols + 1
        for barcode in subset['barcode'].unique():
            barcode_data = subset[subset['barcode'] == barcode]
            fig.add_trace(go.Scatter(x=barcode_data['coverage_mean_read_length'],
                                     y=barcode_data['mean_full_matches'],
                                     mode='lines+markers',
                                     name=str(barcode),
                                     legendgroup=str(barcode),
                                     showlegend=(i == 1),
                                     line=dict(color=barcode_color_map[barcode])),
                                     row=row,
                                     col=col)
    fig.update_xaxes(title_text='Coverage, Mean Read Length', title_font=dict(size=8), title_standoff=5)
    fig.update_yaxes(title_text='Mean Full Matches', title_font=dict(size=8), title_standoff=5)
    fig.update_layout(title_text='Target-specific Full Matches (by barcode)', showlegend=True)
    fig.write_html(output_html_full)
    fig.write_image(output_path_full, scale=3, width=1200, height=1200)

    # Interactive plot for otb matches
    fig = make_subplots(rows=rows, cols=cols, shared_yaxes=True, subplot_titles=filtered_otb['id'].unique())
    for i, unique_id in enumerate(filtered_otb['id'].unique(), start=1):
        subset = otb_stats[otb_stats['id'] == unique_id]
        row = (i - 1) // cols + 1
        col = (i - 1) % cols + 1
        for barcode in subset['barcode'].unique():
            barcode_data = subset[subset['barcode'] == barcode]
            fig.add_trace(go.Scatter(x=barcode_data['coverage_mean_read_length'],
                                     y=barcode_data['mean_otb_matches'],
                                     mode='lines+markers',
                                     name=str(barcode),
                                     legendgroup=str(barcode),
                                     showlegend=(i == 1),
                                     line=dict(color=barcode_color_map[barcode])),
                                     row=row,
                                     col=col)
    fig.update_xaxes(title_text='Coverage, Mean Read Length', title_font=dict(size=8), title_standoff=5)
    fig.update_yaxes(title_text='Mean OTB Matches', title_font=dict(size=8), title_standoff=5)
    fig.update_layout(title_text='Target-specific OTB (by barcode)', showlegend=True)
    fig.write_html(output_html_otb)
    fig.write_image(output_path_otb, scale=3, width=1200, height=1200)

    print(f"Partial Overalps Lineplot panel saved as {output_html_partial}")
    print(f"Full Overlaps Lineplot panel saved as {output_html_full}")
    print(f"OTB Overlaps Lineplot panel saved as {output_html_otb}")
    return str(output_html_full), str(output_html_partial), str(output_html_otb)

#%%
# log file plot
def parse_log(log_file):
    """Parses a log file to extract timestamps, memory usage, execution times, and iteration labels for specific iterations."""
    log_data = {}
    with open(log_file, 'r') as file:
        for line in file:
            process_match = re.search(r'Process: (\w+)', line)
            if not process_match:
                continue
            timestamp_match = re.search(r'^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d+)', line)
            cpu_match = re.search(r'CPU: (\d+\.\d+)%', line)
            memory_match = re.search(r'Memory: (\d+\.\d+)%', line)
            if timestamp_match:
                timestamp = datetime.strptime(timestamp_match.group(1).split(',')[0], "%Y-%m-%d %H:%M:%S")
                formatted_timestamp = timestamp.strftime("%Y-%m-%d-%H-%M-%S")
                if formatted_timestamp not in log_data:
                    log_data[formatted_timestamp] = {'memory_usage': None, 'cpu_usage': None, 'label': None}
            if memory_match:
                log_data[formatted_timestamp]['memory_usage'] = float(memory_match.group(1))
            if cpu_match:
                log_data[formatted_timestamp]['cpu_usage'] = float(cpu_match.group(1))
            log_data[formatted_timestamp]['label'] = process_match.group(1)
    return log_data

#%%
# Extract data for plotting
def plot_log_data(experiment_name, logfile, output_path):

    process_data=parse_log(logfile)

    timestamps = sorted(process_data.keys(), key=lambda x: datetime.strptime(x, "%Y-%m-%d-%H-%M-%S"))
    memory_usage = [process_data[timestamp]['memory_usage'] for timestamp in timestamps]
    cpu_usage = [process_data[timestamp]['cpu_usage'] for timestamp in timestamps]

    # Plot with Plotly and save as HTML
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scatter(x=timestamps,
                             y=memory_usage,
                             mode='lines+markers',
                             name='Memory Usage (%)',
                             line=dict(color='black', width=4),
                             marker=dict(size=8)),
                             secondary_y=False)
    fig.add_trace(go.Scatter(x=timestamps,
                             y=cpu_usage,
                             mode='lines+markers',
                             name='CPU Usage (%)',
                             line=dict(color='red', width=4),
                             marker=dict(size=8)),
                             secondary_y=True)

    # Add semi-transparent red box between 80 and 100
    fig.add_shape(type="rect", x0=timestamps[0], x1=timestamps[-1], y0=80, y1=100,
                  fillcolor="red", opacity=0.1, line_width=0, secondary_y=True)

    # Add semi-transparent yellow box between 60 and 80
    fig.add_shape(type="rect", x0=timestamps[0], x1=timestamps[-1], y0=60, y1=80,
                  fillcolor="yellow", opacity=0.1, line_width=0, secondary_y=False)

    fig.update_xaxes(title_text='Timestamp')
    fig.update_yaxes(title_text='Memory Usage (%)', range=[0, 100], secondary_y=False)
    fig.update_yaxes(title_text='CPU Usage (%)', range=[0, 100], secondary_y=True)
    fig.update_layout(title_text='Memory and CPU Usage Over Time', legend=dict(orientation="h", yanchor="bottom", y=1.05, xanchor="right", x=1))

    html_output_path = os.path.join(output_path, f"{experiment_name}_log_plot.html")
    fig.write_html(html_output_path)

    svg_output_path = os.path.join(output_path, f"{experiment_name}_log_plot.svg")
    fig.write_image(svg_output_path, format='svg')

    print(f"Ressource plot saved as {html_output_path}")
    return html_output_path
#%%

#%%
def generate_html_report(image_paths, config=None, output_html="report.html"):
    """
    Creates an improved HTML report with embedded Plotly HTML plots.

    Parameters:
    - image_paths: List of HTML plot file paths.
    - config: Optional path to a configuration file.
    - output_html: Output HTML file name.
    """

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Analysis Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 25px;
                text-align: center;
                background-color:  #d8dadb;
            }}
            h1, h2, h3 {{ color: #333; }}
            h3 {{ margin-top: 40px; border-bottom: 3px solid #444; padding-bottom: 5px; }}
            .grid-container {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 50px;
                justify-content: center;
                margin: 20px auto;
                max-width: 80%;
            }}
            .grid-container iframe {{
                width: 100%;
                height: 400px;
                border: none;
                border-radius: 25px;
                background: white;
                box-shadow: 4px 4px 10px rgba(0,0,0,0.1);
                transition: transform 0.2s ease-in-out;
            }}
            .config-box {{
                text-align: left;
                background: #fff;
                padding: 15px;
                border-radius: 25px;
                box-shadow: 4px 4px 10px rgba(0,0,0,0.1);
                max-width: 50%;
                margin: 20px auto;
                overflow: auto;
                white-space: pre-wrap;
            }}
        </style>
    </head>
    <body>
        <h1>Analysis Report</h1>
        <p>Results in "{output_html}"...</p>
        <h3>General Overview</h3>
        <div class="grid-container">
            <iframe src="{image_paths[0]}" title="Total Reads"></iframe>
            <iframe src="{image_paths[1]}" title="Percentage Reads"></iframe>
            
        </div>
        <h3>Match Overview</h3>
        <div class="grid-container">
            <iframe src="{image_paths[2]}" title="..."></iframe>
            <iframe src="{image_paths[3]}" title="Full Matches"></iframe>
            <iframe src="{image_paths[4]}" title="Partial Matches"></iframe>
        </div>
        <div class="grid-container">

        </div>
        <h3>Region-Specific Read Overlaps</h3>
        <div class="grid-container">
            <iframe src="{image_paths[5]}" title="Full Match Panel" style="height: 1200px; border: none;"></iframe>
            <iframe src="{image_paths[6]}" title="Partial Match Panel" style="height: 1200px; border: none;"></iframe>
            <iframe src="{image_paths[7]}" title="OTB Match Panel" style="height: 1200px; border: none;"></iframe>
        </div>
    """
    html_content += f"""
    <h3>Coverage Overview</h3>
    <div class="grid-container">
        <button onclick="prevPlot()">Previous</button>
        <button onclick="nextPlot()">Next</button>
    </div>
    <div class="grid-container">
        <iframe id="plotFrame" src="{image_paths[9]}" title="Coverage Plot" style="height: 800px;"></iframe>
    </div>
    <script>
        var plots = {image_paths};
        var currentPlotIndex = 9;

        function showPlot(index) {{
            document.getElementById('plotFrame').src = plots[index];
        }}

        function prevPlot() {{
            if (currentPlotIndex > 9) {{
                currentPlotIndex--;
                showPlot(currentPlotIndex);
            }}
        }}

        function nextPlot() {{
            if (currentPlotIndex < plots.length - 1) {{
                currentPlotIndex++;
                showPlot(currentPlotIndex);
            }}
        }}
    </script>
    """

    html_content += f"""
        <h3>CPU and Memory Usage</h3>
        <div class="grid-container">
            <iframe src="{image_paths[8]}" title="Ressources"></iframe>
        </div>
        </div>
    """

    # Include config file content
    if config:
        with open(config, 'r') as file:
            config_content = file.read()
        html_content += f"""
        <h3>Simulation Configuration Parameters</h3>
        <div class="config-box">
            <pre>{config_content}</pre>
        </div>
        """

    html_content += "</body></html>"

    # Save HTML file
    with open(output_html, "w") as f:
        f.write(html_content)

    print(f"HTML report saved as {output_html}")
