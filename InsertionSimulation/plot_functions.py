#%%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import math
import numpy as np

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

pio.kaleido.scope.mathjax = None
pio.templates.default = "seaborn"

#sns.set_theme(style="ticks")

def get_barcode_color_mapping(barcodes):
    if isinstance(barcodes, (list, pd.Series, np.ndarray)):
        barcodes = pd.Series(barcodes)
    colors = px.colors.qualitative.Plotly  # Use Plotly's qualitative color palette
    barcode_colors = (colors * (len(barcodes) // len(colors) + 1))[:len(barcodes)]  # Repeat colors if needed
    return dict(zip(barcodes, barcode_colors))  # Return barcode-to-color mapping

def read_data(filepath):
    data = pd.read_csv(filepath, sep='\t')
    return data

#data = read_data('../out/tcr_20/tcr_20_matches_table.csv')
#basic_data = read_data('../out/tcr_20/tcr_20_barcode_distribution_table.csv')
#%%
def plot_barcode_distribution(data, output_path):

    #output
    output_path_total = os.path.join(output_path, "barplot_total_reads.svg")
    output_path_perc = os.path.join(output_path, "barplot_percentage_reads.svg")
    
    output_html_total = os.path.join(output_path, "barplot_total_reads.html")
    output_html_perc = os.path.join(output_path, "barplot_percentage_reads.html")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure the folder exists
    
    # Identify barcode columns based on their position before the 'coverage' column
    coverage_index = data.columns.get_loc('coverage')
    barcode_columns = data.columns[:coverage_index]

    # Group by coverage and mean_read_length and calculate mean across iterations
    grouped = data.groupby(['coverage', 'mean_read_length'])[barcode_columns].mean().reset_index()
    grouped['coverage_mean_read_length'] = grouped['coverage'].astype(str) + '_' + grouped['mean_read_length'].astype(str)
    grouped["mean"] = grouped[barcode_columns].mean(axis=1)
    grouped["sum"] = grouped[barcode_columns].sum(axis=1)
    
    # Interactive plot with Plotly
    fig = px.bar(grouped, x='coverage_mean_read_length', y='mean', title='Mean Total Reads by Coverage and Mean Read Length', color_discrete_sequence=['black'])
    fig.update_xaxes(title_text='Coverage and Mean Read Length', title_font=dict(size=12))
    fig.update_yaxes(title_text='Mean Total Reads', title_font=dict(size=12))
    fig.write_html(output_html_total)
    fig.write_image(output_path_total, scale=10)

    #percentage plot
    barcode_color_map = get_barcode_color_mapping(barcode_columns)

    for barcode in barcode_columns:
        grouped[barcode] = grouped[barcode] / grouped["sum"] * 100

    grouped_melted = grouped.melt(id_vars=['coverage_mean_read_length'], value_vars=barcode_columns, var_name='Barcode', value_name='Percentage')
    grouped_melted = grouped_melted.sort_values('Barcode', key=pd.to_numeric)

    # Interactive plot
    fig = px.bar(grouped_melted, x='coverage_mean_read_length', y='Percentage', color='Barcode', title='Percentage Share of Barcodes by Coverage and Mean Read Length', color_discrete_map=barcode_color_map)
    fig.update_xaxes(title_text='Coverage and Mean Read Length', title_font=dict(size=12))
    fig.update_yaxes(title_text='Percentage Share of Barcodes', title_font=dict(size=12))
    fig.write_html(output_html_perc)
    
    #Static plot
    fig.update_layout(legend=dict(font=dict(size=8), orientation="h", yanchor="bottom", y=-0.75, x=0.5, xanchor="center"))
    fig.write_image(output_path_perc, scale=10)

    return str(output_html_total), str(output_html_perc)

#plot_barcode_distribution(basic_data, output_path="./output/tcr_20/")  

# %%
def plot_lineplot(data, output_path):
    
    #output
    output_path_partial = os.path.join(output_path, "lineplot_partial_matches.svg")
    output_path_full = os.path.join(output_path, "lineplot_full_matches.svg")
    
    output_html_partial = os.path.join(output_path, "lineplot_partial_matches.html")
    output_html_full = os.path.join(output_path, "lineplot_full_matches.html")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure the folder exists

    data[["ID", "Barcode", "Iteration"]] = data["Insertion"].str.split("_", expand=True)

    numeric_cols = ["full_matches",
                    "partial_matches", 
                    "mean_read_length", 
                    "coverage",
                    "Barcode", 
                    "Iteration"]

    data = data[["ID"] + numeric_cols]

    data[numeric_cols] = data[numeric_cols].apply(pd.to_numeric, errors="coerce")
    
    # Partial matches plot
    grouped = data.groupby(['coverage', 'mean_read_length', 'ID'])['partial_matches'].mean().reset_index()
    
    # Static plot with seaborn
    plt.figure(figsize=(40, 10))
    sns.lineplot(data=grouped, x='ID', y='partial_matches', 
                 hue='coverage', style='mean_read_length', 
                 markers=True, dashes=False)
    plt.xlabel('ID')
    plt.ylabel('Mean Partial Matches')
    plt.title('Lineplot of Partial Matches by Coverage and Mean Read Length')
    plt.legend(title='Coverage and Mean Read Length')
    plt.xticks(rotation=90)
    plt.savefig(output_path_partial, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Interactive plot with Plotly for partial matches
    fig = px.line(grouped, x='ID', y='partial_matches', color='coverage', line_dash='mean_read_length', markers=True, title='Lineplot of Partial Matches by Coverage and Mean Read Length')
    fig.update_xaxes(title_text='ID', title_font=dict(size=12))
    fig.update_yaxes(title_text='Mean Partial Matches', title_font=dict(size=12))
    fig.write_html(output_html_partial)
    
    
    # Full matches plot
    grouped = data.groupby(['coverage', 'mean_read_length', 'ID'])['full_matches'].mean().reset_index()
    
    # Static plot with seaborn
    plt.figure(figsize=(40, 10))
    sns.lineplot(data=grouped, x='ID', y='full_matches', 
                 hue='coverage', style='mean_read_length', 
                 markers=True, dashes=False)  
    plt.xlabel('ID')
    plt.ylabel('Mean Full Matches')
    plt.title('Lineplot of Full Matches by Coverage and Mean Read Length')
    plt.legend(title='Coverage and Mean Read Length')
    plt.xticks(rotation=90)
    plt.savefig(output_path_full, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Interactive plot with Plotly
    fig = px.line(grouped, x='ID', y='full_matches', color='coverage', line_dash='mean_read_length', markers=True, title='Lineplot of Full Matches by Coverage and Mean Read Length')
    fig.update_xaxes(title_text='ID', title_font=dict(size=12))
    fig.update_yaxes(title_text='Mean Full Matches', title_font=dict(size=12))
    fig.write_html(output_html_full)

    return str(output_html_full), str(output_html_partial)
    
#plot_lineplot(data, output_path="./output/tcr_20/")

#%%
def plot_isolated_lineplot(data, output_path, filter=20, id_list=None):

    if filter > 20:
        print("Individual plots are limited to 20 IDs. Please reduce the filter value or select up to 20 specific IDs.")
        sys.exit()
    
    # Output paths
    output_path_partial = os.path.join(output_path, "panel_lineplot_partial_matches.svg")
    output_path_full = os.path.join(output_path, "panel_lineplot_full_matches.svg")
    output_html_partial = os.path.join(output_path, "panel_lineplot_partial_matches.html")
    output_html_full = os.path.join(output_path, "panel_lineplot_full_matches.html")

    # Use loaded data
    data[["ID", "Barcode", "Iteration"]] = data["Insertion"].str.split("_", expand=True)

    numeric_cols = ["full_matches",
                    "partial_matches", 
                    "mean_read_length", 
                    "coverage",
                    "Barcode", 
                    "Iteration"]

    data = data[["ID"] + numeric_cols]

    data[numeric_cols] = data[numeric_cols].apply(pd.to_numeric, errors="coerce")
    
    grouped_partial = data.groupby(['coverage', 'mean_read_length', 'Barcode', 'Iteration', 'ID'])['partial_matches'].mean().reset_index()
    grouped_full = data.groupby(['coverage', 'mean_read_length', 'Barcode', 'Iteration', 'ID'])['full_matches'].mean().reset_index()
    
    # Create a new column for the combination of coverage and mean_read_length
    grouped_partial['coverage_mean_read_length'] = grouped_partial['coverage'].astype(str) + '_' + grouped_partial['mean_read_length'].astype(str)
    grouped_full['coverage_mean_read_length'] = grouped_full['coverage'].astype(str) + '_' + grouped_full['mean_read_length'].astype(str)
    
    # Filter IDs
    if id_list:
        ids = grouped_partial['ID'].unique()
        filtered_ids = [id for id in ids if id in id_list]
        filtered_partial = grouped_partial[grouped_partial['ID'].isin(filtered_ids)]
        filtered_full = grouped_full[grouped_full['ID'].isin(filtered_ids)]
    else:
        first_ids = grouped_partial['ID'].unique()[:filter]
        filtered_partial = grouped_partial[grouped_partial['ID'].isin(first_ids)]
        filtered_full = grouped_full[grouped_full['ID'].isin(first_ids)]

    # Calculate mean and standard deviation for partial matches
    partial_stats = filtered_partial.groupby(['coverage_mean_read_length', 'ID', 'Barcode']).agg(
        mean_partial_matches=('partial_matches', 'mean'),
        std_partial_matches=('partial_matches', 'std')
    ).reset_index()

    # Calculate mean and standard deviation for full matches
    full_stats = filtered_full.groupby(['coverage_mean_read_length', 'ID', 'Barcode']).agg(
        mean_full_matches=('full_matches', 'mean'),
        std_full_matches=('full_matches', 'std')
    ).reset_index()

    #dimensions
    num_plots = len(filtered_partial['ID'].unique())
    
    if num_plots < 5:
        rows = 1
        cols = num_plots
    else:
        rows = math.ceil(math.sqrt(num_plots))  # At least a square root in rows
        cols = math.ceil(num_plots / rows)  # Distribute plots evenly

    """# Static plot for partial matches
    # Color
    barcode_color_map = get_barcode_color_mapping(grouped_full["Barcode"].unique())
    g = sns.FacetGrid(partial_stats, col='ID', col_wrap=cols, height=4, aspect=1.5)
    g.map(sns.lineplot, 'coverage_mean_read_length', 'mean_partial_matches', 'Barcode', markers=True, dashes=False, palette=barcode_color_map)
    g.add_legend()
    g.set_axis_labels('Coverage and Mean Read Length', 'Mean Partial Matches')
    for ax in g.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(90)
    g.fig.suptitle('Lineplot of Mean Partial Matches by Coverage and Mean Read Length', y=1.02)
    #os.makedirs(os.path.dirname(output_path_partial), exist_ok=True)  # Ensure the folder exists
    #plt.savefig(output_path_partial, dpi=300, bbox_inches="tight")
    #plt.close()
    
    # Static plot for full matches
    g = sns.FacetGrid(full_stats, col='ID', col_wrap=cols, height=4, aspect=1.5)
    g.map(sns.lineplot, 'coverage_mean_read_length', 'mean_full_matches', 'Barcode', markers=True, dashes=False, palette=barcode_color_map)
    g.add_legend()
    g.set_axis_labels('Coverage and Mean Read Length', 'Mean Full Matches')
    for ax in g.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(90)
    g.fig.suptitle('Lineplot of Mean Full Matches by Coverage and Mean Read Length', y=1.05)
    os.makedirs(os.path.dirname(output_path_full), exist_ok=True)  # Ensure the folder exists
    plt.savefig(output_path_full, dpi=300, bbox_inches="tight")
    plt.close()"""

    # Interactive plot for partial matches
    
    #color
    barcode_color_map = get_barcode_color_mapping(filtered_partial["Barcode"].unique())
    fig = make_subplots(rows=rows, cols=cols, shared_yaxes=True, subplot_titles=filtered_partial['ID'].unique())
    for i, unique_id in enumerate(filtered_partial['ID'].unique(), start=1):
        subset = partial_stats[partial_stats['ID'] == unique_id]
        row = (i - 1) // cols + 1
        col = (i - 1) % cols + 1
        for barcode in subset['Barcode'].unique():
            barcode_data = subset[subset['Barcode'] == barcode]
            fig.add_trace(go.Scatter(x=barcode_data['coverage_mean_read_length'], y=barcode_data['mean_partial_matches'], mode='lines+markers', name=str(barcode), legendgroup=str(barcode), showlegend=(i == 1), line=dict(color=barcode_color_map[barcode])), row=row, col=col)
    fig.update_xaxes(title_text='Coverage and Mean Read Length', title_font=dict(size=8), title_standoff=5)
    fig.update_yaxes(title_text='Mean Partial Matches', title_font=dict(size=8), title_standoff=5)
    fig.update_layout(title_text='Lineplot of Mean Partial Matches by Coverage and Mean Read Length', showlegend=True)
    fig.write_html(output_html_partial)
    fig.write_image(output_path_partial, scale=3, width=1200, height=1200)

    # Interactive plot for full matches
    fig = make_subplots(rows=rows, cols=cols, shared_yaxes=True, subplot_titles=filtered_full['ID'].unique())
    for i, unique_id in enumerate(filtered_full['ID'].unique(), start=1):
        subset = full_stats[full_stats['ID'] == unique_id]
        row = (i - 1) // cols + 1
        col = (i - 1) % cols + 1
        for barcode in subset['Barcode'].unique():
            barcode_data = subset[subset['Barcode'] == barcode]
            fig.add_trace(go.Scatter(x=barcode_data['coverage_mean_read_length'], y=barcode_data['mean_full_matches'], mode='lines+markers', name=str(barcode), legendgroup=str(barcode), showlegend=(i == 1), line=dict(color=barcode_color_map[barcode])), row=row, col=col)
    fig.update_xaxes(title_text='Coverage and Mean Read Length', title_font=dict(size=8), title_standoff=5)
    fig.update_yaxes(title_text='Mean Partial Matches', title_font=dict(size=8), title_standoff=5)
    fig.update_layout(title_text='Lineplot of Mean Full Matches by Coverage and Mean Read Length', showlegend=True)
    fig.write_html(output_html_full)
    fig.write_image(output_path_full, scale=3, width=1200, height=1200)

    return str(output_html_full), str(output_html_partial)

#plot_isolated_lineplot(data, output_path="./output/tcr_20/", filter=2, id_list=["TRAC", "TRBC1", "TRBV1", "TRAV3"])
#plot_isolated_lineplot(data, output_path="./output/tcr_20/", filter=20)

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
                margin: 30px;
                text-align: center; 
                background-color: #f4f4f4;
            }}
            h1, h2, h3 {{ color: #333; }}
            h3 {{ margin-top: 40px; border-bottom: 2px solid #444; padding-bottom: 5px; }}
            .grid-container {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(800px, 1fr));
                gap: 50px;
                justify-content: center;
                margin: 20px auto;
                max-width: 90%;
            }}
            .grid-container iframe {{
                width: 100%;
                height: 800px;
                border: 5px solid #ccc;
                border-radius: 8px;
                background: white;
                box-shadow: 2px 2px 10px rgba(0,0,0,0.1);
                transition: transform 0.2s ease-in-out;
            }}
            .grid-container iframe:hover {{
                transform: scale(1.02);
                border-color: #007bff;
            }}
            .config-box {{
                text-align: left;
                background: #fff;
                padding: 15px;
                border-radius: 5px;
                box-shadow: 2px 2px 10px rgba(0,0,0,0.1);
                max-width: 80%;
                margin: 20px auto;
                overflow: auto;
                white-space: pre-wrap;
                border-left: 4px solid #ccc;
            }}
        </style>
    </head>
    <body>
        <h1>Analysis Report</h1>
        
        <h3>Overview Plots</h3>
        <div class="grid-container">
            <iframe src="{image_paths[0]}" title="Total Reads"></iframe>
            <iframe src="{image_paths[1]}" title="Percentage Reads"></iframe>
        </div>

        <h3>Match Statistics</h3>
        <div class="grid-container">
            <iframe src="{image_paths[2]}" title="Partial Matches"></iframe>
            <iframe src="{image_paths[3]}" title="Full Matches"></iframe>
        </div>

        <h3>Region-Specific Read Overlaps</h3>
        <div class="grid-container">
            <iframe src="{image_paths[4]}" title="Full Match Panel"></iframe>
            <iframe src="{image_paths[5]}" title="Partial Match Panel"></iframe>
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

# Run the function
'''
generate_html_report("./output/tcr_20/",
    config="../sim_config_roi.ini", 
    output_html="./output/report.html")


   image_paths=[
        "./tcr_20/barplot_total_reads.html",
        "./tcr_20/barplot_percentage_reads.html",
        "./tcr_20/lineplot_full_matches.html", 
        "./tcr_20/lineplot_partial_matches.html",
        "./tcr_20/panel_lineplot_full_matches.html",
        "./tcr_20/panel_lineplot_partial_matches.html"
    ], 

'''