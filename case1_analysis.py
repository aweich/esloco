#%%
import pandas as pd
from scipy.stats import pearsonr
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

#%% preprocessing the simulation output
raw = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/Case1_bias/Case1_bias2_matches_only.csv", sep="\t", header=0)
print(raw["Insertion"].unique())

# Split using the last two underscores to handle problematic entries like 'GBA_x_0_999'
raw[["gene", "Barcode", "Iteration"]] = raw["Insertion"].str.rsplit("_", n=2, expand=True)

print(raw["gene"].unique())
#%%
print(raw.head())
raw = raw.drop(columns=["Insertion", "Barcode"])
print(raw.head())

raw["Iteration"] = raw["Iteration"].astype(int)
raw["full_matches"] = raw["full_matches"].astype(int)
raw["partial_matches"] = raw["partial_matches"].astype(int)

# Group by gene, mean_read_length, and coverage, then average full_matches and partial_matches across all iterations
averaged_data = raw.groupby(["gene", "mean_read_length", "coverage"], as_index=False)[["full_matches", "partial_matches"]].median()

print(averaged_data.head())
#%%
simulation_data = averaged_data.pivot_table(index=["gene"], columns="coverage", values="full_matches").reset_index()
simulation_data.index = simulation_data["gene"]
simulation_data = simulation_data.drop(columns=["gene"])
print(simulation_data.head())

simulation_data_partial = averaged_data.pivot_table(index=["gene"], columns="coverage", values="partial_matches").reset_index()
simulation_data_partial.index = simulation_data_partial["gene"]
simulation_data_partial = simulation_data_partial.drop(columns=["gene"])
print(simulation_data_partial.head())

#%%
"""
This is the code that works with multi-columned simulation input data
"""

# Load data (assuming two DataFrames: sequencing_data and simulation_data)

sequencing_data = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/full_matches_filtered.bed", sep="\t",  
                              usecols=[3,4], index_col=0, names=["gene", "count"])
sequencing_data_partial = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/partial_matches_filtered_final.bed", sep=" ",  
                              usecols=[3,4], index_col=0, names=["gene", "count"])

#simulation_data = pd.read_csv("/home/weichan/temporary/Data/Simulation/dummy_full_matches.bed", sep="\t",  
#                              usecols=[3,4], index_col=0, names=["gene", "count"])
#simulation_data_partial = pd.read_csv("/home/weichan/temporary/Data/Simulation/dummy_partial_matches_final.bed", sep="\t",  
#                              usecols=[3,4], index_col=0, names=["gene", "count"])

print(simulation_data_partial.head())

# Ensure both have the same order of genes
sequencing_data = sequencing_data.sort_index()
simulation_data = simulation_data.sort_index()
sequencing_data_partial = sequencing_data_partial.sort_index()
simulation_data_partial = simulation_data_partial.sort_index()


# combination of full and partial matches
#%% All vs all: i.e. also correlation between sim partial vs seq full 
# -> if more than one cov column is present, it shows the best correlation
combined_data = pd.concat([sequencing_data, simulation_data], axis=1, keys=["Seq", "Sim"]).fillna(0)
correlation_matrix = combined_data.corr()
print(correlation_matrix)

# Plot heatmap
plt.figure(figsize=(8,6))
sns.heatmap(correlation_matrix, 
            cmap="coolwarm", 
            center=0, 
            annot=True, 
            fmt=".2f", 
            cbar=True,
            vmin=-1.0, 
            vmax=1.0)
plt.title("Correlation Matrix for full matches")
plt.xlabel("")
plt.ylabel("")
plt.show()

# Compute covariate correlations
covariate_corrs = simulation_data.corrwith(sequencing_data["count"], axis=0).to_frame(name="Pearson Correlation")

# Display results
print(covariate_corrs)

# Normalize correlation values for consistent coloring
norm = plt.Normalize(vmin=-1.0, vmax=1.0)
colors = plt.cm.coolwarm(norm(covariate_corrs["Pearson Correlation"]))

# Plot barplot with heatmap colors
plt.figure(figsize=(8,6))
sns.barplot(x=covariate_corrs.index, y=covariate_corrs["Pearson Correlation"], palette=colors)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Pearson Correlation")
plt.title("Covariate-Level Correlation (Full)")
plt.show()

#%% Partial matches
print(sequencing_data_partial.head())

combined_data_partial = pd.concat([sequencing_data_partial, simulation_data_partial], axis=1, keys=["Seq", "Sim"]).fillna(0)
correlation_matrix_partial = combined_data_partial.corr(method="pearson")
print(correlation_matrix_partial)

# Plot heatmap or bar plot
plt.figure(figsize=(8,6))
sns.heatmap(correlation_matrix_partial, 
            cmap="coolwarm", 
            center=0, 
            annot=True, 
            fmt=".2f", 
            cbar=True,
            vmin=-1.0, 
            vmax=1.0)
plt.title("Correlation Matrix for partial matches")
plt.xlabel("")
plt.ylabel("")
plt.show()

# Compute covariate correlations
covariate_corrs = simulation_data_partial.corrwith(sequencing_data_partial["count"], axis=0).to_frame(name="Pearson Correlation")

# Display results
print(covariate_corrs)

# Normalize correlation values for consistent coloring
norm = plt.Normalize(vmin=-1.0, vmax=1.0)
colors = plt.cm.coolwarm(norm(covariate_corrs["Pearson Correlation"]))

# Plot barplot with heatmap colors
plt.figure(figsize=(8,6))
sns.barplot(x=covariate_corrs.index, y=covariate_corrs["Pearson Correlation"], palette=colors)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Pearson Correlation")
plt.title("Covariate-Level Correlation (Partial)")
plt.show()


#%% Band Altmann to assess whetehr there are udnerlyign biases

## if there is no pattern and mean difference close to 0: No bias. If most points within 1.96 SD, good. Others might be outliers influencing the result

full_combined = pd.concat([sequencing_data, simulation_data], axis=1).fillna(0)

mean_counts = (full_combined["count"].values + full_combined.loc[:,15].values) / 2 #18 is coverage
diff_counts = full_combined["count"].values - full_combined.loc[:,15].values

diff_counts = np.cbrt(diff_counts)
mean_counts = np.cbrt(mean_counts)

mean_diff = np.mean(diff_counts)
std_diff = np.std(diff_counts)

upper_limit = mean_diff + 1.96 * std_diff
lower_limit = mean_diff - 1.96 * std_diff

# Count the number of points within the range
within_limits = np.sum((diff_counts >= lower_limit) & (diff_counts <= upper_limit))
total_points = len(diff_counts)

print(within_limits)
print(total_points)
# Calculate the percentage of points within limits
percentage_within_limits = (within_limits / total_points) * 100

print(f"Percentage of points within ±1.96 SD: {percentage_within_limits:.2f}%")

# Identify points outside the limits
outside_points = np.where((diff_counts < lower_limit) | (diff_counts > upper_limit))[0]

# Plot Bland-Altman plot with labels for points outside the limits
plt.figure(figsize=(5, 5))
plt.scatter(mean_counts, diff_counts, alpha=0.5, edgecolors='k', linewidth=0.5, label='Data Points')
plt.axhline(mean_diff, color='red', linestyle='--', label='Mean Difference')
plt.axhline(mean_diff + 1.96 * std_diff, color='blue', linestyle='--', label='+1.96 SD')
plt.axhline(mean_diff - 1.96 * std_diff, color='blue', linestyle='--', label='-1.96 SD')

# Add labels for points outside the limits
for idx in outside_points:
    plt.text(mean_counts[idx], diff_counts[idx], full_combined.index[idx], fontsize=8, color='red')

plt.xlabel("Mean of Counts")
plt.ylabel("Difference (Seq - Sim)")
plt.title("Bland-Altman Plot (Full) (CBRT)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
sys.exit()
#%%
partial_combined = pd.concat([sequencing_data_partial, simulation_data_partial], axis=1).fillna(0)

print(partial_combined.tail())
partial_mean_counts = (partial_combined["count"].values + partial_combined.loc[:,15].values) / 2
partial_diff_counts = partial_combined["count"].values - partial_combined.loc[:,15].values

print(len(partial_diff_counts))
partial_diff_counts = np.cbrt(partial_diff_counts)
partial_mean_counts = np.cbrt(partial_mean_counts)


partial_mean_diff = np.mean(partial_diff_counts)
partial_std_diff = np.std(partial_diff_counts)

print(partial_mean_diff)

# Define upper and lower limits of agreement
upper_limit = partial_mean_diff + 1.96 * partial_std_diff
lower_limit = partial_mean_diff - 1.96 * partial_std_diff

# Count the number of points within the range
within_limits = np.sum((partial_diff_counts >= lower_limit) & (partial_diff_counts <= upper_limit))
total_points = len(partial_diff_counts)
print(within_limits)
print(total_points)

# Calculate the percentage of points within limits
percentage_within_limits = (within_limits / total_points) * 100
print(f"Percentage of points within ±1.96 SD: {percentage_within_limits:.2f}%")

# Identify points outside the limits
outside_points = np.where((partial_diff_counts < lower_limit) | (partial_diff_counts > upper_limit))[0]

# Plot Bland-Altman plot with labels for points outside the limits
plt.figure(figsize=(5, 5))
plt.scatter(partial_mean_counts, partial_diff_counts, alpha=0.5, edgecolors='k', linewidth=0.5, label='Data Points')
plt.axhline(partial_mean_diff, color='red', linestyle='--', label='Mean Difference')
plt.axhline(partial_mean_diff + 1.96 * partial_std_diff, color='blue', linestyle='--', label='+1.96 SD')
plt.axhline(partial_mean_diff - 1.96 * partial_std_diff, color='blue', linestyle='--', label='-1.96 SD')

# Add labels for points outside the limits
for idx in outside_points:
    plt.text(partial_mean_counts[idx], partial_diff_counts[idx], partial_combined.index[idx], fontsize=8, color='red')

plt.xlabel("Mean of Counts")
plt.ylabel("Difference (Seq - Sim)")
plt.title("Bland-Altman Plot (Partial) (CBRT)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

# Print the indices of points outside the limits
print("Indices of points outside the limits:", outside_points)
###
###
###
###
###

#%% lineplot with the columns partial

# looks like the simulation underestimates the partial matches

combined_data_partial = pd.concat([sequencing_data_partial, simulation_data_partial], axis=1)

print(combined_data_partial.head())
# Plot lineplot for each gene
plt.figure(figsize=(15, 10))
for n, gene in enumerate(combined_data_partial.index):
    # Plot the rest of the columns as a connected line
    plt.plot(combined_data_partial.columns.astype(str), np.log10(combined_data_partial.loc[gene]), marker='o', label=gene)

# Plot horizontal line for the first column
for n, gene in enumerate(combined_data_partial.index):
    # Plot horizontal line
    plt.hlines(y=np.log10(combined_data_partial.iloc[n, 0]), xmin=-0.25, xmax=len(combined_data_partial.columns[1:]) + 0.25, 
               colors='gray', linestyles='dashed', label=None, alpha=0.25)
    # Add gene label on the left side of the line
    plt.text(x=len(combined_data_partial.columns[1:]) + 0.5, y=np.log10(combined_data_partial.iloc[n, 0]), s=gene,
             ha='left', va='center', fontsize=10, alpha=0.7)

plt.xlabel("Coverage")
plt.ylabel("Log10 Counts")
plt.title("Lineplot of Partial Counts for Each Gene")
plt.xticks()
plt.legend(title="Gene", bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=8)
plt.tight_layout()
plt.show()

#%%
#%% lineplot with the columns full

combined_data_full = pd.concat([sequencing_data, simulation_data], axis=1)

print(combined_data_full.head())
# Plot lineplot for each gene
plt.figure(figsize=(15, 10))
for n, gene in enumerate(combined_data_full.index):
    # Plot the rest of the columns as a connected line
    plt.plot(combined_data_full.columns.astype(str), combined_data_full.loc[gene], marker='o', label=gene)

# Plot horizontal line for the first column
for n, gene in enumerate(combined_data_full.index):
    # Plot horizontal line
    plt.hlines(y=combined_data_full.iloc[n, 0], xmin=-0.25, xmax=len(combined_data_full.columns[1:]) + 0.25, 
               colors='gray', linestyles='dashed', label=None, alpha=0.25)
    # Add gene label on the left side of the line
    plt.text(x=len(combined_data_full.columns[1:]) + 0.5, y=combined_data_full.iloc[n, 0], s=gene,
             ha='left', va='center', fontsize=10, alpha=0.7)

plt.xlabel("Coverage")
plt.ylabel("Counts")
plt.title("Lineplot of Full Counts for Each Gene")
plt.xticks()
plt.legend(title="Gene", bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=8)
plt.tight_layout()
plt.show()


#%%

data = combined_data_full.reset_index().melt(id_vars='gene', var_name='Coverage', value_name='Full Matches')
# Static plot with seaborn
plt.figure(figsize=(20, 10))
sns.lineplot(data=data,
             x='gene', y='Full Matches', hue='Coverage', markers=True, dashes=False, linewidth=3, alpha=0.7,
             palette=['red' if cov == 'count' else 'grey' for cov in data['Coverage'].unique()])
plt.xlabel('Coverage')
plt.ylabel('Mean Full Matches')
plt.title('Lineplot of Full Matches by Coverage and Gene')
plt.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

data2 = combined_data_partial.reset_index().melt(id_vars='gene', var_name='Coverage', value_name='Partial Matches')
# Static plot with seaborn
plt.figure(figsize=(20, 10))
sns.lineplot(data=data2,
             x='gene', y='Partial Matches', hue='Coverage', markers=True, dashes=False, linewidth=3, alpha=0.7,
             palette=['red' if cov == 'count' else 'grey' for cov in data['Coverage'].unique()])
plt.xlabel('Coverage')
plt.ylabel('Mean Partial Matches')
plt.title('Lineplot of Partial Matches by Coverage and Gene')
plt.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()


#%% Concordance correlation coefficient: Measures agreement and precision (better than pearson)

# Function to calculate CCC
def concordance_correlation_coefficient(y_true, y_pred):
    """Computes Concordance Correlation Coefficient (CCC)."""
    mean_x = np.mean(y_true)
    mean_y = np.mean(y_pred)
    var_x = np.var(y_true)
    var_y = np.var(y_pred)
    covariance = np.cov(y_true, y_pred)[0, 1]

    ccc = (2 * covariance) / (var_x + var_y + (mean_x - mean_y) ** 2)
    return ccc

#%%
# Compute correlation matrix
# Iterate over each column in simulation_data_partial
# Create a grid of subplots for each column

print(simulation_data_partial)
num_columns = len(simulation_data_partial.columns)
fig, axes = plt.subplots(nrows=int(np.ceil(num_columns / 3)), ncols=3, figsize=(15, 5 * int(np.ceil(num_columns / 3))))
axes = axes.flatten()

for i, column in enumerate(simulation_data_partial.columns):
    x = simulation_data_partial[column].values  # Current column for x
    y = sequencing_data_partial["count"].values  # y remains the same

    # Compute CCC for the current column
    ccc_value = concordance_correlation_coefficient(x, y)
    print(f"Concordance Correlation Coefficient (CCC) for {column}: {round(ccc_value, 4)}")

    # Plot scatter and regression line in the corresponding subplot
    sns.regplot(ax=axes[i], x=x, y=y, scatter_kws={'alpha': 0.4}, line_kws={"color": "red"})
    axes[i].set_title(f"{column} (CCC={round(ccc_value, 3)})")
    axes[i].set_xlabel("Simulation Counts")
    axes[i].set_ylabel("Sequencing Counts")

    # Add the true line (1:1 line) in black
    min_val = min(min(x), min(y))
    max_val = max(max(x), max(y))
    axes[i].plot([min_val, max_val], [min_val, max_val], '--', color='black', label="1:1 Line")

# Remove any unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Adjust layout
plt.tight_layout()
plt.show()

# Compute correlation matrix full data
num_columns = len(simulation_data.columns)
fig, axes = plt.subplots(nrows=int(np.ceil(num_columns / 3)), ncols=3, figsize=(15, 5 * int(np.ceil(num_columns / 3))))
axes = axes.flatten()

for i, column in enumerate(simulation_data.columns):
    x = simulation_data[column].values  # Current column for x
    y = sequencing_data["count"].values  # y remains the same

    # Compute CCC for the current column
    ccc_value = concordance_correlation_coefficient(x, y)
    print(f"Concordance Correlation Coefficient (CCC) for {column}: {round(ccc_value, 4)}")

    # Plot scatter and regression line in the corresponding subplot
    sns.regplot(ax=axes[i], x=x, y=y, scatter_kws={'alpha': 0.4}, line_kws={"color": "red"})
    axes[i].set_title(f"{column} (CCC={round(ccc_value, 3)})")
    axes[i].set_xlabel("Simulation Counts")
    axes[i].set_ylabel("Sequencing Counts")

    # Add the true line (1:1 line) in black
    min_val = min(min(x), min(y))
    max_val = max(max(x), max(y))
    axes[i].plot([min_val, max_val], [min_val, max_val], '--', color='black', label="1:1 Line")

# Remove any unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Adjust layout
plt.tight_layout()
plt.show()


# %% bias
import ast

raw = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/Case1_bias/Case1_bias_matches_table.csv", sep="\t", header=0)
# Split using the last two underscores to handle problematic entries like 'GBA_x_0_999'
raw[["gene", "Barcode", "Iteration"]] = raw["Insertion"].str.rsplit("_", n=2, expand=True)
bias_list = []
bias_list.extend(raw["bias"].explode())
print(bias_list)
float_list = [val for sublist in bias_list for val in ast.literal_eval(sublist)]
#%%
print(len(float_list))
print(type(float_list))  # Should be <class 'list'>
print(type(float_list[0]))  # Should be <class 'float'>

# Create a seaborn histogram for all the elements in the bias list
plt.figure(figsize=(10, 6))
sns.histplot([x for x in float_list if -1 <= x <= 1], color="blue", alpha=0.7)
plt.title("Density Plot of Bias")
plt.xlabel("Bias")
plt.ylabel("Density")
plt.show()

# %%
roi = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/roi_info.csv", sep="\t", header=0)

print(roi.head())
plt.figure(figsize=(10, 6))
sns.histplot(roi, x=roi["Lenght"], color="blue", alpha=0.4)
plt.show()
# %%
# Ensure both have the same order of genes
sequencing_data = sequencing_data.sort_index()
simulation_data = simulation_data.sort_index()
sequencing_data_partial = sequencing_data_partial.sort_index()
simulation_data_partial = simulation_data_partial.sort_index()

print(sequencing_data.head())
print(simulation_data.head())
#merged = pd.merge(sequencing_data, simulation_data, left_index=True, right_index=True)
merged = pd.merge(sequencing_data_partial, simulation_data_partial, left_index=True, right_index=True)
#%%
# Create individual panels for each gene
num_genes = len(merged.index)
fig, axes = plt.subplots(nrows=int(np.ceil(num_genes / 5)), ncols=5, figsize=(15, 20))
axes = axes.flatten()

# Determine the global y-axis limits
y_min = merged.min().min()
y_max = merged.max().max()

for i, gene in enumerate(merged.index):
    ax = axes[i]
    gene_data = merged.loc[gene]
    # Plot the first value ("count") as a red dot
    ax.scatter(0, gene_data["count"], color="red", zorder=2, s=40)
    # Plot the other column values as black dots with alpha 0.5
    ax.scatter(range(1, len(gene_data)), gene_data.iloc[1:], color="black", alpha=0.5, s=40, zorder=1)
    # Connect the dots with a line
    ax.plot(range(len(gene_data)), gene_data, color="black", alpha=0.2, zorder=0)
    ax.set_title(gene)
    ax.set_xticks(range(len(gene_data)))
    ax.set_xticklabels(["count"] + list(gene_data.index[1:]), rotation=45, ha="right")
    ax.set_ylabel("Values")
    ax.set_ylim(y_min, y_max)  # Set the same y-axis limits for all subplots

# Remove any unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()

#%%
# Function to plot data for a specific gene
def plot_gene_data(gene_name):
    if gene_name not in merged.index:
        print(f"Gene '{gene_name}' not found in the data.")
        return
    
    gene_data = merged.loc[gene_name]
    plt.figure(figsize=(4, 4))
    # Plot the first value ("count") as a red dot
    plt.scatter(0, gene_data["count"], color="red", zorder=2, s=40, label="Sequencing Count")
    # Plot the other column values as black dots with alpha 0.5
    plt.scatter(range(1, len(gene_data)), gene_data.iloc[1:], color="black", alpha=0.5, s=40, label="Simulation Counts")
    # Connect the dots with a line
    plt.plot(range(len(gene_data)), gene_data, color="black", alpha=0.2, zorder=0)
    plt.title(f"Data for Gene: {gene_name}")
    plt.xticks(range(len(gene_data)), ["count"] + list(gene_data.index[1:]), rotation=45, ha="right")
    plt.ylabel("Values")
    plt.ylim(y_min, y_max)  # Use the same y-axis limits as the global plot
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_gene_data("MAPT")  # Example gene name
# %%
