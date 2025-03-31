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
raw = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/test_Case1/Case1_matches_only.csv", sep="\t", header=0)
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

sequencing_data = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/full_matches.bed", sep="\t",  
                              usecols=[3,4], index_col=0, names=["gene", "count"])
sequencing_data_partial = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/partial_matches_final.bed", sep="\t",  
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

#%%
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

sys.exit()

#%% Band Altmann to assess whetehr there are udnerlyign biases

## if there is no pattern and mean difference close to 0: No bias. If most points within 1.96 SD, good. Others might be outliers influencing the result

full_combined = pd.concat([sequencing_data, simulation_data], axis=1).fillna(0)

mean_counts = (full_combined["count"].values + full_combined.iloc[:,-1].values) / 2
diff_counts = full_combined["count"].values - full_combined.iloc[:,-1].values

diff_counts = np.cbrt(diff_counts)
mean_counts = np.cbrt(mean_counts)

mean_diff = np.mean(diff_counts)
std_diff = np.std(diff_counts)

plt.figure(figsize=(5, 5))
plt.scatter(mean_counts, diff_counts, alpha=0.5, edgecolors='k', linewidth=0.5, label='Data Points')
plt.axhline(mean_diff, color='red', linestyle='--', label='Mean Difference')
plt.axhline(mean_diff + 1.96 * std_diff, color='blue', linestyle='--', label='+1.96 SD')
plt.axhline(mean_diff - 1.96 * std_diff, color='blue', linestyle='--', label='-1.96 SD')
plt.xlabel("Mean of Counts")
plt.ylabel("Difference (Seq - Sim)")
plt.title("Bland-Altman Plot (Full) (CBRT)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

#%%
partial_combined = pd.concat([sequencing_data_partial, simulation_data_partial], axis=1).fillna(0)
print(partial_combined.tail())
partial_mean_counts = (partial_combined["count"].values + partial_combined.iloc[:,-1].values) / 2
partial_diff_counts = partial_combined["count"].values - partial_combined.iloc[:,-1].values

print(len(partial_diff_counts))
partial_diff_counts = np.cbrt(partial_diff_counts)
partial_mean_counts = np.cbrt(partial_mean_counts)


partial_mean_diff = np.mean(partial_diff_counts)
partial_std_diff = np.std(partial_diff_counts)

print(partial_mean_diff)

plt.figure(figsize=(5, 5))
plt.scatter(partial_mean_counts, partial_diff_counts, alpha=0.5, edgecolors='k', linewidth=0.5, label='Data Points')
plt.axhline(partial_mean_diff, color='red', linestyle='--', label='Mean Difference')
plt.axhline(partial_mean_diff + 1.96 * std_diff, color='blue', linestyle='--', label='+1.96 SD')
plt.axhline(partial_mean_diff - 1.96 * std_diff, color='blue', linestyle='--', label='-1.96 SD')
plt.xlabel("Mean of Counts")
plt.ylabel("Difference (Seq - Sim)")
plt.title("Bland-Altman Plot (Partial) (CBRT)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

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
    plt.plot(range(len(combined_data_partial.columns)), np.log10(combined_data_partial.loc[gene]), marker='o', label=gene)

# Plot horizontal line for the first column
for n, gene in enumerate(combined_data_partial.index):
    # Plot horizontal line
    plt.hlines(y=np.log10(combined_data_partial.iloc[n, 0]), xmin=-0.25, xmax=len(combined_data_partial.columns[1:]) + 0.25, 
               colors='gray', linestyles='dashed', label=None, alpha=0.25)
    # Add gene label on the left side of the line
    plt.text(x=len(combined_data_partial.columns[1:]) + 0.5, y=np.log10(combined_data_partial.iloc[n, 0]), s=gene,
             ha='left', va='center', fontsize=10, alpha=0.7)

plt.xlabel("Coverage / 2")
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
    plt.plot(range(len(combined_data_full.columns)), combined_data_full.loc[gene], marker='o', label=gene)

# Plot horizontal line for the first column
for n, gene in enumerate(combined_data_full.index):
    # Plot horizontal line
    plt.hlines(y=combined_data_full.iloc[n, 0], xmin=-0.25, xmax=len(combined_data_full.columns[1:]) + 0.25, 
               colors='gray', linestyles='dashed', label=None, alpha=0.25)
    # Add gene label on the left side of the line
    plt.text(x=len(combined_data_full.columns[1:]) + 0.5, y=combined_data_full.iloc[n, 0], s=gene,
             ha='left', va='center', fontsize=10, alpha=0.7)

plt.xlabel("Coverage / 2")
plt.ylabel("Counts")
plt.title("Lineplot of Full Counts for Each Gene")
plt.xticks()
plt.legend(title="Gene", bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=8)
plt.tight_layout()
plt.show()

# %%
"""
This is the code, if the input data only has one column, i.e. if one specific condition is present that needs to be compared against the True Positive
"""

sequencing_data = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/full_matches.bed", sep="\t",  
                              usecols=[3,4], index_col=0, names=["gene", "count"])
sequencing_data_partial = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/partial_matches_final.bed", sep="\t",  
                              usecols=[3,4], index_col=0, names=["gene", "count"])


simulation_data = pd.read_csv("/home/weichan/temporary/Data/Simulation/dummy_full_matches.bed", sep="\t",  
                              usecols=[3,4], index_col=0, names=["gene", "count"])
simulation_data_partial = pd.read_csv("/home/weichan/temporary/Data/Simulation/dummy_partial_matches_final.bed", sep="\t",  
                              usecols=[3,4], index_col=0, names=["gene", "count"])

# Ensure both have the same order of genes
sequencing_data = sequencing_data.sort_index()
simulation_data = simulation_data.sort_index()
sequencing_data_partial = sequencing_data_partial.sort_index()
simulation_data_partial = simulation_data_partial.sort_index()

print(simulation_data_partial.head())
#%%
# Compute correlation
pearson_corr, _ = pearsonr(sequencing_data.values.flatten(), simulation_data.values.flatten())
pearson_corr_partial, _ = pearsonr(sequencing_data_partial.values.flatten(), simulation_data_partial.values.flatten())

print(f"Pearson Correlation (Full): {pearson_corr:.3f}")
print(f"Pearson Correlation (Partial): {pearson_corr_partial:.3f}")

#%%
plt.figure(figsize=(6,6))
sns.regplot(x=sequencing_data.values.flatten(), y=simulation_data.values.flatten(), scatter_kws={'alpha': 0.5}, line_kws={"color": "red"})
plt.xlabel("Sequencing Counts")
plt.ylabel("Simulation Counts")
plt.title(f"Pearson Correlation (Full): {pearson_corr:.3f}")
plt.show()

plt.figure(figsize=(6,6))
sns.regplot(x=sequencing_data_partial.values.flatten(), y=simulation_data_partial.values.flatten(), scatter_kws={'alpha': 0.5}, line_kws={"color": "red"})
plt.xlabel("Sequencing Counts")
plt.ylabel("Simulation Counts")
plt.title(f"Pearson Correlation (Partial): {pearson_corr_partial:.3f}")
plt.show()
# %%
