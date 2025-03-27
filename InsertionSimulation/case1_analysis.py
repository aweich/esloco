#%%
import pandas as pd
from scipy.stats import pearsonr
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)
#%%
"""
This is the code that works with multi-columned simulation input data
"""

# Load data (assuming two DataFrames: sequencing_data and simulation_data)

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

#%% # dummy data: manually add more columns for now
for i in range(1, 6):
    simulation_data[f"cov{i}"] = np.random.randint(0, 10, size=len(simulation_data))
print(simulation_data.head())

for i in range(1, 6):
    simulation_data_partial[f"cov{i}"] = np.random.randint(0, 1500, size=len(simulation_data_partial))
print(simulation_data_partial.head())

#%%

# combination of full and partial matches



#%% All vs all: i.e. also correlation between sim partial vs seq full 
# -> if more than one cov column is present, it shows the best correlation
combined_data = pd.concat([sequencing_data, simulation_data], axis=1, keys=["Seq", "Sim"])
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

combined_data_partial = pd.concat([sequencing_data_partial, simulation_data_partial], axis=1, keys=["Seq", "Sim"])
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
