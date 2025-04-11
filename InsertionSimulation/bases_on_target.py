
## bases on target calculation
#%%
import pandas as pd
from scipy.stats import pearsonr
import sys
import seaborn as sns
import numpy as np
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

raw = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/Case1_bias/Case1_bias4_matches_table.csv", sep="\t", header=0)
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

# %%

# Sum up the integers in the "overlap" column
raw["bases_on_target"] = raw["overlap"].apply(lambda x: sum(eval(x)) if isinstance(x, str) else sum(x))

# Calculate bases on target for each iteration by dividing by the "Length" column
#raw["bases_on_target_x"] = raw["bases_on_target"] / raw["Length"]

# Group by gene, mean_read_length, and coverage, then calculate the median of bases_on_target_perc across iterations
# Calculate the median bases_on_target_perc for each gene, mean_read_length, and coverage
averaged_bases_on_target = raw.groupby(["gene", "mean_read_length", "coverage"], as_index=False)["bases_on_target"].median()

print(averaged_bases_on_target.head())

# Generate individual plots for each coverage setting
coverage_levels = raw["coverage"].unique()

print(len(averaged_bases_on_target))

# %%
seq = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/seq_bases_on_target.bed", sep="\t", header=None)
print(seq.head())
seq["Length"] = seq[2] - seq[1]
true_positive = seq.groupby([3, "Length"])[10].sum().reset_index()
print(len(true_positive.index))
true_positive.columns = [3, "Length", "bases_on_target"]
true_positive["bases_on_target_x"] = true_positive["bases_on_target"] / true_positive["Length"]
print(true_positive.head())

# %%


statistical_results = []
for coverage in coverage_levels:
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=raw[raw["coverage"] == coverage], x="gene", y="bases_on_target", palette="Set3")
    sns.lineplot(data=true_positive, x=3, y="bases_on_target", color="red", label="True Positive", marker="o", alpha=0.5)
    plt.xlabel("Gene")
    plt.ylabel("Bases on Target (X Coverage)")
    plt.xticks(rotation=90)

    empirical_p_values = []
    temp_results = []

    for gene in raw["gene"].unique():
        gene_data = raw[(raw["coverage"] == coverage) & (raw["gene"] == gene)]["bases_on_target"]
        true_positive_value = true_positive[true_positive[3] == gene]["bases_on_target"].values

        if true_positive_value.size > 0:
            sim_mean = gene_data.mean()
            distance = abs(true_positive_value[0] - sim_mean)

            # Calculate empirical p-value (two-sided)
            more_extreme = np.sum(np.abs(gene_data - sim_mean) >= distance)
            empirical_p = more_extreme / len(gene_data)

            temp_results.append({"gene": gene, "coverage": coverage, "empirical_p": empirical_p})
            empirical_p_values.append(empirical_p)

    # Apply Benjamini-Hochberg correction
    _, corrected_pvals, _, _ = multipletests(empirical_p_values, alpha=0.05, method="fdr_bh")

    # Append corrected results
    significant_count = 0
    for i, result in enumerate(temp_results):
        result["adjusted_p"] = corrected_pvals[i]
        statistical_results.append(result)

        if corrected_pvals[i] <= 0.05:
            significant_count += 1

        x_position = raw["gene"].unique().tolist().index(result["gene"])
        y_position = 4e7  # Customize for your y-axis

        # Round and label
        p_value = round(result["adjusted_p"], 3)
        if p_value > 0.05:
            label = "n.s."
        elif 0.01 < p_value <= 0.05:
            label = "*"
        else:
            label = "**"

        plt.text(x_position, y_position, f"{label}, {p_value}", ha="center", va="bottom", fontsize=6, rotation=90)

    print(f"Coverage {coverage}: {significant_count} significant values")

    # Legend
    '''
    plt.legend(handles=[
        plt.Line2D([0], [0], color='none', label="Significance Levels:"),
        plt.Line2D([0], [0], color='none', label="n.s.: p > 0.05"),
        plt.Line2D([0], [0], color='none', label="*: p ≤ 0.05"),
        plt.Line2D([0], [0], color='none', label="**: p ≤ 0.01")
    ], loc="upper right", frameon=False)
    '''

    plt.tight_layout()
    plt.show()


# %%
def concordance_correlation_coefficient(y_true, y_pred):
    """Computes Concordance Correlation Coefficient (CCC) with consistent ddof=0."""
    mean_x = np.mean(y_true)
    mean_y = np.mean(y_pred)
    var_x = np.var(y_true)  # ddof=0 by default
    var_y = np.var(y_pred)
    covariance = np.mean((y_true - mean_x) * (y_pred - mean_y))  # Manual population covariance

    ccc = (2 * covariance) / (var_x + var_y + (mean_x - mean_y) ** 2)
    return ccc
#%%

ccc_data = averaged_bases_on_target.pivot_table(index=["gene"], columns="coverage", values="bases_on_target").reset_index()
ccc_data.index = ccc_data["gene"]
ccc_data = ccc_data.drop(columns=["gene"])
print(len(ccc_data.index))
print(len(true_positive.index))

# Ensure the column exists or rename it if necessary
# Align indices before concatenation
true_positive_aligned = true_positive.set_index(ccc_data.index)
ccc_data = pd.concat([ccc_data, true_positive_aligned["bases_on_target"]], axis=1)
ccc_data["seq. count"] = ccc_data["bases_on_target"]
ccc_data = ccc_data.drop(columns=["bases_on_target"])
print(ccc_data.head())

### ccc
num_columns = len(ccc_data.columns)
fig, axes = plt.subplots(nrows=int(np.ceil(num_columns / 2)), ncols=6, figsize=(20, 5 * int(np.ceil(num_columns / 3))))
axes = axes.flatten()

for i, column in enumerate(ccc_data.columns):
    x = ccc_data[column].values  # Current column for x
    y = ccc_data["seq. count"].values  # y remains the same

    # Compute CCC for the current column
    ccc_value = concordance_correlation_coefficient(x, y)
    print(f"Concordance Correlation Coefficient (CCC) for {column}: {ccc_value:.3f}")

    # Plot scatter and regression line in the corresponding subplot
    sns.regplot(ax=axes[i], x=x, y=y, scatter_kws={'alpha': 0.4}, line_kws={"color": "red"})
    axes[i].set_title(f"{column} (CCC={ccc_value:.3f})")
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
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.show()

# %%
x = np.array([1, 2, 3, 4, 5])
print("CCC self test:", concordance_correlation_coefficient(x, x))
# %%
#%% Bland altmann plot

full_combined = ccc_data 

mean_counts = (full_combined["seq. count"].values + full_combined.loc[:,26].values) / 2 
diff_counts = full_combined["seq. count"].values - full_combined.loc[:,26].values

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
plt.title("Bland-Altman Plot (CBRT)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
# %%
roi = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/roi_info.csv", sep="\t", header=0)

print(roi.head())
plt.figure(figsize=(10, 6))
sns.histplot(roi, x=roi["Lenght"], color="blue", alpha=0.4, bins=45)
plt.show()

# %%
seqreads = 6763580 #manually from file
reads = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test/Case1_bias/Case1_bias4_barcode_distribution_table.csv", sep="\t", header=0)
print(reads.head())
# Plot the barplot with error bars based on the standard deviation
plt.figure(figsize=(5, 5))
sns.catplot(data=reads, x="coverage", y=reads["Total_Reads"], kind="bar", errorbar="sd", palette="grey")
plt.axhline(seqreads, color='red', linewidth=3, alpha=0.6, linestyle='--', label='Number of Sequencing Reads')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Coverage")
plt.ylabel("Mean Simulated Reads")
plt.show()
# %%
