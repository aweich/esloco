#!/bin/python

#%%

import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests



custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)
sns.set_context("talk", font_scale=1.0)

# Path to your log file
logfile = "/home/weichan/permanent/Projects/VIS/VIS_Simulation/out_rev/ressources/Case1_log.log"

# Regex patterns for the two line types
iteration_pattern = re.compile(
    r"Iteration\s+(\d+):\s+Time=(\d+(?:\.\d+)?)s,\s+Memory Used=(\d+(?:\.\d+)?)MB,\s+CPU=(\d+(?:\.\d+)?)%"
)
coverage_pattern = re.compile(r"Coverage:\s+(\d+)")

data = []
coverage = None

# Read and parse the file
with open(logfile) as f:
    for line in f:
        # Check for iteration info
        if "Iteration" in line:
            match = iteration_pattern.search(line)
            if match:
                iteration = int(match.group(1))
                time_s = float(match.group(2))
                memory_mb = float(match.group(3))
                cpu_percent = float(match.group(4))
                data.append({
                    "iteration": iteration,
                    "time_s": time_s,
                    "memory_mb": memory_mb,
                    "cpu_percent": cpu_percent,
                    "coverage": coverage  # will be filled when found
                })
        # Check for coverage info
        elif "Coverage: " in line:
            match = coverage_pattern.search(line)
            if match:
                coverage = int(match.group(1))

# Convert to DataFrame
df = pd.DataFrame(data)
#%%
# Optional: display summary
print(df.head())

print(df.tail())


# %%
print(df["time_s"])

fig, axes = plt.subplots(3, 1, figsize=(15, 8), sharex=True)

start = int(df["iteration"].min())
end = int(df["iteration"].max())
vlines = list(range(start, end + 1, 50))

# Time subplot
sns.lineplot(data=df, x="iteration", y="time_s", ax=axes[0])
for x in vlines:
    axes[0].axvline(x=x, color="gray", linestyle="--", alpha=0.5)
axes[0].set_xlabel("Iteration")
axes[0].set_ylabel("Time (s)")
axes[0].set_title("Time vs Iteration")
axes[0].set_xticks(vlines)

# Memory subplot
sns.lineplot(data=df, x="iteration", y="memory_mb", ax=axes[1])
for x in vlines:
    axes[1].axvline(x=x, color="gray", linestyle="--", alpha=0.5)
axes[1].set_xlabel("Iteration")
axes[1].set_ylabel("Memory (MB)")
axes[1].set_title("Memory vs Iteration")
axes[1].set_xticks(vlines)

sns.lineplot(data=df, x="iteration", y="cpu_percent", ax=axes[2])
for x in vlines:
    axes[2].axvline(x=x, color="gray", linestyle="--", alpha=0.5)
axes[2].set_xlabel("Iteration")
axes[2].set_ylabel("CPU core occupancy in %")
axes[2].set_title("Memory vs Iteration")
axes[2].set_xticks(vlines)


plt.tight_layout()
plt.show()
# %%
import numpy as np
print(f'Average time per iteration in h: {(np.mean(df["time_s"])/60)/60}')
# 266 min
#number of combinations per iteration

print(f'Max. memory consumed in an iteration: {df["memory_mb"].max()}')
print(df["memory_mb"])

#%%
from pathlib import Path
import seaborn as sns

data_dir = Path("/home/weichan/permanent/Projects/VIS/VIS_Simulation/out_rev/")
pattern = "*.npy"
#%%
records = []

for file in data_dir.glob(pattern):
    values = np.load(file).flatten()
    # store values in a list with filename label
    records.append(pd.DataFrame({"value": values, "dataset": file.stem}))
#%%

records = []

# === Load data ===
for file in data_dir.glob(pattern):
    data = np.load(file, allow_pickle=True)

    # If each entry is a list/array of lengths, expand those
    expanded = []
    for entry in data:
        if isinstance(entry, (list, np.ndarray)):
            expanded.extend(entry)
        else:
            expanded.append(entry)

    # Store values with filename label
    file_id = f"Sigma: {file.stem.split("_")[2]}"
    records.append(pd.DataFrame({"value": expanded, "dataset": file_id}))

# Combine
df = pd.concat(records, ignore_index=True)
print(df)
#%%
sns.set_context("talk", font_scale=0.8)
#df["logvalue"] =np.log10(df["value"])
# === Plot ===



plt.figure(figsize=(5, 5))
# Plot log-transformed values to highlight differences in sigma
sns.kdeplot(
    data=df,
    x=np.log(df["value"]),
    hue="dataset",
    common_norm=False,
    alpha=0.6,
    linewidth=4
)
sns.despine()
plt.xlabel("log(Read Length)")
plt.ylabel("Density")
#plt.legend(title="")
plt.title("Log-transformed distributions of Simulated Reads")
plt.tight_layout()
plt.show()
# %%

#speed benchmarking


data_dir = Path("/home/weichan/permanent/Projects/VIS/VIS_Simulation/out_rev/speed_benchmark/")
pattern = "*.log"

data = []
for file in data_dir.glob(pattern):
    with open(file) as f:
        for line in f:
            # Check for iteration info
            if "Total execution time:" in line:
                print(line)
                match = re.search(r"Total execution time:\s+(\d+(?:\.\d+)?)", line)
                print(match.group(1))
                if match:
                    # Extract number from dataset name
                    num_match = re.search(r"(\d+)", file.stem)
                    dataset_num = int(num_match.group(1)) if num_match else None
                    
                    data.append({
                        "execution_time_s": float(match.group(1)),
                        "dataset": file.stem,
                        "number": int(dataset_num)
                    })

df = pd.DataFrame(data)
#%%
print(df.tail())
# %%

df["label"] = df["dataset"].str.split("_").str[-2].astype("category")
df["category"] = df["label"].str.replace(r"(\d+)", "", regex=True)
#df["number"] = df["dataset"].str.extract(r'(\d+)').astype(int)  # Extract and convert to int
#%%
df.sort_values(["category", "number"], inplace=True)#, kind="stable")
print(df)

print(df.head())
print(type(df["category"].iloc[0]))
#%%

plt.figure(figsize=(5, 5))
sns.boxplot(data=df, 
            x="category", 
            y="execution_time_s", 
            hue="category",
            linewidth=2,

            showfliers=False)#, errorbar="sd", capsize=0.2)
sns.despine()
plt.xticks(rotation=90)
plt.xlabel("Category")
plt.ylabel("Category-level Execution Time (s)")

plt.figure(figsize=(12, 4))
sns.barplot(data=df, 
            x="dataset", 
            y="execution_time_s", 
            errorbar=None, 
            hue="category",
            width=0.75,
            linewidth=2, 
            edgecolor="black")
sns.despine()
plt.xticks(rotation=90)
plt.xlabel("Dataset Size")
# %%
import scipy.stats as stats
from sklearn.metrics import r2_score

categories = df["category"].unique()

# Determine subplot grid size automatically (2 columns)
n = len(categories)
ncols = 2
nrows = int(np.ceil(n / ncols))
df["log_t"] = np.log(df["execution_time_s"])

fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), sharey=True)
axes = axes.flatten()  # make indexing simpler
print(categories)
categories = ['barcodes','weight','coverage','iterations','masking','targetlengthkb']

for i, category in enumerate(categories):
    subset = df[df["category"] == category]
    ax = axes[i]

    # Scatter points
    ax.scatter(
        subset["number"], subset["execution_time_s"],
        s=80, edgecolor="black", linewidth=1.5, marker="o", alpha=0.8, color="white", zorder=5
    )

    x = subset["number"].values
    y = subset["execution_time_s"].values

    # --- Linear model ---
    slope, intercept, r_lin, p, se = stats.linregress(x, y)
    y_pred_lin = intercept + slope * x
    r2_lin = r2_score(y, y_pred_lin)

    # --- Exponential model: y = a * exp(bx) ---
    # (fit log(y) = log(a) + b*x)
    mask = y > 0  # avoid log(0)
    b, log_a = np.polyfit(x[mask], np.log(y[mask]), 1)
    a = np.exp(log_a)
    y_pred_exp = a * np.exp(b * x)
    r2_exp = r2_score(y[mask], y_pred_exp[mask])

    # --- Quadratic model: y = a*x^2 + b*x + c ---
    a2, b2, c2 = np.polyfit(x, y, 2)
    y_pred_quad = a2*x**2 + b2*x + c2
    r2_quad = r2_score(y, y_pred_quad)

    # --- Choose best model ---
    models = {
        "Linear": r2_lin,
        "Exponential": r2_exp,
        "Quadratic": r2_quad
    }
    best_model = max(models, key=models.get)

    print(f"\nModel comparison for category '{category}':")
    for name, r2 in models.items():
        print(f"  {name:10s} R² = {r2:.4f}")
    print(f"→ Best model: {best_model}")

    # --- Plot best model ---
    x_vals = np.linspace(x.min(), x.max(), 200)

    for model in models:
        if model == "Linear":
            y_vals = intercept + slope * x_vals
        elif model == "Exponential":
            y_vals = a * np.exp(b * x_vals)
        elif model == "Quadratic":
            y_vals = a2*x_vals**2 + b2*x_vals + c2

        ax.plot(x_vals, y_vals, linestyle="--", linewidth=2.5, alpha=0.7, label=f"{model} fit")
        ax.legend(frameon=False, fontsize=10)

    #x_vals = np.linspace(subset["number"].min(), subset["number"].max(), 100)
    #y_vals = intercept + slope * x_vals
   #ax.plot(x_vals, y_vals, linestyle="--", linewidth=2.5, alpha=0.7, color="C1")


    # Annotate regression equation
    #linear
    leq_text = f"Linear\n$y = {intercept:.2f} + {slope:.2f}x$\n$R^2 = {r2_lin:.3f}$"
    ax.text(
        1.05, 0.95, leq_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="round,pad=0.3")
        )
    #exponential y = a * exp(bx)
    eeq_text = f"Exponential\n$y = {a:.2f} \\cdot e^{{{b:.5f}x}}$\n$R^2 = {r2_exp:.3f}$"
    ax.text(
        1.05, 0.65, eeq_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="round,pad=0.3")
        )
    #a2*x**2 + b2*x + c2
    qeq_text = f"Quadratic\n$y = {a2:.5f}x² + {b2:.2f}x + {c2:.2f}$\n$R^2 = {r2_quad:.3f}$"
    ax.text(
        1.05, 0.35, qeq_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="round,pad=0.3")
        )

    ax.axvspan(xmin=ax.get_xlim()[0], xmax=50, 
               color="lightgrey", alpha=0.6, zorder=0)
    ax.set_title(category.upper(), fontstyle="oblique")
    ax.set_ylabel("Execution Time (s)")
    ax.set_xlabel("")
    sns.despine(ax=ax)

# Hide unused subplots if number of categories is odd
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
#plt.savefig("/home/weichan/temporary/Data/Simulation/RevisionPlots/SpeedBenchmark.svg", format="svg")
plt.show()


#%%
#create normalization for the soft masking

df = pd.read_csv("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/SRR8955270_coverage_100kb.bed", sep="\t", header=None, names=["chrom", "start", "end", "mean_cov"])

print(df.head())
# handle variable bin lengths (if needed)
df["length"] = df["end"] - df["start"]
equal_bins = df["length"].nunique() == 1
print(df.head())
df["mean_cov"] = pd.to_numeric(df["mean_cov"], errors='coerce')
print(len(df.index))

#%%

print(len(df[df["length"] != 100000].index))
print(len(df[df["length"] != 100000].index) / len(df.index) * 100) #only 1.4% of the bins are not 100k, difference can be ignored for now

eps = 1e-6      # pseudocount to avoid zeros; increase if your coverages are integers and small
alpha = 1.0     # power transform: alpha>1 sharpen, alpha<1 flatten
use_log_softmax = False  # set True to use log+softmax approach

# if any NaNs
df["mean_cov"] = df["mean_cov"].fillna(0.0)

# handle variable bin lengths (if needed)
df["length"] = df["end"] - df["start"]
equal_bins = df["length"].nunique() == 1

# Option A: weight proportional to mean_cov (equal bins)
if equal_bins:
    if not use_log_softmax:
        weights = (df["mean_cov"].values + eps) ** alpha
    else:
        scores = np.log(df["mean_cov"].values + eps)
        # beta controls sharpness, here use alpha as beta
        weights = np.exp(alpha * (scores - np.nanmax(scores)))
# Option B: if bins have unequal length -> convert to expected reads (mean_cov * length)
else:
    # total expected reads per bin (if mean_cov is per-base avg)
    expected_reads = (df["mean_cov"].values + eps) * df["length"].values
    if not use_log_softmax:
        weights = expected_reads ** alpha
    else:
        scores = np.log(expected_reads + eps)
        weights = np.exp(alpha * (scores - np.nanmax(scores)))

# Normalize to probabilities
probs = weights / np.sum(weights)
df["prob"] = probs

# Quick diagnostics
print("Total bins:", len(df))
print("Sum probs:", df["prob"].sum())
print("Top bins (by prob):")
print(df.sort_values("prob", ascending=False).head())

# Plot distribution of mean_cov and of prob
plt.figure(figsize=(12,4))
plt.subplot(1,2,1)
sns.histplot(np.log10(df["mean_cov"]+1), bins=100, kde=False)
plt.yscale("log")
plt.xlabel("Log10(mean_coverage+1)")
plt.title("Raw mean coverage distribution")

plt.subplot(1,2,2)
sns.histplot(np.log10(df["prob"]+1), bins=100, kde=False)
plt.yscale("log")
plt.xlabel("Sampling probability per bin")
plt.title("Normalized sampling probability")
plt.tight_layout()
plt.show()

# %%
def coverage_to_blocking(
    cov_array,
    method="linear",
    eps=1e-9,
    alpha=1.0,
    k=1.0,
    center=None,
    clip=True,
):
    """
    Map mean coverage values to blocking probabilities in [0,1] with cov==0 -> 1 and high cov -> 0.
    
    Parameters
    ----------
    cov_array : array-like
        Mean coverage values (non-negative).
    method : {"linear","power","log","sigmoid"}
        Mapping method.
    eps : float
        Small pseudocount to avoid log(0) / division-by-zero.
    alpha : float
        Power exponent used for "power" method (alpha>1 steepens).
    k : float
        Steepness parameter for sigmoid.
    center : float or None
        Midpoint for sigmoid. If None, set to median or mean of cov_array.
    clip : bool
        Clip output to [0,1].
    
    Returns
    -------
    blocking : np.ndarray
        Blocking probabilities (same shape as cov_array).
    """
    cov = np.asarray(cov_array, dtype=float)
    # ensure no negative values
    cov = np.where(np.isnan(cov), 0.0, cov)
    cov[cov < 0] = 0.0

    if method == "linear":
        # linear invert using min/max. If min is nonzero, we usually want min=0 so cov=0 -> 1.
        cov_min = 0.0  # force 0 -> maps to 1
        cov_max = cov.max() + eps
        scaled = (cov - cov_min) / (cov_max - cov_min)
        blocking = 1.0 - scaled

    elif method == "power":
        cov_max = cov.max() + eps
        scaled = (cov + eps) / (cov_max)
        blocking = 1.0 - (scaled ** alpha)

    elif method == "log":
        # log scaling - good if coverage spans orders of magnitude
        cov_log = np.log(cov + eps)
        cov_log_min = np.log(eps)
        cov_log_max = cov_log.max() if cov_log.max() > cov_log_min else cov_log_min + eps
        scaled = (cov_log - cov_log_min) / (cov_log_max - cov_log_min)
        blocking = 1.0 - scaled

    elif method == "sigmoid":
        if center is None:
            # set center to the 75th percentile by default (tune if you like)
            center = np.percentile(cov, 50)
        s = 1.0 / (1.0 + np.exp(-k * (cov - center)))
        blocking = 1.0 - s

    else:
        raise ValueError("Unknown method: choose 'linear','power','log','sigmoid'")

    if clip:
        blocking = np.clip(blocking, 0.0, 1.0)

    return blocking


methods = ["linear", "power", "log", "sigmoid"]
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for i, method in enumerate(methods):
    if method == "power":
        b = coverage_to_blocking(df["mean_cov"], method=method, alpha=0.25)
    elif method == "sigmoid":
        b = coverage_to_blocking(df["mean_cov"], method=method, k=0.2, center=0.25)
    else:
        b = coverage_to_blocking(df["mean_cov"], method=method)
    
    axes[i].scatter(df["mean_cov"], b, marker='o', linewidth=0, alpha=0.8)
    axes[i].set_xscale("log")
    axes[i].set_xlabel("Mean coverage")
    axes[i].set_ylabel("Blocking probability")
    axes[i].set_title(f"{method.capitalize()} method")
    axes[i].grid(True, alpha=0.3)

plt.suptitle("Mean Coverage to Weighted Blocking (methods)", fontsize=14, y=1.00)
plt.tight_layout()
plt.show()

# %%
def coverage_to_blocking_linear(cov):
    cov = np.asarray(cov, dtype=float)
    cov = np.clip(cov, 0, None)  # prevent negatives
    cov_max = cov.max()
    if cov_max == 0:
        return np.ones_like(cov)
    blocking = 1 - (cov / cov_max)
    return blocking


blocking = coverage_to_blocking_linear(np.log10(df["mean_cov"]+1))
blocking = np.round(blocking , decimals=2)
#%%
plt.figure(figsize=(6,4))
sns.scatterplot(x=df["mean_cov"], y=blocking, marker='o', size=120, sizes=(120,120), alpha=0.3)
plt.xlabel("log10(Mean coverage + 1)")
plt.ylabel("Blocking probability")
plt.xscale("log")
plt.legend([],[], frameon=False)
plt.title("Linear inverse mapping: 0→1, max→0")
#plt.grid(True)
plt.show()
# %%
sns.histplot(blocking, bins=50, kde=True)
plt.yscale("log")
print(blocking)

#%%
def coverage_to_blocking_sigmoid(cov, k=None, x0=None):
    cov = np.asarray(cov, dtype=float)
    cov = np.clip(cov, 0, None)
    if x0 is None:
        x0 = np.median(cov)
    if k is None:
        k = 0.5 / x0  # moderate nonlinearity

    blocking = 1 / (1 + np.exp(k * (cov - x0)))
    return blocking


blocking_sig = coverage_to_blocking_sigmoid(df["mean_cov"], k=0.05)
blocking_sig = np.round(blocking_sig, decimals=2)

plt.figure(figsize=(6,4))
sns.scatterplot(x=df["mean_cov"], y=blocking_sig, marker='o', size=120, sizes=(120,120), alpha=0.3)
plt.xlabel("")
plt.xscale("log")
plt.ylabel("Blocking probability")
plt.legend([],[], frameon=False)
plt.title("Linear inverse mapping: 0→1, max→0")
#plt.grid(True)
plt.show()


# %%



# Plot distribution of mean_cov and of prob
plt.figure(figsize=(12,4))
plt.subplot(1,3,1)
sns.histplot(np.log10(df["mean_cov"]+1), bins=100, kde=False)
plt.yscale("log")
plt.xlabel("Log10(mean_coverage+1)")
plt.title("Raw mean coverage distribution")

plt.subplot(1,3,2)#, sharey=plt.gca().axes)
sns.histplot(blocking, bins=50, kde=True)
#plt.yscale("log")
plt.xlabel("Masking probability")
plt.title("Linear Soft masking")
plt.tight_layout()


plt.subplot(1,3,3)#, sharey=plt.gca().axes)
sns.histplot(blocking_sig, bins=50, kde=True)
#plt.yscale("log")
plt.xlabel("Masking probability")
plt.title("Sigmoidal soft masking")
plt.tight_layout()
plt.show()


#%%

print(blocking_sig)

df["sigmoidal_blocking"] = np.round(blocking_sig, decimals=2)
df["linear_blocking"] = np.round(blocking , decimals=2)

print(df.head())

#%%
#file export

df_sig = df[["chrom", "start", "end", "sigmoidal_blocking"]]
#df_sig.to_csv("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/Case1_sigmoidal_mask.bed", sep="\t", header=False, index=False)

df_lin = df[["chrom", "start", "end", "linear_blocking"]]
#df_lin.to_csv("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/Case1_linear_mask.bed", sep="\t", header=False, index=False)
# %%


# %%

# Plotting the influence of different number sof iterations

raw = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test_reproduction/Case1/Case1_matches_table.csv", sep="\t", header=0)

raw["iteration"] = raw["iteration"].astype(int)

seq = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test_reproduction/fromsra/SRR8955270_GroundTruth_OTBs.bed", sep="\t", header=None)
print(seq.head())
seq["Length"] = seq[2] - seq[1]

#%%
genelengths = seq[[3, "Length"]].drop_duplicates()
genelengths = genelengths.rename(columns={3: "gene"})
print(genelengths.head())
raw = pd.merge(raw, genelengths, how="left", left_on="target", right_on="gene")
print(raw.head())
raw["otb_norm"] = raw["on_target_bases"] / raw["Length"]
print(raw.head())

# %%

collection = []

for target in raw["target"].unique():
    subset= raw[raw["target"] == target]
    subset = subset[["mean_read_length", "coverage", "on_target_bases", "iteration", "target", "Length"]]

    subset["iteration_group"] = pd.cut(
        subset["iteration"],
        bins=[0, 10, 100, np.inf],
        labels=["<10", "10–100", ">100"]
    )

    results = subset.groupby(["target", "iteration_group", "coverage", "Length"], as_index=False).agg(
        mean_bases_on_target=('on_target_bases', 'mean'),
        sem_bases_on_target=('on_target_bases', 'sem')
    )

    results["norm_SEM"] = results["sem_bases_on_target"] / results["Length"]
    results["norm_mean"] = results["mean_bases_on_target"] / results["Length"]

    collection.append(results)

#%%
import random 
combined = pd.concat(collection, ignore_index=True)

combined = combined[combined["coverage"] == 26]

# Optional: sort by length for smoother lines
combined = combined.sort_values("Length")

# Create the plot
plt.figure(figsize=(12, 6))
sns.lineplot(
    data=combined,
    x="Length",
    y="norm_SEM",
    hue="iteration_group",
    style="iteration_group",
    markers=True,
    dashes=False,
    linewidth=3,
    alpha=0.7,
    palette=["red", "orange", "blue"],
    markeredgecolor="black",
    markeredgewidth=1.5,
    markersize=10,
    markerfacecolor="white" 
    
)

# Add labels to the points
print(combined["target"].unique())
for i in combined["target"].unique():
    print(combined["Length"].loc[combined["target"] == i].unique()[0])
    plt.text(combined["Length"].loc[combined["target"] == i].unique()[0], 
             2.1,
             i, 
             fontsize=11, ha="center", va="bottom", rotation=45)
plt.xscale("log")
plt.title("Effect of Target Length on SEM across Iteration Groups (Coverage = 27x)", pad=40)
plt.xlabel("Target Length (bp)")
plt.ylabel("SEM (Target length normalized)")
plt.legend(title="Iteration Group")
sns.despine()
plt.tight_layout()
#plt.savefig("/home/weichan/temporary/Data/Simulation/RevisionPlots/TargetLength_SEM.svg", format="svg", bbox_inches='tight')
plt.show()


#%% comparison PacBio vs ONT OTBs

pb = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test_reproduction/fromsra/SRR8955270_GroundTruth_OTBs.bed", sep="\t", header=None)
print(pb.head())
pb["Length"] = pb[2] - pb[1]
pb = pb.groupby([3, "Length"])[10].sum().reset_index()

pb.columns = [3, "Length", "bases_on_target"]
pb["norm_OTB"] = pb["bases_on_target"] / pb["Length"]
print(pb.head())
#%%
ont = pd.read_csv("/home/weichan/temporary/Data/Simulation/ONT_benchmark/ONT_GroundTruth_OTBs.bed", sep="\t", header=None)
ont["Length"] = ont[2] - ont[1]
ont = ont.groupby([3, "Length"])[10].sum().reset_index()

ont.columns = [3, "Length", "bases_on_target"]
ont["norm_OTB"] = ont["bases_on_target"] / ont["Length"]
print(ont.head())

#%%
#combine the data

both = pd.merge(pb, ont, how="inner", on=[3, "Length"], suffixes=("_pb", "_ont"))
print(both.head())

#%%
raw = pd.read_csv("/home/weichan/temporary/Data/Simulation/ONT_benchmark/out_full/Case1_ONT_matches_table.csv", sep="\t", header=0)

print(raw.head())

# Split using the last two underscores to handle problematic entries like 'GBA_x_0_999'
raw[["gene", "Barcode", "Iteration"]] = raw["target_region"].str.rsplit("_", n=2, expand=True)

raw = raw[(raw["coverage"] == 55)]
raw = raw.drop(columns=["target_region", "Barcode"])
print(raw.head())

raw["Iteration"] = raw["Iteration"].astype(int)
raw["full_matches"] = raw["full_matches"].astype(int)
raw["partial_matches"] = raw["partial_matches"].astype(int)

# %%

#repeat with PB data

pbsim = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test_reproduction/Case1/Case1_matches_table.csv", sep="\t", header=0)
#pbsim = pd.read_csv("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/Case1_logistic_matches_table.csv", sep="\t", header=0)
pbsim[["gene", "Barcode", "Iteration"]] = pbsim["target_region"].str.rsplit("_", n=2, expand=True)

pbsim = pbsim[(pbsim["coverage"] == 26)]
pbsim = pbsim.drop(columns=["target_region", "Barcode"])

pbsim["Iteration"] = pbsim["Iteration"].astype(int)
pbsim["full_matches"] = pbsim["full_matches"].astype(int)
pbsim["partial_matches"] = pbsim["partial_matches"].astype(int)
pbsim["pb_otbs"] = pbsim["on_target_bases"]
print(pbsim.head())


#%%
print(raw.head())
#%%
# ont merge with both
both = both.rename(columns={3: "target"})

# Merge simulation data (ONT simulation)
combined = pd.merge(raw, both, on="target", how="inner")
combined["norm_OTB_ontsim"] = (combined["on_target_bases"] / combined["Length"]) / 54 #coverage

# Merge PacBio simulation data
pbsim = pbsim.rename(columns={"target": "target"})  # ensure consistent naming
combined = pd.merge(pbsim, combined, on="target", how="inner")
combined["norm_OTB_pbsim"] = (combined["pb_otbs"] / combined["Length"]) / 26

print("Merged Data:")
print(combined.head())
#%%
plot_df = pd.melt(
    combined,
    id_vars=["target"],
    value_vars=["norm_OTB_ontsim", "norm_OTB_pbsim"],
    var_name="Platform",
    value_name="Normalized_OTB"
)

# Map readable names
platform_labels = {
    "norm_OTB_ontsim": "ONT",
    "norm_OTB_pbsim": "PB"
}
plot_df["Platform"] = plot_df["Platform"].map(platform_labels)


#%%

results = []  # store dicts with gene, platform, empirical_p, N_sim, etc.

# compute empirical p-values
for platform in plot_df["Platform"].unique():
    df_plat = plot_df[plot_df["Platform"] == platform]
    #N_sims_per_gene = df_plat.groupby("target")["Normalized_OTB"].count().to_dict()

    for gene in df_plat["target"].unique():
        gene_data = df_plat[df_plat["target"] == gene]["Normalized_OTB"].values
        N = len(gene_data)
        if N == 0:
            continue

        # Fetch the corresponding true value (make sure you divide correctly)
        if platform == "ONT":
            true_val = both.loc[both["target"] == gene, "norm_OTB_ont"].values
            denom = 54.0
        else:
            true_val = both.loc[both["target"] == gene, "norm_OTB_pb"].values
            denom = 26.0

        if true_val.size == 0:
            continue

        true_v = float(true_val[0]) / denom  # keep your scaling if desired

        center = np.median(gene_data)

        # two-sided empirical p (distance from center)
        distances = np.abs(gene_data - center)
        obs_distance = abs(true_v - center)
        more_extreme = np.sum(distances >= obs_distance)

        empirical_p = more_extreme/ N

        results.append({
            "gene": gene,
            "platform": platform,
            "empirical_p": empirical_p,
            "N": N,
            "true_v": true_v,
            "sim_center": center,
            "more_extreme": int(more_extreme)
        })

# convert to DataFrame
res_df = pd.DataFrame(results)

# Apply BH correction per platform
res_df["adjusted_p"] = np.nan
for platform in res_df["platform"].unique():
    mask = res_df["platform"] == platform
    pvals = res_df.loc[mask, "empirical_p"].values
    if len(pvals) == 0:
        continue
    _, adj_p, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    res_df.loc[mask, "adjusted_p"] = adj_p

# optional: also compute global correction across all tests
_, adj_p_global, _, _ = multipletests(res_df["empirical_p"].values, alpha=0.05, method="fdr_bh")
res_df["adjusted_p_global"] = adj_p_global

print(res_df.sort_values(["platform","adjusted_p"]).head(20))

print(res_df[res_df["adjusted_p"] <= 0.05].groupby("platform").size())

#%%
gene_order = sorted(plot_df["target"].unique().tolist())

# Convert target columns to ordered categorical
plot_df["target"] = pd.Categorical(plot_df["target"], categories=gene_order, ordered=True)
both["target"] = pd.Categorical(both["target"], categories=gene_order, ordered=True)
sns.set_context("talk", font_scale=1.1)
plt.figure(figsize=(18, 6))

sns.boxplot(
    data=plot_df,
    x="target",
    y="Normalized_OTB",
    hue="Platform",
    dodge=True,
    palette={"ONT": "#00A2FF", "PB": "#FF0080"},
    showfliers=False,
    linecolor="black",
    width=0.8,
    fill=False,
    orient="v",
    linewidth=2,
    #saturation=0.5,
    #native_scale=True,
    gap=0,
    order=gene_order  # explicitly specify order
)

# Get numeric positions for categorical x-axis
x_positions = np.arange(len(gene_order))
offset = 0.25  # adjust to align with boxplot position

# Overlay mean reference lines from `both`
sns.scatterplot(
    data=both,
    x=x_positions - offset,
    y=both["norm_OTB_ont"]/54,
    linewidth=2,
    edgecolor="#000000",
    #s=120,
    alpha=1,
    marker="o",
    label="ONT (Sequencing)",
    color="#00A2FF",
    zorder=3
    #sort=False  # preserve order
)


sns.scatterplot(
    data=both,
    x=x_positions + offset,  # shift right to align with PB boxplot
    y=both["norm_OTB_pb"]/26,
    linewidth=2,
    alpha=1,
    edgecolor="#000000",
    marker="o",
    label="PacBio (Sequencing)",
    color="#FF0080",
    zorder=3
)


# choose which adjusted p to show (platform-level)
for _, row in res_df.iterrows():
    gene = row["gene"]
    plat = row["platform"]
    adjp = row["adjusted_p"]
    x_pos = gene_order.index(gene)
    p_value = round(adjp, 3)

    # choose y coordinate carefully (above plotted points)
    y_pos = 1.6 # tweak; or compute per-gene max + margin

    if adjp <= 0.001:
        label = "***"
    elif adjp <= 0.01:
        label = "**"
    elif adjp <= 0.05:
        label = "*"
    else:
        label = "n.s."
    
    if label != "n.s.":
        plt.text(x_pos + (0.25 if plat == "PB" else -0.25), y_pos, f"{label}, {p_value} ({plat})",ha="center", va="bottom", fontsize=14, rotation=90)


plt.xticks(rotation=90)
plt.xlabel("Target Region")
plt.ylabel("Normalized On-Target Bases")
plt.legend(title="Platform", frameon=False)
sns.despine()
plt.tight_layout()
plt.show()

# %%

# plot different distributions
import matplotlib.pyplot as plt

def load_lengths(path, max_reads=1_000_000):
    # Stream file and stop after max_reads entries
    lens = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= max_reads:
                break
            lens.append(float(line.strip()))
    return np.array(lens)

lens1 = load_lengths("/home/weichan/temporary/Data/Simulation/ONT_benchmark/ONT_lengths.txt")
lens2 = load_lengths("/home/weichan/temporary/Data/Simulation/ROI_test_reproduction/fromsra/PB_lengths.txt")

# --- Combine ---
df_lens = pd.DataFrame({
    "length": np.concatenate([lens1, lens2]),
    "Platform": ["ONT"] * len(lens1) + ["PB"] * len(lens2)
})

# --- Filter invalids ---
df_lens = df_lens[df_lens["length"] > 0]

# --- Plot histogram (much faster than KDE) ---
plt.figure(figsize=(3, 3))
sns.histplot(
    data=df_lens,
    x="length",
    hue="Platform",
    palette={"ONT": "#00A2FF", "PB": "#FF0080"},
    log_scale=True,      # logarithmic x-axis
    element="step",      # outline hist for clarity
    stat="density",      # comparable to KDE
    common_norm=False,
    bins=200,            # adjust granularity
    alpha=0.75
)

plt.xlabel("")
plt.ylabel("")
plt.yticks([0, 1])
plt.title("")
plt.legend([],[], frameon=False)
sns.despine()
plt.tight_layout()
#plt.savefig("/home/weichan/temporary/Data/Simulation/RevisionPlots/Distributions.png", format="png", bbox_inches='tight', dpi=500)
plt.show()
# %%
print(combined.head())
#ccc plot for specific coverage
def concordance_correlation_coefficient(y_true, y_pred):
    """Computes Concordance Correlation Coefficient (CCC) with consistent ddof=0."""
    mean_x = np.mean(y_true)
    mean_y = np.mean(y_pred)
    var_x = np.var(y_true)  # ddof=0 by default
    var_y = np.var(y_pred)
    covariance = np.mean((y_true - mean_x) * (y_pred - mean_y))  # Manual population covariance

    ccc = (2 * covariance) / (var_x + var_y + (mean_x - mean_y) ** 2)
    return ccc


#%% ONT CCC
sns.set_context("talk", font_scale=0.8)

ccc_data = combined[["target", "on_target_bases_y", "iteration_y"]]
ccc_data = ccc_data.groupby("target").agg({"on_target_bases_y": "mean"}).reset_index()
ccc_data.index = ccc_data["target"]
ccc_data = ccc_data.drop(columns=["target"])
print(ccc_data.head())

# Ensure the column exists or rename it if necessary
# Align indices before concatenation

# Set index on both DataFrames to align by target/gene
both_indexed = both.set_index("target")
ccc_data_indexed = ccc_data.reset_index().set_index("target")

# Concatenate along columns (axis=1) with aligned indices

ccc_data_indexed = ccc_data_indexed.rename(columns={"on_target_bases_y": "sim"})
ccc_data = pd.concat([ccc_data_indexed, both_indexed["bases_on_target_ont"]], axis=1)
print(ccc_data.head())

ccc_data["seq"] = ccc_data["bases_on_target_ont"]
ccc_data = ccc_data.drop(columns=["bases_on_target_ont"])
print(ccc_data.head())

fig, axes = plt.subplots(2, 1, figsize=(4, 8))
x = ccc_data["sim"].values
y = ccc_data["seq"].values

# Compute CCC
ccc_value = concordance_correlation_coefficient(x, y)
print(f"Concordance Correlation Coefficient (CCC) for 54: {ccc_value:.3f}")

# Left plot: Regression
sns.regplot(
    ax=axes[0],
    x=x,
    y=y,
    color='#00A2FF',
    scatter_kws={'alpha': 0.7, 's': 90, "zorder":10, "edgecolor":"black", "linewidths":1.5},
    line_kws={"color": "black", "linestyle": "-"},
    label='ONT Data Points',
)
axes[0].set_title(f"54x (CCC={ccc_value:.4f})")
axes[0].set_xlabel("ONT Simulation OTBs")
axes[0].set_ylabel("ONT Sequencing OTBs")
axes[0].legend()

# Right plot: Bland-Altman
mean_counts = (y + x) / 2 
diff_counts = y - x

diff_counts = np.cbrt(diff_counts)
mean_counts = np.cbrt(mean_counts)

mean_diff = np.mean(diff_counts)
std_diff = np.std(diff_counts)

upper_limit = mean_diff + 1.96 * std_diff
lower_limit = mean_diff - 1.96 * std_diff

within_limits = np.sum((diff_counts >= lower_limit) & (diff_counts <= upper_limit))
total_points = len(diff_counts)
percentage_within_limits = (within_limits / total_points) * 100

print(f"Percentage of points within ±1.96 SD: {percentage_within_limits:.2f}%")

outside_points = np.where((diff_counts < lower_limit) | (diff_counts > upper_limit))[0]

axes[1].scatter(mean_counts, diff_counts, alpha=0.7, edgecolors='k', s=90, linewidth=1.5, label='Data Points', color='#00A2FF')
axes[1].axhline(mean_diff, color='red', linestyle='--', label='Mean Difference', linewidth=3)
axes[1].axhline(upper_limit, color='black', linestyle='--', label='+1.96 SD', linewidth=3)
axes[1].axhline(lower_limit, color='black', linestyle='--', label='-1.96 SD', linewidth=3)

for idx in outside_points:
    axes[1].text(mean_counts[idx], diff_counts[idx], ccc_data.index[idx], fontsize=10, color='red')

axes[1].set_xlabel("Mean of OTBs")
axes[1].set_ylabel("Difference (Seq - Sim)")
axes[1].set_title("Bland-Altman Plot (CBRT)")
axes[1].legend()

plt.tight_layout()
plt.show()
#%%
cov26 = combined[combined["coverage_x"] == 26]
ccc_data = cov26[["target", "on_target_bases_x", "iteration_x"]]

ccc_data = ccc_data.groupby("target").agg({"on_target_bases_x": "mean"}).reset_index()
ccc_data.index = ccc_data["target"]
ccc_data = ccc_data.drop(columns=["target"])
print(ccc_data.head())

# Ensure the column exists or rename it if necessary
# Align indices before concatenation

# Set index on both DataFrames to align by target/gene
both_indexed = both.set_index("target")
ccc_data_indexed = ccc_data.reset_index().set_index("target")

# Concatenate along columns (axis=1) with aligned indices

ccc_data_indexed = ccc_data_indexed.rename(columns={"on_target_bases_x": "sim"})
ccc_data = pd.concat([ccc_data_indexed, both_indexed["bases_on_target_pb"]], axis=1)
print(ccc_data.head())

ccc_data["seq"] = ccc_data["bases_on_target_pb"]
ccc_data = ccc_data.drop(columns=["bases_on_target_pb"])
print(ccc_data.head())

fig, axes = plt.subplots(2, 1, figsize=(4, 8))
x = ccc_data["sim"].values
y = ccc_data["seq"].values

# Compute CCC
ccc_value = concordance_correlation_coefficient(x, y)
print(f"Concordance Correlation Coefficient (CCC) for 26: {ccc_value:.3f}")

# Left plot: Regression
sns.regplot(
    ax=axes[0],
    x=x,
    y=y,
    color='#FF0080',
    scatter_kws={'alpha': 0.7, 's': 90, "zorder":10, "edgecolor":"black", "linewidths":1.5},
    line_kws={"color": "black", "linestyle": "-"},
    label='PB Data Points',
)
axes[0].set_title(f"26x (CCC={ccc_value:.4f})")
axes[0].set_xlabel("PB Simulation OTBs")
axes[0].set_ylabel("PB Sequencing OTBs")
axes[0].legend()

# Right plot: Bland-Altman
mean_counts = (y + x) / 2 
diff_counts = y - x

diff_counts = np.cbrt(diff_counts)
mean_counts = np.cbrt(mean_counts)

mean_diff = np.mean(diff_counts)
std_diff = np.std(diff_counts)

upper_limit = mean_diff + 1.96 * std_diff
lower_limit = mean_diff - 1.96 * std_diff

within_limits = np.sum((diff_counts >= lower_limit) & (diff_counts <= upper_limit))
total_points = len(diff_counts)
percentage_within_limits = (within_limits / total_points) * 100

print(f"Percentage of points within ±1.96 SD: {percentage_within_limits:.2f}%")

outside_points = np.where((diff_counts < lower_limit) | (diff_counts > upper_limit))[0]

axes[1].scatter(mean_counts, diff_counts, alpha=0.7, edgecolors='k', s=90, linewidth=1.5, label='Data Points', color='#FF0080')
axes[1].axhline(mean_diff, color='red', linestyle='--', label='Mean Difference', linewidth=3)
axes[1].axhline(upper_limit, color='black', linestyle='--', label='+1.96 SD', linewidth=3)
axes[1].axhline(lower_limit, color='black', linestyle='--', label='-1.96 SD', linewidth=3)

for idx in outside_points:
    axes[1].text(mean_counts[idx], diff_counts[idx], ccc_data.index[idx], fontsize=10, color='red')

axes[1].set_xlabel("Mean of OTBs")
axes[1].set_ylabel("Difference (Seq - Sim)")
axes[1].set_title("Bland-Altman Plot (CBRT)")
axes[1].legend()

plt.tight_layout()
plt.show()
# %%
import pyranges as pr

# normalization for the last question

genes = pd.read_csv("/home/weichan/temporary/Data/Simulation/ROI_test_reproduction/fromsra/SRR8955270_GroundTruth_OTBs.bed", sep="\t", header=None)
print(genes.head())

cov = pd.read_csv("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/SRR8955270_coverage_100kb.bed", sep="\t", header=None, names=["chrom", "start", "end", "mean_cov"])
cov["mean_cov"] = pd.to_numeric(cov["mean_cov"], errors='coerce')
cov["mean_cov"] = cov["mean_cov"].fillna(0.0)
# Coverage dataframe
cov_df = cov.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End", "mean_cov": "mean_cov"})

# Gene dataframe — keep only relevant columns (chrom, start, end, gene)
genes_df = genes[[0, 5, 6, 3]].copy()
genes_df.columns = ["Chromosome", "Start", "End", "Gene"]


print(cov_df.head())
print(genes_df.head())

pr_cov = pr.PyRanges(cov_df)
pr_genes = pr.PyRanges(genes_df)

print(type(pr_cov)) #pandas and not pyyRanges 
# Find overlaps
overlaps = pr_cov.overlap(pr_genes)

print(type(overlaps.head()))
test = overlaps.copy()
print(type(test))
# %%
overlap_idx = cov_df.index[cov_df.set_index(["Chromosome", "Start", "End"]).index.isin(
    overlaps[["Chromosome", "Start", "End"]].set_axis(["Chromosome", "Start", "End"], axis=1).apply(tuple, axis=1)
)]
plt.figure(figsize=(14, 5))

# Base coverage line
plt.plot(cov_df["Start"], np.log(cov_df["mean_cov"]), color="lightgray", linewidth=1.5, label="Coverage")

# Highlight overlapping bins
plt.scatter(
    y=cov_df.loc[overlap_idx, "Start"],
    x=np.log(cov_df.loc[overlap_idx, "mean_cov"]),
    color="red",
    s=20,
    label="Panel Genes"
)

plt.xlabel("Genomic position (bp)")
plt.ylabel("Mean coverage")
plt.title("Coverage profile with gene overlaps")
plt.legend(frameon=False)
plt.tight_layout()
plt.show()
# %%
def coverage_to_blocking_linear(cov):
    cov = np.asarray(cov, dtype=float)
    cov = np.clip(cov, 0, None)  # prevent negatives
    cov_max = cov.max()
    if cov_max == 0:
        return np.ones_like(cov)
    blocking = 1 - (cov / cov_max)
    return blocking


transformed = np.log(cov_df["mean_cov"]+1)
# Scale by median
transformed = transformed / np.median(transformed)
# Z-score normalization
transformed = (transformed - np.mean(transformed)) / np.std(transformed)
cov_df["transformed"] = transformed
blocking = coverage_to_blocking_linear(transformed)
blocking = np.round(blocking , decimals=2)

print(cov_df.head())
#%%
plt.figure(figsize=(6,4))
sns.scatterplot(x=cov_df["transformed"], y=blocking, marker='o', size=120, sizes=(120,120), alpha=0.3)
sns.scatterplot(
    y=blocking[overlap_idx],
    x=cov_df.loc[overlap_idx, "transformed"],
    color="red",
    s=100,
    label="Panel Genes"
)
plt.xlabel("log10(Mean coverage + 1)")
plt.ylabel("Blocking probability")
#plt.xscale("log")
#plt.legend([],[], frameon=False)
plt.title("Linear inverse mapping: 0→1, max→0")
#plt.grid(True)
plt.show()
# %% 

from scipy.stats import rankdata
# new strategies
from scipy.stats.mstats import winsorize
def winsorize_series(s, lower_pct, upper_pct):
    lo = s.quantile(lower_pct)
    hi = s.quantile(upper_pct)
    return s.clip(lower=lo, upper=hi)


# 3) Logistic mapping for smooth middle behavior (after robust scaling)
def map_block_prob_logistic(cov_df, log_offset=1, lower_pct=0.1, upper_pct=0.9, k=10, mid=0.6):
    cov = cov_df['mean_cov'].astype(float).copy()
    cov_w = winsorize_series(cov, lower_pct, upper_pct)
    #cov_w = winsorize(cov, limits=[lower_pct, None])
    logcov = np.log10(cov_w)# + log_offset)
    # scale to [0,1]
    lo, hi = np.percentile(logcov, 1), np.percentile(logcov, 99)
    s = (logcov - lo) / (hi - lo)
    s = np.clip(s, 0, 1)
    # logistic: p(block) = 1 / (1 + exp(k*(s - mid)))
    # large k = sharper threshold, smaller k = smoother curve
    block_prob = 1.0 / (1.0 + np.exp(k * (s - mid)))
    return block_prob

cov_df['block_prob_logistic'] = map_block_prob_logistic(cov_df)

cov_df['m04'] = map_block_prob_logistic(cov_df, mid=0.4)
cov_df['m05'] = map_block_prob_logistic(cov_df, mid=0.5)
cov_df['m06'] = map_block_prob_logistic(cov_df, mid=0.6)

print(cov_df.head())
cov_df["label"] = cov_df["Chromosome"].astype(str) + ":" + cov_df["Start"].astype(str)

#plt.figure(figsize=(5, 8))
sns.lineplot(x=cov_df["mean_cov"], y=cov_df["block_prob_logistic"], label="Logistic", color="red")
for i in overlap_idx:
    plt.text(cov_df.loc[i, "mean_cov"], cov_df.loc[i, "block_prob_logistic"], 
             cov_df.loc[i, "label"], fontsize=8, ha='left', rotation=90)

# Color chr1 dots in blue, others in red
chr1_mask = (cov_df.loc[overlap_idx, "Chromosome"] == "chr1") & (cov_df.loc[overlap_idx, "Start"] >= 155000000)
chr1_idx = overlap_idx[chr1_mask]
other_idx = overlap_idx[~chr1_mask]
other_idx = overlap_idx[~chr1_mask]

# Plot chr1 dots in blue
if len(chr1_idx) > 0:
    sns.scatterplot(
        y=cov_df["block_prob_logistic"][chr1_idx],
        x=cov_df.loc[chr1_idx, "mean_cov"],
        color="blue",
        s=100,
        label="Panel Genes (chr1)",
        zorder=3
    )

# Plot other chromosomes in red
if len(other_idx) > 0:
    sns.scatterplot(
        y=cov_df["block_prob_logistic"][other_idx],
        x=cov_df.loc[other_idx, "mean_cov"],
        color="red",
        s=100,
        label="Panel Genes (other)"
    )

plt.xscale("log")

#%%
df_logistic = cov_df[["Chromosome", "Start", "End", "block_prob_logistic"]]
#df_logistic.to_csv("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/Case1_logistic90k10m06_mask.bed", sep="\t", header=False, index=False)

# %%

maskpath = Path("/home/weichan/temporary/Data/Simulation/mask_test/")
plt.figure(figsize=(4,4))
for file in maskpath.glob("Case1_logistic_90_k10*.bed"):
    maskfile = pd.read_csv(file, sep="\t")
    print(maskfile.head())
    case1, logistic, cutoff, k, m, mask = file.stem.split("_")
    label = f"{m}"
    
    # Color mapping based on cutoff value
    color_map = {"m04": "blue", "m05": "red", "m06": "orange"}
    color = color_map.get(m, "black")  # default to black if cutoff not in map
    
    sns.lineplot(x=cov_df["mean_cov"], y=maskfile["weight"], label=label, linewidth=4, alpha=0.6, color=color)
    sns.scatterplot(
    y=cov_df[m][other_idx],
    x=cov_df.loc[other_idx, "mean_cov"],
    color=color,
    s=60,
    alpha=0.8,
    edgecolor="black",
    label=""
)
    plt.xscale("log")
plt.xlabel("Mean coverage (binned 100kb)")
plt.ylabel("Blocking probability / mask weight")
#plt.savefig("/home/weichan/temporary/Data/Simulation/RevisionPlots/SoftMask_Lineplot.svg", format="svg")
plt.show()

# %%

'''outpath = Path("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/logistic/refined/")
for file in outpath.glob("*matches_table.csv"):
    out = pd.read_csv(file, sep="\t", header=0)
    case1, logistic, cutoff, k, m, matches, table = file.stem.split("_")
    label = f"{cutoff}_{k}_{m}"

    color_map = {"m04": "blue", "m05": "red", "m06": "orange"}
    color = color_map.get(m, "black")  # default to black if cutoff not in map
    
    
    out[["gene", "Barcode", "Iteration"]] = out["target_region"].str.rsplit("_", n=2, expand=True)

    out = out.drop(columns=["target_region", "Barcode"])

    out["Iteration"] = out["Iteration"].astype(int)
    out[f"{label}_otbs"] = out["on_target_bases"]

    out = out.rename(columns={"target": "target"})  # ensure consistent naming
    combined = pd.merge(out, both, on="target", how="inner")
    combined[f"norm_{label}_OTB"] = (combined[f"{label}_otbs"] / combined["Length"]) / 26

    #
    results = []
    for gene in combined["target"].unique():
        gene_data = combined[combined["target"] == gene][f"norm_{label}_OTB"].values
        N = len(gene_data)
        if N == 0:
            continue

        true_val = both.loc[both["target"] == gene, "norm_OTB_pb"].values
        denom = 26.0

        if true_val.size == 0:
            continue

        true_v = float(true_val[0]) / denom  # keep your scaling if desired

        center = np.median(gene_data)

        # two-sided empirical p (distance from center)
        distances = np.abs(gene_data - center)
        obs_distance = abs(true_v - center)
        more_extreme = np.sum(distances >= obs_distance)

        empirical_p = more_extreme/ N

        results.append({
            "gene": gene,
            "empirical_p": empirical_p,
            "N": N,
            "true_v": true_v,
            "sim_center": center,
            "more_extreme": int(more_extreme)
        })
    res_df = pd.DataFrame(results)

    # Apply BH correction per platform
    res_df["adjusted_p"] = np.nan
    _, adj_p, _, _ = multipletests(res_df["empirical_p"], alpha=0.05, method="fdr_bh")
    res_df["adjusted_p"] = adj_p
    print(len(res_df[res_df["adjusted_p"] <= 0.05]))
    
    
    gene_order = sorted(combined["target"].unique().tolist())

    # Convert target columns to ordered categorical
    combined["target"] = pd.Categorical(combined["target"], categories=gene_order, ordered=True)
    both["target"] = pd.Categorical(both["target"], categories=gene_order, ordered=True)
    sns.set_context("talk", font_scale=1.1)
    plt.figure(figsize=(18, 6))

    sns.boxplot(
        data=combined,
        x="target",
        y=f"norm_{label}_OTB",
        showfliers=False,
        linecolor="black",
        width=0.8,
        fill=False,
        color=color,
        orient="v",
        linewidth=2,
        #saturation=0.5,
        #native_scale=True,
        gap=0,
        order=gene_order  # explicitly specify order
    )

    # Get numeric positions for categorical x-axis
    x_positions = np.arange(len(gene_order))
    offset = 0.0 # adjust to align with boxplot position

    sns.lineplot(
        data=both,
        x=x_positions + offset,  # shift right to align with PB boxplot
        y=both["norm_OTB_pb"]/26,
        linewidth=2,
        alpha=1,
        #edgecolor="#000000",
        marker="o",
        #label="PacBio (Sequencing)",
        color="#000000",
        zorder=3
    )

    for _, row in res_df.iterrows():
        gene = row["gene"]
        adjp = row["adjusted_p"]
        x_pos = gene_order.index(gene)
        p_value = round(adjp, 3)

        # choose y coordinate carefully (above plotted points)
        y_pos = 1.6 # tweak; or compute per-gene max + margin

        if adjp <= 0.001:
            label = "***"
        elif adjp <= 0.01:
            label = "**"
        elif adjp <= 0.05:
            label = "*"
        else:
            label = "n.s."
        
        #if label != "n.s.":
        plt.text(x_pos, y_pos, f"{label}, {p_value}",ha="center", va="bottom", fontsize=14, rotation=90)

    plt.text(x_positions[0], 0.25, f"Significant: {len(res_df[res_df['adjusted_p'] <= 0.05])}", fontsize=24, weight='bold')
    plt.xticks(rotation=90)
    #plt.title()
    plt.xlabel("Target Region")
    plt.ylabel("Normalized On-Target Bases")
    plt.legend(frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.show()
'''
# %%
outpath = Path("/home/weichan/temporary/Data/Simulation/Revision_SoftMasking_Case1/logistic/refined/")

# accumulate simulation values and per-gene empirical p-results
long_records = []
pvals_results = []

for file in outpath.glob("*matches_table.csv"):
    out = pd.read_csv(file, sep="\t", header=0)
    case1, logistic, cutoff, k, m, matches, table = file.stem.split("_")
    label = f"{cutoff}_{k}_{m}"  # mask id used as hue in the final plot

    out[["gene", "Barcode", "Iteration"]] = out["target_region"].str.rsplit("_", n=2, expand=True)
    out = out.drop(columns=["target_region", "Barcode"])
    out["Iteration"] = out["Iteration"].astype(int)
    out = out.rename(columns={"target": "target"})  # keep naming consistent

    # merge with sequencing truth table `both`
    merged = pd.merge(out, both, on="target", how="inner")
    merged[f"norm_{label}_OTB"] = (merged["on_target_bases"] / merged["Length"]) / 26.0

    # record long-format simulation values for plotting
    tmp = merged[["target", f"norm_{label}_OTB"]].copy()
    tmp = tmp.rename(columns={f"norm_{label}_OTB": "norm_OTB"})
    tmp["mask"] = label
    long_records.append(tmp)

    # compute empirical p-values per gene (same approach as original loop)
    for gene in merged["target"].unique():
        gene_data = merged[merged["target"] == gene][f"norm_{label}_OTB"].values
        N = len(gene_data)
        if N == 0:
            continue

        true_val = both.loc[both["target"] == gene, "norm_OTB_pb"].values
        if true_val.size == 0:
            continue

        true_v = float(true_val[0]) / 26.0
        center = np.median(gene_data)
        distances = np.abs(gene_data - center)
        obs_distance = abs(true_v - center)
        more_extreme = np.sum(distances >= obs_distance)
        empirical_p = more_extreme / N

        pvals_results.append({
            "gene": gene,
            "mask": label,
            "m":m,
            "empirical_p": empirical_p,
            "N": N
        })

# make combined DataFrames
sim_df = pd.concat(long_records, ignore_index=True)
pvals_df = pd.DataFrame(pvals_results)

# BH correction per mask
pvals_df["adjusted_p"] = np.nan
for mask_label, group in pvals_df.groupby("mask"):
    pvals = group["empirical_p"].values
    if len(pvals) == 0:
        continue
    _, adj, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    pvals_df.loc[pvals_df["mask"] == mask_label, "adjusted_p"] = adj

# Prepare plotting order and palette
gene_order = sorted(sim_df["target"].unique().tolist())
sim_df["target"] = pd.Categorical(sim_df["target"], categories=gene_order, ordered=True)
both["target"] = pd.Categorical(both["target"], categories=gene_order, ordered=True)

masks = sorted(sim_df["mask"].unique().tolist())
nm = len(masks)
palette = {'90_k10_m04': "blue", '90_k10_m05': "red", '90_k10_m06': "orange"}

# compute max y per gene to place annotations
max_per_gene = sim_df.groupby("target")["norm_OTB"].max()
# consider also PB truth values when computing annotation height
pb_vals = both.set_index("target")["norm_OTB_pb"] / 26.0
max_y = pd.concat([max_per_gene, pb_vals.rename("pb")], axis=1).max(axis=1)
y_margin = (max_y.max() - max_y.min()) * 0.08 if max_y.max() != max_y.min() else 0.05

plt.figure(figsize=(20, 6))
sns.boxplot(
    data=sim_df,
    x="target",
    y="norm_OTB",
    hue="mask",
    palette=palette,
    showfliers=False,
    linewidth=2,
    fill=False,
    width=0.8
)

# overlay PB sequencing reference points
x_positions = np.arange(len(gene_order))
sns.lineplot(
    data=both,
    x=x_positions,  # shift right to align with PB boxplot
    y=both["norm_OTB_pb"]/26,
    linewidth=2,
    alpha=1,
    #edgecolor="#000000",
    marker="o",
    #label="PacBio (Sequencing)",
    color="#000000",
    zorder=3
)

# annotate significant mask/gene pairs with stars
# compute horizontal offsets for masks so stars don't overlap
if nm > 1:
    offsets = np.linspace(-0.25, 0.25, nm)
else:
    offsets = [0.0]
mask_offset = {mask: offsets[i] for i, mask in enumerate(masks)}

for _, row in pvals_df.iterrows():
    adjp = row["adjusted_p"]
    if np.isnan(adjp):
        continue
    if adjp <= 0.001:
        star = "***"
    elif adjp <= 0.01:
        star = "**"
    elif adjp <= 0.05:
        star = "*"
    else:
        continue  # only annotate significant tests

    gene = row["gene"]
    mask = row["mask"]
    m = row["m"]
    x_pos = gene_order.index(gene) + mask_offset[mask]
    y_pos = max_y.loc[gene] + y_margin
    plt.text(x_pos, y_pos, f"{star} {adjp} ({m})", ha="center", va="bottom", fontsize=12, color="black", rotation=90)
     
plt.xticks(rotation=90)
plt.xlabel("Target Region")
plt.ylabel("Normalized On-Target Bases")
plt.legend(frameon=False)
sns.despine()
plt.tight_layout()
#plt.savefig("/home/weichan/temporary/Data/Simulation/RevisionPlots/SoftMask_Boxplots.svg", format="svg")
plt.show()
#%%
plt.figure(figsize=(4, 4))
sns.barplot(y=["Unmasked", "m04", "m05", "m06"], x=[8, 6,4,4], palette={"Unmasked": "lightgrey", "m04": "blue", "m05": "red", "m06": "orange"},
            edgecolor="black", linewidth=2)
plt.xticks([0,2,4,6,8,10])
plt.xlabel("Number of significant genes")
#plt.savefig("/home/weichan/temporary/Data/Simulation/RevisionPlots/SoftMask_Barplots.svg", format="svg")


#%%
import pysam

bam_path = "/home/weichan/temporary/Data/Simulation/ROI_test_reproduction/fromsra/SRR8955270.bam"

# Region of interest
chrom = "chr1"
# use 155,000,000-range coordinates to match GBA1 location
start = 155200000
end   = 155300000

bam = pysam.AlignmentFile(bam_path, "rb")

# returns four arrays (A,C,G,T counts)
cov = bam.count_coverage(contig=chrom, start=int(start), stop=int(end), quality_threshold=0)

print(cov)

# total coverage per base
total_cov = np.sum(cov, axis=0)
print(total_cov)
# genomic coordinate array
positions = np.arange(start, end)

#%%155234452	155244699
gba1_start = 155234452
gba1_end = 155244699

# plot coverage and highlight GBA1 using absolute genomic coordinates (same scale as `positions`)
plt.figure(figsize=(16,4))
plt.plot(positions, total_cov, linewidth=2, color="grey", label="Read Coverage")
plt.axvspan(gba1_start, gba1_end, color="red", alpha=0.3, label="GBA1 gene region")
plt.xlim(start, end)
plt.xlabel("Genomic position")
plt.ylabel("Read Coverage")
#plt.title("Coverage over region")
plt.legend(frameon=False)
plt.tight_layout()
#plt.savefig("/home/weichan/temporary/Data/Simulation/RevisionPlots/SoftMask_GBA1Coverage.svg", format="svg")

plt.show()

# %%
