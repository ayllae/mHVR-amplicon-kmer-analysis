#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list
from natsort import natsorted


input_file = "cluster_relative_frequencies_per_sample.csv"
output_file = "heatmap_braycurtis_blue_annotated.png"
cluster_samples = True   # True = hierarchical clustering

# ----------------------------------------------------------
# Load data
# ----------------------------------------------------------
data = pd.read_csv(input_file, index_col=0)
data = data.apply(pd.to_numeric)

# Normalize rows to sum to 1 (important for proportions)
data = data.div(data.sum(axis=1), axis=0)

# ----------------------------------------------------------
# Bray–Curtis dissimilarity
# ----------------------------------------------------------
dist_array = squareform(pdist(data.values, metric="braycurtis"))

dist_df = pd.DataFrame(
    dist_array,
    index=data.index,
    columns=data.index
)

# ----------------------------------------------------------
# Ordering samples
# ----------------------------------------------------------
if cluster_samples:
    linkage_matrix = linkage(
        pdist(data.values, metric="braycurtis"),
        method="average"
    )
    ordered_indices = leaves_list(linkage_matrix)
    ordered_labels = dist_df.index[ordered_indices]
else:
    ordered_labels = natsorted(dist_df.index)

dist_df = dist_df.loc[ordered_labels, ordered_labels]

# ----------------------------------------------------------
# Plot heatmap (single blue gradient)
# ----------------------------------------------------------
n = dist_df.shape[0]
fig_size = max(8, n * 0.7)

plt.figure(figsize=(fig_size, fig_size))

ax = sns.heatmap(
    dist_df,
    cmap="vlag_r",                 # dark = more similar
    annot=True,
    fmt=".2f",
    annot_kws={"size": 12},
    square=True,
    linewidths=0.5,
    cbar_kws={"shrink": 0.5},
    vmin=0.05,
    vmax=0.30                          # Bray–Curtis is bounded [0,1]
)

#increase bar font
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=16)   # <-- change size here

plt.title("Bray–Curtis distance", fontsize=25)
plt.xticks(rotation=45, ha="right", fontsize=18)
plt.yticks(rotation=0, fontsize=18)

plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()

print("Annotated blue heatmap saved to:", output_file)