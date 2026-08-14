#!/usr/bin/env python3

# Part III: Bray–Curtis distance analysis
# Calculates pairwise Bray–Curtis distances, performs
# hierarchical clustering, generates a heatmap, and saves
# the distance matrix for downstream analysis in R.
#
# Code by Aylla von Ermland
###############################################

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list

input_file = "cluster_relative_frequencies_per_sample.csv"
output_file = "heatmap_braycurtis_blue_annotated.png"
distance_output_file = "braycurtis_distance_matrix.csv"

#Load relative k-mer cluster frequencies
data = pd.read_csv(
    input_file,
    index_col=0
)

data = data.apply(pd.to_numeric)
#Calculate pairwise Bray-Curtis dissimilarities between samples.
bray_dist = pdist(
    data.values,
    metric="braycurtis"
)
#Convert distances into a square matrix.
dist_array = squareform(
    bray_dist
)
dist_df = pd.DataFrame(
    dist_array,
    index=data.index,
    columns=data.index
)
#Remove row/column names such as "sample_ids".
dist_df.index.name = None
dist_df.columns.name = None
#Save Bray-Curtis distance matrix
dist_df.to_csv(
    distance_output_file
)
print(
    "Bray-Curtis distance matrix saved to:",
    distance_output_file
)
#Hierarchical clustering
#Perform average-linkage hierarchical clustering using the Bray-Curtis dissimilarities.
linkage_matrix = linkage(
    bray_dist,
    method="average"
)
#Obtain sample order from the hierarchical clustering.
ordered_indices = leaves_list(
    linkage_matrix
)
ordered_labels = dist_df.index[
    ordered_indices
]
#Reorder rows and columns of the distance matrix according to the hierarchical clustering.
dist_df = dist_df.loc[
    ordered_labels,
    ordered_labels
]
#Plot Bray-Curtis heatmap
n = dist_df.shape[0]

fig_size = max(
    8,
    n * 0.7
)

plt.figure(
    figsize=(
        fig_size,
        fig_size
    )
)

ax = sns.heatmap(
    dist_df,
    cmap="vlag_r",
    annot=True,
    fmt=".2f",
    annot_kws={
        "size": 12
    },
    square=True,
    linewidths=0.5,
    cbar_kws={
        "shrink": 0.5
    },
    vmin=0,
    vmax=0.4
)
#Figure formatting
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(
    labelsize=16
)

ax.set_xlabel("")
ax.set_ylabel("")

plt.title(
    "Bray-Curtis distance",
    fontsize=25
)

plt.xticks(
    rotation=45,
    ha="right",
    fontsize=18
)

plt.yticks(
    rotation=0,
    fontsize=18
)

plt.tight_layout()

#Save heatmap
plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    "Bray-Curtis heatmap saved to:",
    output_file
)