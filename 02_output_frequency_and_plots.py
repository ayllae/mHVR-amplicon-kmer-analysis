#!/usr/bin/env python3

#Part 2: Visualize and summarize clustering results. Steps:
#This code loads the clustered reads generated in Part 1, computes UMAP coordinates,
#creates PCA and UMAP visualizations, and outputs tables showing the number
#and relative frequency of reads assigned to each cluster per sample.

#Code by Aylla Ermland (ayllaermland@gmail.com)
######################

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

# Load clustering results stored in "adata_clusters.h5ad" file from Part 1.
print("Loading clustering results...")
adata=sc.read_h5ad("adata_clusters.h5ad")

#Compute UMAP coordinates from the neighbor graph. UMAP is only used for visualization and does not affect previous clustering. 
#Clustering was performed in previous 01_cluster_reads.py by Leiden clustering performed on the PCA-reduced data.
print("Computing UMAP...")
sc.tl.umap(adata)

#Count how many clusters exist.
n_clusters=len(
    adata.obs["cluster"]
    .cat
    .categories
)

#high resolution figures.
mpl.rcParams["figure.dpi"]=300
mpl.rcParams["savefig.dpi"]=300

#Remove plot frame and use white background.
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False
)

#Assign a unique color to each cluster so PCA and UMAP use a consistent color scheme.
cmap=plt.get_cmap("tab20")
colors=[
    mpl.colors.to_hex(cmap(i%20))
    for i in range(n_clusters)
]

#Store colors so Scanpy automatically uses them.
adata.uns["cluster_colors"]=colors

#PCA plot showing relationships among reads in principal component space.
print("Generating PCA plot...")
sc.pl.pca(
    adata,
    color="cluster",
    palette=colors,
    size=25,
    save="_pca_clusters.pdf",
    show=False
)

#Global UMAP showing all reads colored by cluster.
print("Generating UMAP plot...")
sc.pl.umap(
    adata,
    color="cluster",
    palette=colors,
    title="UMAP of k-mer frequencies (Leiden clustering)",
    size=25,
    save="_umap_clusters.pdf",
    show=False
)

#Generate one UMAP plot per sample.
#Coordinates remain identical to the global UMAP; Only reads from a single sample are displayed.
print("Generating sample-specific UMAPs...")

for sample in adata.obs["sample"].unique():
    # Select reads belonging to one sample.
    subset=adata[
        adata.obs["sample"]==sample,
        :
    ]
    # Plot only that sample's reads.
    sc.pl.umap(
        subset,
        color="cluster",
        palette=colors,
        title=f"UMAP - Sample: {sample}",
        size=25,
        save=f"_umap_clusters_{sample}.pdf",
        show=False
    )
# Count how many reads from each sample belong
# to each cluster.
print("Generating cluster tables...")
cluster_table=(
    adata.obs.groupby(
        ["sample","cluster"],
        observed=False
    )
    .size()
    .unstack(fill_value=0)
)
#Save raw read counts.
cluster_table.to_csv(
    "cluster_counts_per_sample.csv"
)
#Calculate total reads present in each sample.
sample_read_counts=(
    adata.obs["sample"]
    .value_counts()
    .sort_index()
)
#Convert raw counts into relative frequencies, based on total number of reads in each sample. 
relative_freq_table=cluster_table.div(
    sample_read_counts,
    axis=0
)

#Save relative frequencies.
relative_freq_table.to_csv(
    "cluster_relative_frequencies_per_sample.csv"
)

print("\nFinished Part 2.")