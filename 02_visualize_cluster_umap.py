#!/usr/bin/env python3

#Part 2: Visualize and summarize clustering results. Steps:

#I. This code loads the clustered reads generated in Part 1, computes UMAP coordinates,
#II. Creates PCA and UMAP visualizations, and outputs tables showing the number
# and relative frequency of reads assigned to each cluster per sample.

#Code by Aylla von Ermland (ayllaermland@gmail.com)
######################

#Import modules
import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

#Load clustering results stored in "adata_clusters.h5ad" file from Part 1.
print("Loading clustering results...")
adata=sc.read_h5ad("adata_clusters.h5ad")

#Create directory where all figures will be saved
figure_dir="Figures"
os.makedirs(figure_dir,exist_ok=True)

#Compute UMAP coordinates from the neighbor graph. UMAP is only used for visualization and does not affect previous clustering.
#Clustering was performed in previous 01_cluster_reads.py by Leiden clustering performed on the PCA-reduced data.
print("Computing UMAP...")
sc.tl.umap(adata)

#Store global UMAP limits so all sample UMAPs use the same scale.
xmin,xmax=adata.obsm["X_umap"][:,0].min(),adata.obsm["X_umap"][:,0].max()
ymin,ymax=adata.obsm["X_umap"][:,1].min(),adata.obsm["X_umap"][:,1].max()

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
    frameon=True
)

#Assign a unique color to each cluster so PCA and UMAP use a consistent color scheme.
colors=sc.pl.palettes.godsnot_102[:n_clusters]

#Create explicit cluster -> color mapping.
#This guarantees that cluster IDs always have the same color,
#even when plotting individual samples that contain only a subset of clusters.
cluster_categories=adata.obs["cluster"].cat.categories

color_dict={
    cluster:colors[i]
    for i,cluster in enumerate(cluster_categories)
}

#Store colors so Scanpy automatically uses them.
adata.uns["cluster_colors"]=colors

#PCA plot showing relationships among reads in principal component space.
print("Generating PCA plot...")

fig=sc.pl.pca(
    adata,
    color="cluster",
    palette=color_dict,
    size=15,
    show=False,
    return_fig=True
)

ax=fig.axes[0]

#axis labels
ax.set_xlabel("PC1",fontsize=10)
ax.set_ylabel("PC2",fontsize=10)

#title font
ax.set_title("All samples",fontsize=10)

#show border
for spine in ax.spines.values():
    spine.set_visible(True)

#add margin around points
ax.margins(0.1)

#tick labels
ax.tick_params(labelsize=8)

fig.savefig(
    os.path.join(
        figure_dir,
        "PCA_clusters.pdf"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

#Global UMAP showing all reads colored by cluster.
print("Generating UMAP plot...")

fig=sc.pl.umap(
    adata,
    color="cluster",
    palette=color_dict,
    size=15,
    show=False,
    return_fig=True
)

ax=fig.axes[0]

#axis labels
ax.set_xlabel("UMAP1",fontsize=10)
ax.set_ylabel("UMAP2",fontsize=10)

#x-axis title
ax.set_title("All samples",fontsize=10)

#show border
for spine in ax.spines.values():
    spine.set_visible(True)

#add whitespace around points
ax.margins(0.10)

#tick labels
ax.tick_params(labelsize=8)

#keep same coordinate scale
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

fig.savefig(
    os.path.join(
        figure_dir,
        "UMAP_clusters.pdf"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

#Generate one UMAP plot per sample.
#Coordinates remain identical to the global UMAP; Only reads from a single sample are displayed.
print("Generating sample-specific UMAPs...")

for sample in adata.obs["sample"].unique():

    # select reads from one sample
    subset=adata[
        adata.obs["sample"]==sample,
        :
    ].copy()

    #Keep all original cluster categories.
    #This prevents Scanpy from reordering or dropping cluster IDs.
    subset.obs["cluster"]=(
        subset.obs["cluster"]
        .cat
        .set_categories(cluster_categories)
    )

    #Keep the same cluster colors used in the full dataset.
    subset.uns["cluster_colors"]=adata.uns["cluster_colors"]

    # generate figure object instead of letting scanpy save automatically
    fig=sc.pl.umap(
        subset,
        color="cluster",
        palette=color_dict,
        size=15,
        show=False,
        return_fig=True
    )

    ax=fig.axes[0]

    ax.set_title(
        sample,
        fontsize=10
    )

    # add axis labels
    ax.set_xlabel(
        "UMAP1",
        fontsize=10
    )

    ax.set_ylabel(
        "UMAP2",
        fontsize=10
    )

    # show border around plot
    for spine in ax.spines.values():
        spine.set_visible(True)

    # add whitespace around points
    ax.margins(0.10)

    # smaller tick labels
    ax.tick_params(
        labelsize=8
    )

    # keep same coordinate scale
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    fig.savefig(
        os.path.join(
            figure_dir,
            f"UMAP_{sample}.pdf"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

# Count how many reads from each sample belong to each cluster.
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