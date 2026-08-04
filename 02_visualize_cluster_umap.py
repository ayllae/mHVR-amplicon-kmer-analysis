#!/usr/bin/env python3

#Part 2: Visualize and summarize clustering results. Steps:

#I. This code loads the clustered reads generated in Part 1, computes UMAP coordinates,
#II. Creates PCA and UMAP visualizations, and outputs tables showing the number
# and relative frequency of reads assigned to each cluster per sample.

#Code by Aylla von Ermland (ayllaermland@gmail.com)
######################

# Import modules
import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl


# Load clustering results stored in the AnnData file created in Part 1.
print("Loading clustering results...")
adata=sc.read_h5ad("adata_clusters.h5ad")

# Confirm that sample and date metadata are present.
required_columns={"sample","date"}
missing_columns=required_columns-set(adata.obs.columns)

if missing_columns:
    raise ValueError(
        f"Missing metadata columns in AnnData: {sorted(missing_columns)}"
    )

# Combine sample ID and collection date into one label.
# Example:
# sample = Sample1
# date = January_2025
# sample_label = Sample1_January_2025

adata.obs["sample_label"]=(
    adata.obs["sample"].astype(str)
    +"_"
    +adata.obs["date"].astype(str)
)

# Create directory where all figures will be saved.
figure_dir="Figures"
os.makedirs(figure_dir,exist_ok=True)

# Compute UMAP coordinates from the neighbor graph.
# UMAP is used only for visualization and does not change the Leiden
# clustering performed in Part 1.
print("Computing UMAP...")
sc.tl.umap(adata)

# Store global UMAP limits so all sample-specific UMAP plots use
# the same coordinate scale.
xmin=adata.obsm["X_umap"][:,0].min()
xmax=adata.obsm["X_umap"][:,0].max()
ymin=adata.obsm["X_umap"][:,1].min()
ymax=adata.obsm["X_umap"][:,1].max()

# Count how many clusters are present.
n_clusters=len(
    adata.obs["cluster"]
    .cat
    .categories
)

# Set high-resolution figure output.
mpl.rcParams["figure.dpi"]=300
mpl.rcParams["savefig.dpi"]=300

# Use a white background and visible plot frame.
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=True
)

# Assign a unique color to each cluster.
colors=sc.pl.palettes.godsnot_102[:n_clusters]

# Store cluster categories in their existing order.
cluster_categories=adata.obs["cluster"].cat.categories

# Create an explicit cluster-to-color mapping.
# This guarantees that each cluster keeps the same color in every plot,
# including sample-specific plots that contain only some clusters.

color_dict={
    cluster:colors[i]
    for i,cluster in enumerate(cluster_categories)
}

# Store colors in the AnnData object so Scanpy uses them consistently.
adata.uns["cluster_colors"]=colors

# Generate PCA plot showing all reads colored by cluster.
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

# Add axis labels.
ax.set_xlabel("PC1",fontsize=10)
ax.set_ylabel("PC2",fontsize=10)

# Add plot title.
ax.set_title("All samples",fontsize=10)

# Show border around the plot.
for spine in ax.spines.values():
    spine.set_visible(True)

# Add whitespace around points.
ax.margins(0.10)

# Set tick label size.
ax.tick_params(labelsize=8)

# Save PCA plot.
fig.savefig(
    os.path.join(
        figure_dir,
        "PCA_clusters.pdf"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

# Generate global UMAP showing all reads colored by cluster.
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

# Add axis labels.
ax.set_xlabel("UMAP1",fontsize=10)
ax.set_ylabel("UMAP2",fontsize=10)

# Add plot title.
ax.set_title("All samples",fontsize=10)

# Show border around the plot.
for spine in ax.spines.values():
    spine.set_visible(True)

# Add whitespace around points.
ax.margins(0.10)

# Set tick label size.
ax.tick_params(labelsize=8)

# Keep the global coordinate scale.
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

# Save global UMAP plot.
fig.savefig(
    os.path.join(
        figure_dir,
        "UMAP_clusters.pdf"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

# Generate one UMAP plot for each sample and collection date.
# Coordinates remain identical to the global UMAP. Only reads belonging
# to one sample and date are displayed in each plot.
print("Generating sample-specific UMAPs...")

for sample_label in sorted(adata.obs["sample_label"].unique()):

    # Select reads from one sample and collection date.
    subset=adata[
        adata.obs["sample_label"]==sample_label,
        :
    ].copy()

    # Keep all original cluster categories.
    #
    # This prevents Scanpy from reordering or dropping cluster IDs.
    subset.obs["cluster"]=(
        subset.obs["cluster"]
        .cat
        .set_categories(cluster_categories)
    )

    # Keep the same cluster colors used in the full dataset.
    subset.uns["cluster_colors"]=adata.uns["cluster_colors"]

    # Generate the sample-specific UMAP.
    fig=sc.pl.umap(
        subset,
        color="cluster",
        palette=color_dict,
        size=15,
        show=False,
        return_fig=True
    )

    ax=fig.axes[0]

    # Use sample ID and date as the plot title.
    ax.set_title(
        sample_label,
        fontsize=10
    )

    # Add axis labels.
    ax.set_xlabel(
        "UMAP1",
        fontsize=10
    )

    ax.set_ylabel(
        "UMAP2",
        fontsize=10
    )

    # Show border around the plot.
    for spine in ax.spines.values():
        spine.set_visible(True)

    # Add whitespace around points.
    ax.margins(0.10)

    # Set tick label size.
    ax.tick_params(
        labelsize=8
    )

    # Keep the same coordinate scale for every sample-specific plot.
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # Save the sample-specific UMAP.
    fig.savefig(
        os.path.join(
            figure_dir,
            f"UMAP_{sample_label}.pdf"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


# Count how many reads from each sample and date belong to each cluster.
print("Generating cluster tables...")

cluster_table=(
    adata.obs.groupby(
        ["sample_label","cluster"],
        observed=False
    )
    .size()
    .unstack(fill_value=0)
)


# Save raw read counts for each sample-date combination.
cluster_table.to_csv(
    "cluster_counts_per_sample.csv"
)


# Calculate the total number of reads in each sample-date combination.
sample_read_counts=(
    adata.obs["sample_label"]
    .value_counts()
    .sort_index()
)


# Convert raw cluster counts into relative frequencies.
# Each cluster count is divided by the total number of reads from the
# corresponding sample and collection date.
relative_freq_table=cluster_table.div(
    sample_read_counts,
    axis=0
)

# Save relative cluster frequencies.
relative_freq_table.to_csv(
    "cluster_relative_frequencies_per_sample.csv"
)

print("\nFinished Part 2.")
