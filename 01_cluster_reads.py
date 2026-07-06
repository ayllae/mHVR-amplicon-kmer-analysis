#!/usr/bin/env python3

# Part 1: This code Cluster reads based on 7-mer frequency composition. Steps: 
# I. Reads are converted into 7-mer frequency profiles,
# II. clustered using Leiden clustering in PCA space,
# III. and saved as adata_clusters.h5ad for downstream analysis.

#Code by Aylla Ermland (ayllaermland@gmail.com)
###############################################

#Import modules
import os,re
import numpy as np
from multiprocessing import Pool
from Bio import SeqIO
import scanpy as sc
import anndata as ad
from scipy import sparse

# Input FASTA containing all reads from all samples pulled together.
fasta_file="all_sequences.fasta"

# K-mer size. Each read will be represented by its 7-mer composition.
k=7

# Total possible DNA 7-mers = 4^7 = 16384.
n_kmers=4**k

#Converting bases to integers for faster performance
# A=0, C=1, G=2, T=3.
base_map=np.full(256,255,dtype=np.uint8)
base_map[ord('A')]=0
base_map[ord('C')]=1
base_map[ord('G')]=2
base_map[ord('T')]=3

def kmer_freq_fast(seq):

    # Convert sequence string into bytes for fast processing.
    seq=seq.encode("ascii")

    # Convert bases into integer representation.
    seq_array=base_map[np.frombuffer(seq,dtype=np.uint8)]

    # Remove any non-ACGT characters.
    seq_array=seq_array[seq_array<4]

    # Create vector that stores counts for every possible 7-mer.
    counts=np.zeros(n_kmers,dtype=np.float32)

    # Reads shorter than 7 bp cannot contain a 7-mer.
    if len(seq_array)<k:
        return counts

    # Build first k-mer using rolling 2-bit encoding.
    kmer=0
    for i in range(k):
        kmer=((kmer<<2)|int(seq_array[i]))

    # Count first k-mer.
    counts[kmer]+=1

    # Bit mask keeps only current k-mer information.
    mask=(1<<(2*k))-1

    # Slide across sequence one base at a time and count all k-mers.
    for i in range(k,len(seq_array)):
        kmer=(((kmer<<2)|int(seq_array[i]))&mask)
        counts[kmer]+=1
    return counts

def extract_sample_id(header):
    # Extract sample name from headers formatted as example: "Sample1_Read12345", and returns "Sample1".
    match=re.match(r'^(\S+?)_Read',header)
    return match.group(1) if match else "unknown"

def process_record(record):
    # Process one read and return: k-mer count vector, read ID and sample ID
    return(
        kmer_freq_fast(str(record.seq)),
        record.id,
        extract_sample_id(record.id)
    )
if __name__=="__main__":
    # Load all sequences into memory.
    print("Parsing FASTA...")
    records=list(SeqIO.parse(fasta_file,"fasta"))
    # Number of reads.
    n_reads=len(records)
    print(f"Total reads: {n_reads}")
    # Create matrix where each rows is a read and each column 16384 possible 7-mers
    X=np.zeros((n_reads,n_kmers),dtype=np.float32)
    # Store metadata for each read.
    seq_ids=[]
    sample_ids=[]

    # Use SLURM core count if running on cluster, otherwise use all available CPU cores.
    n_cores=int(
        os.environ.get(
            "SLURM_CPUS_PER_TASK",
            os.cpu_count()
        )
    )
    print(f"Using {n_cores} cores")
    # Count k-mers in parallel. Each worker processes a batch of reads and returns: counts, read ID and sample ID.
    with Pool(n_cores) as pool:
        for i,(counts,rid,sid) in enumerate(
            pool.imap(
                process_record,
                records,
                chunksize=500
            )
        ):
            # Store k-mer profile.
            X[i]=counts
            # Store metadata.
            seq_ids.append(rid)
            sample_ids.append(sid)

    # Convert dense matrix into sparse format. Most possible k-mers are absent from most reads, so sparse matrices save a large amount of RAM.
    X_sparse=sparse.csr_matrix(X)
    # Remove dense matrix from memory.
    del X

    # Create Scanpy AnnData object.
    # This becomes the main data structure for clustering.
    adata=ad.AnnData(X_sparse)
    # Attach metadata to each read.
    adata.obs["read_id"]=seq_ids
    adata.obs["sample"]=sample_ids
    # Normalize each read to the same total count. This reduces the effect of different read lengths.
    sc.pp.normalize_total(
        adata,
        target_sum=1e4
    )

    # Log-transform values to reduce dominance of extremely abundant k-mers.
    sc.pp.log1p(adata)

    # PCA reduces 16384 dimensions into 30 principal components. These components capture the major variation among reads.
    sc.pp.pca(
        adata,
        n_comps=30,
        svd_solver="randomized"
    )

    # Build nearest-neighbor graph in PCA space. Reads with similar k-mer composition become connected.
    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=15
    )
    # Leiden community detection identifies groups of reads with similar k-mer composition.
    sc.tl.leiden(
        adata,
        resolution=0.5
    )
    # Rename clusters from: 0,1,2,3 ... to 1,2,3,4...
    adata.obs["cluster"]=(
        adata.obs["leiden"]
        .astype(int)
        .add(1)
        .astype(str)
        .astype("category")
    )

    # Preserve numerical cluster ordering.
    adata.obs["cluster"]=(
        adata.obs["cluster"]
        .cat
        .as_ordered()
    )

    # Print number of reads assigned to each cluster.
    print("\nCluster sizes:")
    print(
        adata.obs["cluster"]
        .value_counts()
    )

    # Save complete clustering object.
    # This file contains k-mer matrix; PCA coordinates; neighbor graph; Leiden clusters; sample metadata and is loaded by Code 2 for downstream analysis.
    adata.write("adata_clusters.h5ad")
    print("\nFinished clustering.")
    
    print("\nCluster sizes:")
    print(
        adata.obs["cluster"]
        .value_counts()
    )

    # Filter out clusters with less than 1% of total reads.
    cluster_counts = adata.obs["cluster"].value_counts()
    min_reads = 0.01 * adata.n_obs

    clusters_to_keep = cluster_counts[cluster_counts >= min_reads].index

    print("\nKeeping clusters with >=1% of reads:")
    print(cluster_counts.loc[clusters_to_keep])

    print("\nRemoving clusters with <1% of reads:")
    print(cluster_counts[cluster_counts < min_reads])

    adata = adata[adata.obs["cluster"].isin(clusters_to_keep)].copy()
    adata.obs["cluster"] = adata.obs["cluster"].cat.remove_unused_categories()

    print(f"\nReads remaining after filtering: {adata.n_obs}")

    adata.write("adata_clusters.h5ad")
    print("\nFinished clustering.")