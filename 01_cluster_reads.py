#!/usr/bin/env python3

# Cluster minicircle hypervariable region amplicon read sequences based on their kmer frequency composition. 

#Steps: 
# I. Convert each read into a vector containing the frequency of every possible DNA k-mer;
# II. Normalize the k-mer (k=7) counts and reduce dimensionality using PCA.
#	Each read is scaled so that its total k-mer count equals 10,000. 
#   This reduces the effect of differences in read length and allows reads to be compared 
#   mainly by their relative k-mer composition. 
# III. Build a nearest-neighbor graph in PCA space and identify groups of similar reads
#	using Leiden clustering. 
# IV. Remove low-frequency clusters containing less than 1% of all reads;
# V. Save the final AnnData object for downstream analysis. 

#Code by Aylla von Ermland (ayllaermland@gmail.com)
###############################################

#Import modules
import os,re
import numpy as np
from multiprocessing import Pool
from Bio import SeqIO
import scanpy as sc
import anndata as ad
from scipy import sparse

# Input FASTA containing all reads from all hypervariable amplicon samples pulled together.
# each read header should be formatted as >sampleID_month_year_Read#
fasta_file="all_sequences.fasta"

# K-mer size. Each read will be represented by the frequencies of all possible 7-mer.
k=7

# Total possible DNA 7-mers. Because each position can contain A, C, G or T:
# 4^7 = 16,384 possible 7-mer. 
n_kmers=4**k

#Converting bases to integers for faster performance.
# A=0, C=1, G=2, T=3; all other byte values are initially assigned 255.
base_map=np.full(256,255,dtype=np.uint8)
base_map[ord('A')]=0
base_map[ord('C')]=1
base_map[ord('G')]=2
base_map[ord('T')]=3

#Function that converts one DNA sequence into a vector of 7-mer counts. 
#The returned vector has 16,384 positions, with one position for each possible DNA 7-mer.
#Each value records how many times that 7-mer occurs in the sequence.
def kmer_freq_fast(seq):
    seq=seq.encode("ascii") #Convert sequence string into bytes for fast processing.
    
    seq_array=base_map[np.frombuffer(seq,dtype=np.uint8)] #Convert bases into integer representation.
    
    seq_array=seq_array[seq_array<4] #Remove any non-ACGT characters.

    #Create vector that stores counts for every possible 7-mer.
    counts=np.zeros(n_kmers,dtype=np.float32)
    
    #Reads shorter than 7 bp cannot contain a 7-mer.
    if len(seq_array)<k:
        return counts
    #Loop through all reads, and build first k-mer using rolling 2-bit encoding.
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
    return counts #return the complete 7-mer count profile for this read.

def extract_read_metadata(read_id):
    # Expected format:
    # Sample1_January_2025_Read1234
    parts=read_id.rsplit("_",3)

    if len(parts)!=4:
        raise ValueError(
            f"Unexpected read ID format: {read_id}"
        )

    sample_id,month,year,read_number=parts
    sample_date=f"{month}_{year}"

    return sample_id,sample_date

def process_record(record):
    # Process one read and return:
    # k-mer profile, complete read ID, sample ID and sample date.
    read_id=record.id
    sample_id,sample_date=extract_read_metadata(read_id)

    return(
        kmer_freq_fast(str(record.seq)),
        read_id,
        sample_id,
        sample_date
    )
    
if __name__=="__main__":
    print("Parsing FASTA...") #Load all sequences into memory.
    records=list(SeqIO.parse(fasta_file,"fasta"))
    n_reads=len(records) #Number of reads.
    print(f"Total reads: {n_reads}")
    #Create matrix where each rows is a read and each column 16384 possible 7-mers
    X=np.zeros((n_reads,n_kmers),dtype=np.float32)
    #Store metadata for each read.
    seq_ids=[]
    sample_ids=[]
    sample_dates=[]
    
    #Use SLURM core count if running code on server, otherwise use all available CPU cores.
    n_cores=int(
        os.environ.get(
            "SLURM_CPUS_PER_TASK",
            os.cpu_count()
        )
    )
    print(f"Using {n_cores} cores")
    #Count k-mers in parallel. Each worker processes a batch of reads and returns: counts, read ID and sample ID.
    with Pool(n_cores) as pool:
        for i,(counts,rid,sid,sdate) in enumerate(
            pool.imap(
                process_record,
                records,
                chunksize=500
            )
        ):
            X[i]=counts  #Store k-mer profile.
            seq_ids.append(rid)  #Store metadata.
            sample_ids.append(sid)
            sample_dates.append(sdate)

    #Convert dense matrix into sparse format. Most possible k-mers are absent from most reads, so sparse matrices save a large amount of RAM.
    X_sparse=sparse.csr_matrix(X)
    #Remove dense matrix from memory.
    del X

    # Create the Scanpy AnnData object.
    # Each matrix row represents one read.
    adata=ad.AnnData(X_sparse)

    # Use the complete FASTA read IDs as the AnnData row names.
    adata.obs_names=seq_ids
    adata.obs_names.name="read_id"

    # Store the sample ID associated with each read.
    adata.obs["sample"]=sample_ids
    adata.obs["date"]=sample_dates
    
    #Normalize each read to the same total count. This reduces the effect of different read lengths.
    sc.pp.normalize_total(
        adata,
        target_sum=1e4 #rescale each read so the total k-mer count equals 10,000
    )

    #Log-transform values to reduce dominance of extremely abundant k-mers.
    sc.pp.log1p(adata)

    #PCA reduces 16,384 dimensions into 30 principal components. These components capture the major variation among reads.
    sc.pp.pca(
        adata,
        n_comps=30,
        svd_solver="randomized"
    )

    #Build nearest-neighbor graph in PCA space. Reads with similar k-mer composition become connected.
    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=15
    )
    #Leiden community detection identifies groups of reads with similar k-mer composition.
    sc.tl.leiden(
        adata,
        resolution=0.5
    )
    #Rename clusters from: 0,1,2,3 ... to 1,2,3,4...
    adata.obs["cluster"]=(
        adata.obs["leiden"]
        .astype(int)
        .add(1)
        .astype(str)
        .astype("category")
    )

    #Preserve numerical cluster ordering.
    adata.obs["cluster"]=(
        adata.obs["cluster"]
        .cat
        .as_ordered()
    )

    #Print number of reads assigned to each cluster.
    print("\nCluster sizes:")
    print(
        adata.obs["cluster"]
        .value_counts()
    )

    #Save complete clustering object.
    #This file contains k-mer matrix; PCA coordinates; neighbor graph; Leiden clusters; sample metadata and is loaded by Code 2 for downstream analysis.
    adata.write("adata_clusters.h5ad")
    print("\nFinished clustering.")
    
    print("\nCluster sizes:")
    print(
        adata.obs["cluster"]
        .value_counts()
    )

    #Filter out clusters with less than 1% of total reads.
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