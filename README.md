# *Trypanosoma cruzi* Minicircle Hypervariable Region (mHVR) Sequence Analysis Pipeline

This repository contains the pipeline used to compare *Trypanosoma cruzi* kDNA mHVR amplicon samples based on k-mer composition. The main idea is to represent each read by its 7-mer composition, cluster reads with similar composition, and then compare samples based on the frequency of those clusters. 

Deep amplicon sequencing were performed using Oxford Nanopore sequencing of pooled PCR products targeting mHVR amplicons of multiple sizes. A target depth of approximately 6,000 reads per sample is recommended, although the pipeline also performs well with as few as 3,000 reads per sample. The pipeline has not yet been validated for other amplicon sequencing chemistries, such as Illumina.

Steps:
1. Filter reads by length (200–500 bp)
2. Convert each read into a 7-mer frequency vector
3. Normalize k-mer counts and apply `log(1+x)`
4. Perform PCA
5. Build a nearest-neighbor graph and run Leiden clustering
6. Calculate the frequency of each cluster in each sample
7. Calculate Bray-Curtis distances between samples

UMAP is used to visualize the read clusters.

```text
Raw sequencing files
        ↓
Process_and_clean_files.py
Read filtering
        ↓
01_cluster_reads.py
7-mer composition → normalization → PCA → Leiden clustering
        ↓
02_visualize_cluster_umap.py
UMAP visualization of read clusters
        ↓
03_bray_curtis_matrix.py
Cluster frequencies per sample → Bray-Curtis distance matrix
        ↓
04_Within_dog_analysis.r
Within-dog distance analysis
```

## References

Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*. 2019;9:5233.

Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15.

## Author
Aylla S. K. von Ermland
Tarleton Lab
University of Georgia
