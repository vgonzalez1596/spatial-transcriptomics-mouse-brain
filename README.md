# Project #2: Spatial Transcriptomics Analysis of the Mouse Brain.

This repository contains a spatial transcriptomics analysis of 10x Genomics Visium HD mouse brain data using Seurat and BANKSY to identify spatially restricted gene expression patterns and anatomical domains.

---

## Project Goal

Determine what insights can be gained regarding gene expression in distinct regions of the mammalian brain. Comparison of brain regions can yield information regarding their unique functions. 

---

## Approach

- Carry out proper importing, quality control, and processing of this large Visium HD dataset. As this is a Visium HD spatial dataset, we will perform unsupervised clustering using Seurat’s sketch-based workflow.
- Segment the brain into biologically relevant spatial tissue domains with BANKSY.
- Utilize differential expression analysis between hippocampus and white matter to identify region-specific marker genes.
- Conduct KEGG pathway enrichment analysis on the hippocampus and white matter to identify functional differences between the two tissues. 

---

## Results

- Spatially restricted clusters correspond to known anatomical regions.
- White matter shows strong enrichment of myelination-associated genes (e.g., Plp1).
- Hippocampus clusters are enriched for neuronal differentiation and synaptic signaling pathways.
- KEGG enrichment highlights region-specific biological processes, with 

---

## Example Outputs
KEGG enrichment highlights region-specific biological processes, with 
[Variance Decomposition Analysis](results/figures/03_plot_variance_explained_by_factor.png)

---

## Tools & Packages

R, Seurat, Banksy, msigdbr, ggplot2. 

---

## Repository Structure

`spatial_transcriptomics_Visium_mouse_brain.R` analysis script

`results/figures/` contains exported analysis figures

`environment/` installs required packages and contains session information

---

## Data availability

The data used in this project is from the 10x Genomics Visium HD Mouse Brain (8µm bin size) dataset. This dataset must be [downloaded separately here](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he).
