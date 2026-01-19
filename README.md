# Ditadi_test
scRNA-seq Analysis Pipeline: From Raw Fastqs to Cell Identity
This repository contains a modular, end-to-end bioinformatics pipeline for processing single-cell RNA-sequencing (scRNA-seq) data. The workflow is designed for High-Performance Computing (HPC) environments using the SLURM scheduler and integrates Seurat (R) with CellTypist (Python).

Pipeline Overview
The pipeline is organized into five distinct stages:

1. Raw Data Processing & Denoising
13_Cell_ranger_count_for_GEX_libraries.sh: Align reads and generate count matrices using cellranger count (v8.0.1) against the GRCh38 human reference.

7_CellBender_v_scratch.sh: Removes ambient RNA ("soup") and background noise from the raw Cell Ranger output using a deep generative model.

2. Quality Control & Individual Sample Processing
490_Seurat_first_pass_only_scRNAseq.R: Performs initial filtering. Default thresholds:

Minimum features: 500

Maximum mitochondrial content: 10%

491_Seurat_second_pass_only_scRNAseq.R: Secondary processing of filtered objects to prepare for integration.

3. Integration & Merging
492_Merge_samples_only_scRNAseq_v2.R: Merges individual samples into a unified Seurat object. This stage includes an initial doublet filtering step.

493_Clustering_of_merged_samples_only_scRNAseq.R: Performs standard dimensionality reduction (PCA, UMAP) and high-resolution clustering to identify detailed sub-populations.

4. Refinement & Cross-Language Export
494_merged_clustering_at_low_res.R: Re-clusters the merged dataset at a lower resolution (e.g., 0.5) to define stable cell lineages.

513_Export_RNA_modality...R: Extracts the RNA count matrices and metadata.

3_export_to_h5ad_v2.py: Converts the R-based matrices into the Python-standard AnnData (.h5ad) format, enabling compatibility with the Scanpy/CellTypist ecosystem.

5. Cell Type Annotation & Validation
Cell_Typist_triple_prediction_cell_identity.ipynb: Utilizes CellTypist to automatically assign cell identities using high-level immune models.

Final_QC_in_the_merged_object.ipynb: Performs post-clustering QC, including scDblFinder scoring to validate cluster-level doublet enrichment.

mapping_cell_types.ipynb: Maps the automated predictions from CellTypist back onto the Seurat UMAP for final visualization and manual validation.

# Datasets

Immune_All_Low: This is the default high-resolution model. It is trained on a comprehensive cross-tissue atlas comprising approximately 360,000 cells from 20 tissues across 19 studies. It can identify 90–98 detailed immune cell types and subtypes.
Alsinet: https://www.nature.com/articles/s41467-022-30557-4 Robust temporal map of human in vitro myelopoiesis using single-cell genomics
FBMA: https://www.nature.com/articles/s41586-021-03929-x Blood and immune development in human fetal bone marrow and Down syndrome

Requirements
R Packages
Seurat (v4 or v5)

Signac

scDblFinder

patchwork

optparse

Python Environment
Scanpy

CellTypist

anndata

rpy2 (for object conversion)