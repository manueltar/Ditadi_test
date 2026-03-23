# Ditadi_test
scRNA-seq Analysis Pipeline: From Raw Fastqs to Cell Identity
This repository contains a modular, end-to-end bioinformatics pipeline for processing single-cell RNA-sequencing (scRNA-seq) data. The workflow is designed for High-Performance Computing (HPC) environments using the SLURM scheduler and integrates Seurat (R) with CellTypist (Python).

# Datasets

FBMA: https://www.nature.com/articles/s41586-021-03929-x Blood and immune development in human fetal bone marrow and Down syndrome

# Pipeline

--------------> Path backed up: /group/soranzo/manuel.tardaguila/2026_Ditadi_test/

-------------- Path in scratch: /scratch/manuel.tardaguila/2026_Ditadi_test



# 1. cellranger mapping and counts


$ sbatch ~/Scripts/sbatch/13_Cell_ranger_count_for_GEX_libraries.sh /scratch/manuel.tardaguila/2026_Ditadi_test MCO_01381_3GEX /group/soranzo/manuel.tardagu\
ila/2026_Ditadi_test/fastq_raw/

$ sbatch ~/Scripts/sbatch/13_Cell_ranger_count_for_GEX_libraries.sh /scratch/manuel.tardaguila/2026_Ditadi_test MCO_01382_3GEX /group/soranzo/manuel.tardagu\
ila/wt_day_30_MCO_01382/fastq_raw/



# 2. cellbender correction

$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2026_Ditadi_test/ MCO_01381_3GEX

$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2026_Ditadi_test/ MCO_01382_3GEX



# 3. Seurat First pass

$ bash ~/Scripts/Wraper_scripts/188_Seurat_First_pass_only_scRNAseq_v_Ditadi_test.sh /scratch/manuel.tardaguila/2026_Ditadi_test/ processing_outputs


# 4. Seurat Second pass

$ bash ~/Scripts/Wraper_scripts/189_Seurat_Second_pass_only_scRNAseq_v_Ditadi_test.sh /scratch/manuel.tardaguila/2026_Ditadi_test/

# 5. Merge samples in 1 object

$ bash ~/Scripts/Wraper_scripts/190_Seurat_merge_samples_only_scRNAseq_v_Ditadi_test.sh /scratch/manuel.tardaguila/2026_Ditadi_test/ processing_outputs

# 6. QC

----> Jupyter notebook: Final_QC_in_the_merged_object.ipynb

# 7. Recluster and export h5ad for Cell typist

$ bash ~/Scripts/Wraper_scripts/191_Recluster_at_low_res_and_export_h5ad.sh /scratch/manuel.tardaguila/2026_Ditadi_test/ processing_outputs


# 8. Cell typist prediction

----> Jupyter notebook: Cell_Typist_triple_prediction_cell_identity.ipynb
----> Jupyter notebook: explore_PRE_RPCA_integration.ipynb

# 9. Rpca integration

$ bash ~/Scripts/Wraper_scripts/192_RPCA_and_clustering_v2.sh /scratch/manuel.tardaguila/2026_Ditadi_test/ processing_outputs

# 10. Explore post RPCA and perform manual annotation

----> Jupyter notebook: explore_POST_RPCA_integration.ipynb

----> Jupyter notebook: mapping_cell_types.ipynb

----> Jupyter notebook: Subclustering.ipynb



