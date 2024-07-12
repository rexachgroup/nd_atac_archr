# Introduction
This repo contains code for the scATAC analysis of Rexach et al. (2024). For the scRNA analysis see [rexachgroup/snseq_nd_liger](https://github.com/rexachgroup/snseq_nd_liger).

# Installation / Running
Dependencies for this code are specified using a conda environment; create and activate the environment by running
```
conda create -p conda/ -f environment.yml
conda activate conda/
```

Some steps require R packages that are not available from `conda-forge`  or `bioconda`. Install the dependencies by running `install_archr.R` after the conda environment is created.

# analysis/ folder organization
Subcomponents of the analysis are organized in run order under `src/`:

## `00_import_transfac.R` `00_create_scenic_transfac_db.R`
Import regulons from SCENIC analysis and convert TRANSFAC PRO database.

## `01_cellranger_count`, `02b_cellranger_aggr.R`
Cellranger-ATAC 1.2.0 run.

## `03_archr_import.R` `04_archr_doublet.R` `05_archr_qc.R` `06_archr_harmony_clustering.R`
Import fragment / bed files to ArchR, including ArchR doublet filtering.
ArchR clustering of peak matrix.

## `07_archr_astrocyte_cluster_composition.R` `07_archr_astrocyte_cluster_plotting.R`
Plotting of astrocyte subgroup (C18-C25).

## `08_archr_liger_integration_insula.R` `08_archr_liger_integration_precg.R`
Per-region annotation of gene activity score matrix with snRNAseq data.

## `09_archr_cluster_peak_calls.R`
Peak re-calling using MACS2.

## `10_archr_transcription_factor_transfac_c2.R` `10_archr_transcription_factor_transfac_c7.R`
Correlation of excitatory and microglia cluster peaks with TFs.

## `15_heterochromatin/`
Heterochromatin analysis using ALLCools.
