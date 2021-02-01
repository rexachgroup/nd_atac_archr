# Root directory for fastqs. Each library should be in a separate, named directory.
FASTQ_DIR <- normalizePath("../fastqs/PreCG_final_ATAC2020/DG_32_10X_ATAC_S-20-1357_GAP233/")

CELLRANGER_ATAC_COUNT_DIR <- normalizePath("../data/cellranger-atac-count/precg-atac-2020")
CELLRANGER_ATAC_COUNT_BATCHTOOLS <- paste0(CELLRANGER_ATAC_COUNT_DIR, "-batchtools")

CELLRANGER_ATAC_AGGR_ID <- "precg-atac-2020"
CELLRANGER_ATAC_AGGR_DIR <- normalizePath(file.path("../data/cellranger-atac-aggr"))
CELLRANGER_ATAC_AGGR_BATCHTOOLS <- paste0(file.path(CELLRANGER_ATAC_AGGR_DIR, CELLRANGER_ATAC_AGGR_ID), "-batchtools")

CELLRANGER_ATAC_BIN <- "/geschwindlabshares/lchenprj01/software/seq/cellranger-atac-1.2.0/cellranger-atac"
CELLRANGER_ATAC_REFERENCE <- "/geschwindlabshares/lchenprj01/software/seqdata/cellranger-atac/refdata-cellranger-atac-GRCh38-1.2.0"

GIGA_TO_GIBI <- (1000 ** 3) / (1024 ** 3)
