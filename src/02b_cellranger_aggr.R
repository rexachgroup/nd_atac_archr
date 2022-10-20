library(batchtools)
library(tidyverse)

FASTQ_DIR <- c(normalizePath("../fastqs/PreCG_final_ATAC2020/DG_32_10X_ATAC_S-20-1357_GAP233/"), normalizePath("../fastqs/PreCG_final_ATAC2020/DG_S-20-0434_1_pool_GAP169"))
FASTQ_META <- normalizePath("../data/snATAC_metadata_summary_2021_d.xlsx")

CELLRANGER_ATAC_AGGR_ID <- "precg-atac-2020-nonormalize"
CELLRANGER_ATAC_AGGR_DIR <- normalizePath(file.path("../data/cellranger-atac-aggr"))
CELLRANGER_ATAC_AGGR_BATCHTOOLS <- paste0(file.path(CELLRANGER_ATAC_AGGR_DIR, CELLRANGER_ATAC_AGGR_ID), "-batchtools")
CELLRANGER_ATAC_BIN <- "/geschwindlabshares/lchenprj01/software/seq/cellranger-atac-1.2.0/cellranger-atac"
CELLRANGER_ATAC_REFERENCE <- "/geschwindlabshares/lchenprj01/software/seqdata/cellranger-atac/refdata-cellranger-atac-GRCh38-1.2.0"

# count_reg <- loadRegistry(CELLRANGER_ATAC_COUNT_BATCHTOOLS)
fastqs <- list.dirs(FASTQ_DIR, recursive = FALSE)
aggregation_table <- tibble(fastqs) %>%
    mutate(
        library_id = basename(fastqs),
        fragments = file.path(CELLRANGER_ATAC_COUNT_DIR, library_id, "outs", "fragments.tsv.gz"),
        cells = file.path(CELLRANGER_ATAC_COUNT_DIR, library_id, "outs", "singlecell.csv"),
    ) %>%
    filter(file.exists(fragments), file.exists(cells))

if (dir.exists(CELLRANGER_ATAC_AGGR_BATCHTOOLS)) {
    aggr_reg <- loadRegistry(CELLRANGER_ATAC_AGGR_BATCHTOOLS, writeable = TRUE)
} else {
    aggr_reg <- makeRegistry(CELLRANGER_ATAC_AGGR_BATCHTOOLS)
}

aggregation_table_path <- file.path(CELLRANGER_ATAC_AGGR_BATCHTOOLS, "aggregation.csv")
write_csv(aggregation_table, aggregation_table_path)

RESOURCES <- list(
    ncpus = 64,
    memory = 320, # slurm/sge interprets as *giga*bytes, cellranger as *gibi*bytes
    walltime = 172800
)

cmd <- str_glue("cd {CELLRANGER_ATAC_AGGR_DIR} && {CELLRANGER_ATAC_BIN} aggr",
    "--id={CELLRANGER_ATAC_AGGR_ID}",
    "--csv={aggregation_table_path}",
    "--reference={CELLRANGER_ATAC_REFERENCE}",
    "--nosecondary",
    "--normalize=none",
    "--localcores={RESOURCES$ncpus}",
    "--localmem={floor(RESOURCES$memory * GIGA_TO_GIBI * 0.9)}",
    .sep = " ")

clearRegistry()
ids <- batchMap(system, list(cmd))
setJobNames(ids, c(CELLRANGER_ATAC_AGGR_ID))
submitJobs(ids, resources = RESOURCES)
waitForJobs()
