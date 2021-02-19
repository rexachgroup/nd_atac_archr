source("00_CONSTANTS.R")
library(batchtools)
library(tidyverse)
library(readxl)
FASTQ_DIR <- c(normalizePath("../fastqs/PreCG_final_ATAC2020/DG_32_10X_ATAC_S-20-1357_GAP233/"), normalizePath("../fastqs/PreCG_final_ATAC2020/DG_S-20-0434_1_pool_GAP169"))

if (!dir.exists(CELLRANGER_ATAC_COUNT_BATCHTOOLS)) {
    reg <- makeRegistry(CELLRANGER_ATAC_COUNT_BATCHTOOLS)
} else {
    reg <- loadRegistry(CELLRANGER_ATAC_COUNT_BATCHTOOLS, writeable = TRUE)
}
reg$packages = c("batchtools", "tidyverse")

fastqs <- list.dirs(FASTQ_DIR, recursive = FALSE)
fastq_meta <- read_xlsx(FASTQ_META)

# FORCE_CELLS <- 10000
# Resource spec for batchtools. Use --localcores and --localmem to limit cellranger to job limits
RESOURCES <- list(
    ncpus = 32,
    memory = 128, # slurm/sge interprets as *giga*bytes, cellranger as *gibi*bytes
    walltime = 86400
)

fastq_tb <- tibble(path = fastqs, sample_id = basename(fastqs)) %>%
    inner_join(fastq_meta, by = c("sample_id" = "ATAC_fastq_name"))

fastq_tb <- fastq_tb %>%
    mutate(
        cmd = pmap(., function(...) {
            row <- list(...)
            cmd <- str_glue(
                "cd {CELLRANGER_ATAC_COUNT_DIR} && {CELLRANGER_ATAC_BIN} count",
                "--id={row$sample_id}",
                "--sample={row$sample_id}",
                "--reference={CELLRANGER_ATAC_REFERENCE}",
                "--fastqs={row$path}",
                "--localcores={RESOURCES$ncpus}",
                "--localmem={floor(RESOURCES$memory * GIGA_TO_GIBI * 0.9)}",
                "--force-cells={row$loaded}",
                .sep = " "
            )
        })
    )



check_paths <- map(fastqs, function(path) {
    sample_id <- basename(path)
    file.path(CELLRANGER_ATAC_COUNT_DIR, sample_id, "outs", "fragments.tsv.gz")
})

check_counts <- function(sample_id) {
    return(file.exists(file.path(CELLRANGER_ATAC_COUNT_DIR, sample_id, "outs", "fragments.tsv.gz")))
}

clearRegistry()
batchExport(list(check_counts = check_counts, CELLRANGER_ATAC_COUNT_DIR = CELLRANGER_ATAC_COUNT_DIR))
ids <- batchMap(function(cmd, sample_id) {
    system(cmd)
    if (!check_counts(sample_id)) {
        stop(str_glue("[EE] cellranger did not complete -- output file not found: {sample_id}/outs/fragments.tsv"))
    }
}, args = list(cmd = fastq_tb$cmd, sample_id = fastq_tb$sample_id))
setJobNames(ids, fastq_tb$sample_id)

done_jobs <- check_counts(fastq_tb$sample_id)

ids <- ids[!done_jobs,]
submitJobs(ids, resources = RESOURCES)
waitForJobs()
