source("00_CONSTANTS.R")
library(batchtools)
library(tidyverse)

if (!dir.exists(CELLRANGER_ATAC_COUNT_BATCHTOOLS)) {
    reg <- makeRegistry(CELLRANGER_ATAC_COUNT_BATCHTOOLS)
} else {
    reg <- loadRegistry(CELLRANGER_ATAC_COUNT_BATCHTOOLS, writeable = TRUE)
}


fastqs <- list.dirs(FASTQ_DIR, recursive = FALSE)

FORCE_CELLS <- 10000
# Resource spec for batchtools. Use --localcores and --localmem to limit cellranger to job limits
RESOURCES <- list(
    ncpus = 32,
    memory = 128, # slurm/sge interprets as *giga*bytes, cellranger as *gibi*bytes
    walltime = 86400
)

cmdlist <- map(fastqs, function(path) {
    sample_id <- basename(path)
    cmd <- str_glue("cd {CELLRANGER_ATAC_COUNT_DIR} && {CELLRANGER_ATAC_BIN} count",
                    "--id={sample_id}",
                    "--sample={sample_id}",
                    "--reference={CELLRANGER_ATAC_REFERENCE}",
                    "--fastqs={path}",
                    "--localcores={RESOURCES$ncpus}",
                    "--localmem={floor(RESOURCES$memory * GIGA_TO_GIBI * 0.9)}",
                    .sep = " ")
    if (exists("FORCE_CELLS")) {
        cmd <- str_glue(cmd, " --force-cells={FORCE_CELLS}")
    }
    return(cmd)
})

check_paths <- map(fastqs, function(path) {
    sample_id <- basename(path)
    file.path(CELLRANGER_ATAC_COUNT_DIR, sample_id, "outs", "fragments.tsv.gz")
})

check_counts <- function(sample_id) {
    return(file.exists(file.path(CELLRANGER_ATAC_COUNT_DIR, sample_id, "outs", "fragments.tsv.gz")))
}

clearRegistry()
ids <- batchMap(function(cmd, check_path) {
    system(cmd)
    if (!file.exists(check_path)) {
        stop(str_glue("[EE] cellranger did not complete -- output file not found: {check_path}"))
    }
}, args = list(cmd = cmdlist, check_path = check_paths))
setJobNames(ids, basename(fastqs))

done_jobs <- map_lgl(basename(fastqs), function(sample_id) {
    check_counts(sample_id)
})

ids <- ids[!done_jobs,]
submitJobs(ids, resources = RESOURCES)
waitForJobs()
