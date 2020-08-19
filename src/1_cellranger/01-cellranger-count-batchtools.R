lapply(c("rprojroot",
         "batchtools",
         "tidyverse"),
       function(x) suppressPackageStartupMessages(require(x, character.only = TRUE)))

root <- has_file(".rprojroot")

# cellranger config
source(root$find_file("src/1_cellranger/00_CONSTANTS.R"))
FORCE_CELL_COUNT <- 7000
MEM_LIMIT_TOTAL <- 128
NCORES <- 16
MEM_LIMIT_HDATA_GB <- paste0(ceiling(MEM_LIMIT_TOTAL / NCORES), "G")

# output
OUT_ROOT <- normalizePath(root$find_file("data/01-cellranger-individual"))

# batchtools
REGISTRY <- normalizePath(root$find_file("batchtools/01-cellranger-count-individual"))

init_registry <- function(reg_path, clear = TRUE) {
    if (dir.exists(reg_path)) {
        reg <- loadRegistry(reg_path, writeable = TRUE)
    } else {
        reg <- makeRegistry(reg_path)
    }
    if (clear) clearRegistry()
    reg$packages <- c("batchtools")
    reg$cluster.functions <- makeClusterFunctionsSGE(template =
        as.character(root$find_file("src/rscript_hoffman.tmpl")),
        fs.latency = 1)
    return(reg)
}


libraries <- tibble(dir = list.dirs(FASTQ_DIR, recursive = FALSE))
libraries$name <- basename(libraries$dir)

cellranger_count_cmds <- libraries %>% pmap(function(name, dir) {
    out_dir <- file.path(OUT_ROOT, name)
    stringr::str_glue("rm -rf {out_dir} &&",
                      "cd {OUT_ROOT} &&",
                      "{CELLRANGER_BIN} count",
                      "--id={name}",
                      "--transcriptome={CELLRANGER_REFERENCE}",
                      "--fastqs={dir}",
                      "--sample={name}",
                      "--localcores={NCORES}",
                      "--localmem={MEM_LIMIT_TOTAL}",
                      "--force-cells={FORCE_CELL_COUNT}", .sep = " ")
})

cellranger_count_cmds <- setNames(cellranger_count_cmds, nm = libraries$name)

reg <- init_registry(REGISTRY, clear = TRUE)

ids <- batchMap(cellranger_count_cmds, fun = function(cmd) {
    system(cmd)
}, reg = reg)

getJobTable() %>% 
    pwalk(function(job.id, ...) {
        tags <- getUsedJobTags(job.id, reg)
        if (!libraries$name[[job.id]] %in% tags) {
            addJobTags(job.id, libraries$name[[job.id]], reg)
        }
    })
{}

ids <- getJobTable()[tags %in% c("PFC3_3")]

submitJobs(ids = ids,
           reg = reg,
           resources = list(
                h_data_str = MEM_LIMIT_HDATA_GB,
                h_rt = "24:00:00",
                pe_shared = NCORES,
                measure.memory = TRUE
            ))
waitForJobs(reg = reg)
# reduceResultsDataTable()
# getJobTable()
# failed <- getJobTable()
# submitJobs(ids = failed$job.id, reg = reg,
#            resources = list(h_data_str = MEM_LIMIT_HDATA_GB, h_rt = "24:00:00", pe_shared = NCORES, measure.memory = TRUE))
# waitForJobs()
