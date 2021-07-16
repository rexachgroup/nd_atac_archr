# export peak calls to BED format for external use.

liblist <- c("tidyverse", "Seurat", "ArchR", "rtracklayer")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))

base_dir <- normalizePath("../data/archr/atac-2020-all")
out_dir <- file.path(base_dir, "11_archr_peaks_export")

project_names <- c("precg_C7", "insula_C7", "precg_C2", "insula_C2")
archr_project <- setNames(file.path(base_dir, paste0("peak_calling_", project_names)), project_names)

RESOURCES <- list(
    ncpus = 1,
    memory = 4,
    walltime = 86400,
    partition = "bigmem"
)


main <- function() {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    addArchRThreads(16)
    pmap(list(archr_project, project_names), export_peaks)
    #export_peaks(archr_project[[2]], project_names[[2]])
}

export_peaks <- function(archr_project, project_name) { 
    dir.create(file.path(out_dir), showWarnings = FALSE, recursive = TRUE)
    project <- loadArchRProject(archr_project)
    rtracklayer::export(getPeakSet(project), file.path(out_dir, str_glue("{project_name}.bed")), "BED")
}

if (!interactive()) {
    main()
}
