liblist <- c("tidyverse", "batchtools", "readxl", "ArchR", "Seurat")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

CELLRANGER_COUNT_DIRS <- c(
    "/geschwindlabshares/lchenprj01/nucseq_combined/data/cellranger-atac-count/precg-atac-2020", 
    "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/ATAC"
)
SAMPLE_META <- "/geschwindlabshares/lchenprj01/nucseq_combined/data/snATAC_metadata_summary_2021_d.xlsx"
out_dir <- "../data/archr/atac-2020-all"
data_dir <- normalizePath(file.path(out_dir, "data"))
plot_dir <- normalizePath(file.path(out_dir, "plot"))

main <- function() {
    addArchRThreads(threads = 32)
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    sample_meta <- read_xlsx(SAMPLE_META)
    sample_dirs <- list.dirs(CELLRANGER_COUNT_DIRS, recursive = FALSE)
    sample_tb <- tibble(sample_dirs) %>%
        mutate(
            library_id = basename(sample_dirs),
            fragments = file.path(sample_dirs, "outs", "fragments.tsv.gz"),
            cells = file.path(sample_dirs, "outs", "singlecell.csv"),
        ) %>%
        filter(file.exists(fragments), file.exists(cells))

    orig_dir <- getwd()
    setwd(data_dir)
    addArchRGenome("hg38")
    sample_tb <- mutate(sample_tb,
        arrow = file.path(data_dir, paste0(library_id, ".arrow"))
    ) %>%
    filter(file.exists(arrow))
    doublet_scores <- addDoubletScores(
        input = sample_tb$arrow,
        k = 10,
        knnMethod = "UMAP",
        LSIMethod = 1
    )

    project_dir <- file.path(out_dir, "preprocess")
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    project <- ArchRProject(
        ArrowFiles = sample_tb$arrow,
        outputDirectory = project_dir,
        copyArrows = FALSE
    )
    saveArchRProject(ArchRProj = project, outputDirectory = project_dir)
    saveRDS(sample_tb, file.path(project_dir, "sample_tb.rds"))
}

if(!interactive()) {
    main()
}
