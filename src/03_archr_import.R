liblist <- c("tidyverse", "batchtools", "readxl", "ArchR")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

CELLRANGER_COUNT_DIRS <- c(
    "/geschwindlabshares/lchenprj01/nucseq_combined/data/cellranger-atac-count/precg-atac-2020", 
    "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/ATAC"
)
SAMPLE_META <- "../data/snATAC_metadata_summary_2021_f.xlsx"
out_dir <- "../data/archr/atac-2020-all"
data_dir <- file.path(out_dir, "data")
plot_dir <- file.path(out_dir, "plot")

main <- function() {
    addArchRThreads(16)
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
    sample_arrowfiles <- createArrowFiles(
        inputFiles = sample_tb$fragments,
        sampleNames = sample_tb$library_id,
        minTSS = 2,
        minFrags = 1000,
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
        subThreading = FALSE,
        force = TRUE
    )
}

if (!interactive()) {
    main()
}
