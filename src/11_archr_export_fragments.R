# export fragments for external processing.

liblist <- c("tidyverse", "readxl", "Seurat", "ArchR")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))

base_dir <- normalizePath("../data/archr/atac-2020-all")
out_dir <- file.path(base_dir, "11_archr_fragments_export")

CELLRANGER_COUNT_DIRS <- c(
    "/geschwindlabshares/lchenprj01/nucseq_combined/data/cellranger-atac-count/precg-atac-2020", 
    "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/ATAC"
)
SAMPLE_META <- "/geschwindlabshares/lchenprj01/nucseq_combined/data/snATAC_metadata_summary_2021_d.xlsx"
SUBSET_BIN <- system("which subset-bam_linux", intern = TRUE)
SAMTOOLS_BIN <- system("which samtools", intern = TRUE)

project_names <- c("precg_C7", "insula_C7")
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
    subset_bams_dx(archr_project[[2]], project_names[[2]])
}

subset_bams_dx <- function(archr_project, project_name) { 
    dir.create(file.path(out_dir, project_name), showWarnings = FALSE)
    project <- loadArchRProject(archr_project)
    meta <- project@cellColData %>%
        as_tibble(rownames = "cell_ids")

    dx_meta <- meta %>%
        filter(Clusters == "C7") %>%
        group_by(Clinical.Dx, Sample) %>%
        group_nest %>%
        mutate(barcodes_path = file.path(out_dir, project_name, str_glue("{project_name}_{Clinical.Dx}_barcodes.csv")))

    pmap(dx_meta, function(...) {
        cr <- list(...)
        cell_ids_strip <- str_extract(cr$data$cell_ids, "(?<=#)(.+)")
        write_csv(as.data.frame(cell_ids_strip), cr$barcodes_path)
    })
    
    sample_meta <- read_xlsx(SAMPLE_META)
    sample_dirs <- list.dirs(CELLRANGER_COUNT_DIRS, recursive = FALSE)
    sample_tb <- tibble(sample_dirs) %>%
        mutate(
            library_id = basename(sample_dirs),
            fragments = file.path(sample_dirs, "outs", "fragments.tsv.gz"),
            bam = file.path(sample_dirs, "outs", "possorted_bam.bam"),
            cells = file.path(sample_dirs, "outs", "singlecell.csv"),
        ) %>%
        filter(file.exists(fragments), file.exists(cells))

    args_tb <- left_join(dx_meta, sample_tb, by = c("Sample" = "library_id")) %>%
        mutate(
            subset_bam = file.path(out_dir, project_name, str_glue("{project_name}_{Sample}_{Clinical.Dx}_subset.bam")),
            subset_cmd = str_glue("{SUBSET_BIN} --bam {bam} --cell-barcodes {barcodes_path} --out-bam {subset_bam} --cores 8")
        )

    walk(args_tb$subset_cmd, function(x) {
        print(x)
        if (system(x) != 0) { stop() }
    })
    
    args_tb <- args_tb %>%
        group_by(Clinical.Dx) %>%
        group_nest %>%
        mutate(dx_aggr = pmap(., function(...) {
            cr <- list(...)
            return(str_glue(
                "{SAMTOOLS_BIN}",
                "cat",
                "-o {out_dir}/{project_name}_{cr$Clinical.Dx}.bam",
                "--threads 8",
                "{paste(cr$data$subset_bam, collapse = ' ')}",
                .sep = " "
            ))
        }))

    walk(args_tb$dx_aggr, function(x) {
        print(x)
        if (system(x) != 0) { stop() }
    })
}
