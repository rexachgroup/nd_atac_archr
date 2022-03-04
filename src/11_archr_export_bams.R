# export cluster frags for external processing.

liblist <- c("tidyverse", "readxl", "Seurat", "ArchR")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = FALSE)))
base_dir <- normalizePath("../data/archr/atac-2020-all" )
out_dir <- file.path(base_dir, "11_archr_bam_export")

CELLRANGER_COUNT_DIRS <- c(
    "/geschwindlabshares/lchenprj01/nucseq_combined/data/cellranger-atac-count/precg-atac-2020", 
    "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/ATAC"
)

SAMPLE_META <- "/geschwindlabshares/lchenprj01/nucseq_combined/data/snATAC_metadata_summary_2021_d.xlsx"
SUBSET_BIN <- system("which subset-bam_linux", intern = TRUE)
SAMTOOLS_BIN <- system("which samtools", intern = TRUE)

clusters <- tribble(
    ~name,                  ~archr_path,                                        ~subclusters,       ~dx,
    "insula-ftd-c78",       file.path(base_dir, "liger_integration_insula"),    c("C7", "C8"),      "bvFTD",
    "insula-ftd-c10-c12",   file.path(base_dir, "liger_integration_insula"),    c("C10", "C12"),    "bvFTD",
    "precg-ad-c7",          file.path(base_dir, "liger_integration_precg"),     c("C7"),            "AD",
    "precg-ad-c2",          file.path(base_dir, "liger_integration_precg"),     c("C2"),            "AD",
    "precg-ftd-c2",         file.path(base_dir, "liger_integration_precg"),     c("C2"),            "bvFTD"
)

main <- function() {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    clusters <- clusters %>%
        mutate(
            bam_path = file.path(out_dir, str_glue("{name}.bam")),
            subset_dir = file.path(out_dir, str_glue("{name}"))
        ) %>%
        filter(!file.exists(bam_path))
    cluster_tb <- clusters %>% mutate(sample_barcodes = pmap(., prep_barcodes))
    archr_meta_tb <- cluster_tb$sample_barcodes
    sample_tb <- get_sample_tb(SAMPLE_META, CELLRANGER_COUNT_DIRS)
    
    map(archr_meta_tb, subset_bam_by_sample, sample_tb)
    map(archr_meta_tb, join_bams, out_dir)
    cleanup_subset_dir(cluster_tb)
}

prep_barcodes <- function(...) {
    cr <- list(...)
    dir.create(cr$subset_dir, showWarnings = FALSE)
    project <- loadArchRProject(cr$archr_path)
    meta <- project@cellColData %>%
        as_tibble(rownames = "cell_ids")
    
    dx_meta <- meta %>%
        group_by(Sample, Clinical.Dx) %>%
        filter(Clinical.Dx %in% cr$dx, Clusters %in% cr$subclusters) %>%
        group_nest %>%
        mutate(barcodes_path = file.path(cr$subset_dir, str_glue("{cr$name}_{Sample}_{Clinical.Dx}_barcodes.csv")))
    
    pwalk(dx_meta, function(...) {
        cr <- list(...)
        cell_ids_strip <- str_extract(cr$data$cell_ids, "(?<=#)(.+)")
        write_csv(as.data.frame(cell_ids_strip), cr$barcodes_path)
    })

    dx_meta$name <- cr$name

    return(dx_meta)
}

get_sample_tb <- function(meta_xlsx, cellranger_count_dirs) { 
    sample_meta <- read_xlsx(meta_xlsx)
    sample_dirs <- list.dirs(cellranger_count_dirs, recursive = FALSE)
    sample_tb <- tibble(sample_dirs) %>%
        mutate(
            library_id = basename(sample_dirs),
            fragments = file.path(sample_dirs, "outs", "fragments.tsv.gz"),
            bam = file.path(sample_dirs, "outs", "possorted_bam.bam"),
            cells = file.path(sample_dirs, "outs", "singlecell.csv"),
        ) %>%
        filter(file.exists(fragments), file.exists(cells))
    return(sample_tb)
}

subset_bam_by_sample <- function(archr_meta_tb, sample_tb) { 
    args_tb <- left_join(archr_meta_tb, sample_tb, by = c("Sample" = "library_id")) %>%
        mutate(
            subset_bam = file.path(out_dir, name, str_glue("{name}_{Sample}_{Clinical.Dx}_subset.bam")),
            subset_cmd = str_glue("{SUBSET_BIN} --bam {bam} --cell-barcodes {barcodes_path} --out-bam {subset_bam} --cores 8")
        )

    print(args_tb)
    args_tb %>%
        filter(file.exists(subset_bam)) %>%
        mutate(unlink = unlink(subset_bam))

    walk(args_tb$subset_cmd, function(x) {
        print(x)
        if (system(x) != 0) { stop() }
    })
}

join_bams <- function(archr_meta_tb, out_dir) {
    args_tb <- archr_meta_tb %>%
        mutate(subset_bam = file.path(out_dir, name, str_glue("{name}_{Sample}_{Clinical.Dx}_subset.bam"))) %>%
        group_by(Clinical.Dx) %>%
        group_nest %>%
        mutate(dx_aggr = pmap(., function(...) {
            cr <- list(...)
            return(str_glue(
                "{SAMTOOLS_BIN}",
                "cat",
                "-o {out_dir}/{unique(cr$data$name)}.bam",
                "--threads 8",
                "{paste(cr$data$subset_bam, collapse = ' ')}",
                .sep = " "
            ))
        }))
    print(args_tb)
    walk(args_tb$dx_aggr, function(x) {
        print(x)
        #system(x)
        if (system(x) != 0) { stop() }
    })
}

cleanup_subset_dir <- function(cluster_tb) {
    cluster_tb <- cluster_tb %>%
        filter(file.exists(bam_path))
    unlink(cluster_tb$subset_dir, recursive = TRUE)
}

if (!interactive()) main()
