# intersect GRNBoost gene targets with peaks that overlap known tfs.

base_dir <- normalizePath("../data/archr/atac-2020-all")
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "circlize", "scales", "RColorBrewer", "ComplexHeatmap")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))

base_out_dir <- file.path(base_dir, "13_archr_motif_scenic_intersection")
motif_regulon_filter <- normalizePath("../data/microglia_preCGregulons_genes_GRNBoost_weights.csv")
#scenic_motif_overlap_tb <- file.path(base_dir, 

main <- function() { 
    dir.create(base_out_dir, showWarnings = FALSE)
    scenic_regulon_tb <- read_csv(motif_regulon_filter)
    scenic_tfs <- read_csv(file.path(base_dir, "11_archr_peak_overlaps", "motif_overlap_tb.csv"))
    scenic_join <- scenic_tfs %>%
        group_split(archr_project) %>%
        map(~inner_join(., scenic_regulon_tb, by = c("nearestGene" = "gene")))
    map(scenic_join, function(x) {
        join_path <- file.path(base_out_dir, str_glue("{unique(x$archr_project)}.csv"))
        print(join_path)
        write_csv(x, join_path)
    })
    
}

if (!interactive()) {
    main()
}
