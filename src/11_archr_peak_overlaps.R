# extract ddifferentially-available promoter peaks that overlap with gene_list.
base_dir <- normalizePath("../data/archr/atac-2020-all")
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "circlize", "scales", "RColorBrewer", "ComplexHeatmap")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))
archr_project <- file.path(base_dir, c("chrom_var_precg_C2", "chrom_var_precg_C7", "chrom_var_insula_C2", "chrom_var_insula_C7"))

setBold <- system("tput rev", intern = TRUE)
setNorm <- system("tput sgr0", intern = TRUE)

gene_list <- c("RORB", "PU.1", "NR3C1", "RUNX1", "RUNX2")

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(setBold, msg, setNorm), side = "both"))
}

main <- function() {
    dir.create(file.path(base_dir, "11_archr_peak_overlaps"), showWarnings = FALSE)
    args_tb <- tibble(archr_project = archr_project) %>%
        mutate(peak_data = map(archr_project, get_peak_data))
    args_tb <- args_tb %>%
        full_join(tibble(gene = gene_list), by = character())
    glimpse(args_tb)
    args_tb <- mutate(args_tb, motif_overlap_list = pmap(list(args_tb$peak_data, args_tb$gene), peak_motif_overlap))
    args_tb <- mutate(args_tb, motif_overlap_tb = map(motif_overlap_list, filter_promoter_overlaps))

    motif_overlap_tb <- args_tb %>% 
        select(-motif_overlap_list, -peak_data) %>% 
        unnest(motif_overlap_tb) %>% 
        mutate(archr_project = basename(archr_project))

    write_csv(motif_overlap_tb, file.path(base_dir, "11_archr_peak_overlaps", "motif_overlap_tb.csv"))
}
get_peak_data <- function(proj_dir) {
    writeMsg(proj_dir)
    project <- loadArchRProject(proj_dir)
    consensus_peaks <- getPeakSet(project)
    motif_peak_matrix <- assays(getMatches(project))[["matches"]]
    dx_peak_markers <- readRDS(file.path(proj_dir, str_glue("marker_features_dx.rds")))
    return(list(peaks = consensus_peaks, peak_matrix = motif_peak_matrix, peak_markers = dx_peak_markers))
}

peak_motif_overlap <- function(peak_data, gene_name) {
    is_gene_tf <- str_detect(colnames(peak_data$peak_matrix), fixed(gene_name))
    is_overlapping_tf <- which(peak_data$peak_matrix[, is_gene_tf])
    print(length(is_overlapping_tf))
    tfpeaks <- peak_data$peaks[is_overlapping_tf]

    # per-dx differential peaks
    diff_dx_overlaps <- map(peak_data$peak_markers, function(peak_markers) {
        diff_dxpeak <- getMarkers(peak_markers, cutOff = "FDR <= 0.1", returnGR = TRUE)
        overlap <- GenomicRanges::findOverlaps(diff_dxpeak[[1]], tfpeaks)
        diff_tfpeak <- tfpeaks[subjectHits(overlap)]
    })
    setNames(diff_dx_overlaps, c("AD", "bvFTD", "PSP-S"))
}

filter_promoter_overlaps <- function(peak_motif_overlaps) {
    map(peak_motif_overlaps, function(x) {
        x %>%
            as_tibble(rownames = "x") %>%
            filter(peakType == "Promoter")
    }) %>%
    setNames(nm = names(peak_motif_overlaps)) %>%
    bind_rows(.id = "dx_test")
}

if (!interactive()) {
    main()
}
