# peak overlap iwth RORB.
base_dir <- normalizePath("../data/archr/atac-2020-all")
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "circlize", "scales", "RColorBrewer", "ComplexHeatmap")
archr_project <- file.path(base_dir, c("chrom_var_precg_C2", "chrom_var_precg_C7", "chrom_var_insula_C2", "chrom_var_insula_C7"))

setBold <- system("tput rev", intern = TRUE)
setNorm <- system("tput sgr0", intern = TRUE)

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(setBold, msg, setNorm), side = "both"))
}

main <- function() {
    dir.create(file.path(base_dir, "11_archr_peak_overlaps"), showWarnings = FALSE)
    peak_motif_overlaps <- pmap(list(archr_project), peak_motif_overlap)
    motif_overlap_tb <- pmap(list(archr_project, peak_motif_overlaps), filter_overlaps) %>%
        setNames(nm = basename(archr_project)) %>%
        bind_rows(.id = "project_name")

    write_csv(motif_overlap_tb, file.path(base_dir, "11_archr_peak_overlaps", "rorb_overlap_tb.csv"))
}

peak_motif_overlap <- function(proj_dir) {
    writeMsg(proj_dir)
    project <- loadArchRProject(proj_dir)
    consensus_peaks <- getPeakSet(project)
    motif_peak_matrix <- assays(getMatches(project))[["matches"]]
    is_overlapping_rorb <- motif_peak_matrix[, "RORB_443"]
    print(length(is_overlapping_rorb))
    tfpeaks <- consensus_peaks[is_overlapping_rorb]
    dx_peak_markers <- readRDS(file.path(proj_dir, str_glue("marker_features_dx.rds")))

    # per-dx differential peaks
    diff_dx_overlaps <- map(dx_peak_markers, function(peak_markers) {
        diff_dxpeak <- getMarkers(peak_markers, cutOff = "FDR <= 0.1", returnGR = TRUE)
        overlap <- GenomicRanges::findOverlaps(diff_dxpeak[[1]], tfpeaks)
        diff_tfpeak <- tfpeaks[subjectHits(overlap)]
    })
    setNames(diff_dx_overlaps, c("AD", "bvFTD", "PSP-S"))
}

filter_overlaps <- function(proj_dir, peak_motif_overlaps) {
    map(peak_motif_overlaps, function(x) {
        x %>%
            as_tibble(rownames = "x") %>%
            filter(peakType == "Promoter")
    }) %>%
    setNames(nm = names(peak_motif_overlaps)) %>%
    bind_rows(.id = "dx_test")
}
