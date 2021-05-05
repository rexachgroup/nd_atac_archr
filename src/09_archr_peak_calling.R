liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "ragg")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
setBold <- system("tput rev", intern = TRUE)
setNorm <- system("tput sgr0", intern = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, c("liger_integration_precg", "liger_integration_insula"))
out_archr_project <- file.path(base_dir, c("peak_calling_precg", "peak_calling_insula"))
celltype_markers <- normalizePath("../ext/20210128_cell_markers_noependymial.csv")

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(setBold, msg, setNorm), side = "both"))
}

main <- function(proj_dir, out_dir) {
    addArchRGenome("hg38")
    addArchRThreads(16)

    writeMsg(str_glue("copy {proj_dir} {out_dir}"))
    project <- loadArchRProject(path = proj_dir)
    project <- saveArchRProject(project, out_dir, load = TRUE)
    plot_dir <- file.path(out_dir, "Plots")
    orig_dir <- getwd()
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(out_dir)
    marker_tb <- read_csv(celltype_markers)

    writeMsg("group coverages")
    project <- addGroupCoverages(project, groupBy = "Clusters")

    macs2_path <- findMacs2()
    writeMsg(str_glue("exec {macs2_path}"))
    project <- addReproduciblePeakSet(
        project,
        groupBy = "Clusters",
        pathToMacs2 = macs2_path,
        genomeSize = 2467481108
    )

    writeMsg("add peaks")
    project <- addPeakMatrix(project)

    writeMsg("find peak markers")
    peakMarkers <- getMarkerFeatures(
        project,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )

    peakMarkerList <- getMarkers(peakMarkers, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
    saveRDS(peakMarkers, file.path(out_dir, "peakMarkers.rds"))
    saveRDS(peakMarkerList, file.path(out_dir, "peakMarkerList.rds"))
    writeMsg("heatmap peaks")
    heatmapPeaks <- plotMarkerHeatmap(seMarker = peakMarkers, cutOff = "FDR <= 0.01 & Log2FC >= 1", transpose = TRUE)
    plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = project, addDOC = FALSE)

    writeMsg("genome track peaks")
    peakTracks <- plotBrowserTrack(
        project,
        groupBy = "Clusters",
        geneSymbol = marker_tb$gene_symbol,
        features = peakMarkerList,
        upstream = 50000,
        downstream = 50000
    )
    plotPDF(peakTracks, name = "peak_marker_track", width = 8, height = 6, ArchRProj = project, addDOC = FALSE)

    saveArchRProject(project, load = FALSE)
    setwd(orig_dir)
}

pwalk(list(archr_project, out_archr_project), function(proj_dir, out_dir) {
    writeMsg(str_glue("main({proj_dir}, {out_dir})"))
    main(proj_dir, out_dir)
})
