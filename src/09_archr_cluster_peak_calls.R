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
    addArchRThreads(6)
    start_time <- Sys.time() 
    project <- loadArchRProject(proj_dir)
    subset_cells <- project@cellColData %>%
        as_tibble(rownames = "cell_ids") %>%
        filter(!Sample %in% c("I1_7", "i3_6_at", "P1_7_at1_7")) %>%
        pluck("cell_ids")
    project <- project[subset_cells, ]
    marker_tb <- read_csv(celltype_markers)

    cluster_counts <- as_tibble(project@cellColData, rownames = "cell_ids") %>%
        group_by(Clusters, Clinical.Dx) %>%
        summarize(n = n())

    minReplicates <- 4
    maxReplicates <- 10
    group_coverage_params <- cluster_counts %>%
        dplyr::rename(test_cluster = Clusters) %>%
        filter(test_cluster %in% c("C2", "C7")) %>%
        group_by(test_cluster) %>%
        summarize(minCell = floor(min(n) / minReplicates), maxCell = floor(max(n) / 2)) %>%
        mutate(
            out_dir = file.path(paste0(out_dir, "_", test_cluster)),
            minReplicates = minReplicates,
            maxReplicates = maxReplicates,
            groupBy = "Clinical.Dx"
        )

    print(group_coverage_params)

    pseudobulk_projects <- pmap(group_coverage_params, function(...) {
        print(list(...))
        return(tryCatch(group_coverage_worker(proj_dir, ...), error = print))
    })
    
    per_dx_features <- pmap(list(pseudobulk_projects, group_coverage_params$maxCell), function(project, maxCell) {
        print(maxCell)
        tryCatch(per_dx_marker_features(project, maxCell), error = print)
    })
    pwalk(list(pseudobulk_projects, per_dx_features), function(project, marker_list) {
        tryCatch(plot_dx_markers(project, marker_list, marker_tb), error = print)
    })

    per_dx_feature_marker_list <- pmap(list(pseudobulk_projects, per_dx_features), function(project, marker_list) {
        tryCatch(save_dx_features(project, marker_list), error = print)
    })

    end_time <- Sys.time() 
    print(end_time - start_time)
}

group_coverage_worker <- function(proj_dir, test_cluster, minCell, maxCell, minReplicates, maxReplicates, groupBy, out_dir) {
    writeMsg("group_coverage_worker")
    writeMsg(str_glue("load {proj_dir}"))
    project <- loadArchRProject(proj_dir, showLogo = FALSE)
    meta <- as_tibble(project@cellColData, rownames = "cell_ids") %>% filter(Clusters == test_cluster)
    project <- project[meta$cell_ids, ]
    project <- saveArchRProject(project, out_dir, dropCells = FALSE)

    writeMsg("group coverages")
    project <- addGroupCoverages(project, groupBy = groupBy, minCell = minCell, maxCell = maxCell, minReplicates = minReplicates, maxReplicates = maxReplicates)

    macs2_path <- findMacs2()
    writeMsg(str_glue("exec {macs2_path}"))
    project <- addReproduciblePeakSet(
        project,
        groupBy = groupBy,
        pathToMacs2 = macs2_path,
        genomeSize = 2467481108
    )

    writeMsg("add peaks")
    project <- addPeakMatrix(project)
    saveArchRProject(project, out_dir, load = FALSE)
    return(project)
}

per_dx_marker_features <- function(project, maxCells) { 
    proj_dir <- project@projectMetadata$outputDirectory
    writeMsg(str_glue("output to {proj_dir}"))
    dx_wk <- function(dx) {
        writeMsg(str_glue("per_dx_marker_features {dx}"))
        peakMarkers <- getMarkerFeatures(
            project,
            useMatrix = "PeakMatrix",
            useGroups = dx,
            bgdGroups = "Control",
            groupBy = "Clinical.Dx",
            maxCells = maxCells,
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = "wilcoxon"
        )
    }
    dx_features <- map(c("AD", "bvFTD", "PSP-S"), ~tryCatch(dx_wk(.), error = print))
    saveRDS(dx_features, file.path(proj_dir, "marker_features_dx.rds"))
    setNames(dx_features, c("AD", "bvFTD", "PSP-S"))
}

plot_dx_markers <- function(project, dx_marker_list, marker_tb) {
    proj_dir <- project@projectMetadata$outputDirectory
    cutOff <- "FDR <= 0.1 & abs(Log2FC) >= 0.3"
    map(dx_marker_list, function(markers) {
        peakMarkerList <- getMarkers(markers, cutOff = cutOff, returnGR = TRUE)
        if (length(peakMarkerList) > 0) {    
            dx <- names(peakMarkerList)
            writeMsg(str_glue("volcano {dx}"))
            tryCatch({
                volcanoMarkers <- plotMarkers(seMarker = markers, name = dx, cutOff = cutOff, plotAs = "Volcano")
                plotPDF(volcanoMarkers, name = str_glue("Peak-Marker-Volcano_{basename(proj_dir)}_{dx}"), width = 8, height = 6, ArchRProj = project, addDOC = FALSE)
            }, error = print)
            writeMsg(str_glue("heatmap {dx}"))
            tryCatch({
                heatmapPeaks <- plotMarkerHeatmap(seMarker = markers, cutOff = cutOff, transpose = TRUE, plotLog2FC = TRUE)
                plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width - 8, height = 6, ArchRProj = project, addDOC = FALSE)

            }, error = print)
            writeMsg(str_glue("genome track peaks {dx}"))
            tryCatch({
                peakTracks <- plotBrowserTrack(
                    project,
                    groupBy = "Clinical.Dx",
                    geneSymbol = marker_tb$gene_symbol,
                    features = peakMarkerList,
                    upstream = 50000,
                    downstream = 50000
                )
                plotPDF(peakTracks, name = str_glue("peak_marker_track_{dx}"), width = 8, height = 6, ArchRProj = project, addDOC = FALSE)
            }, error = print)
        }
    })
    
}

save_dx_features <- function(project, dx_marker_list) {
    proj_dir <- project@projectMetadata$outputDirectory
    cutOff <- "FDR <= 0.1 & abs(Log2FC) >= 0.3"
    map(dx_marker_list, function(markers) {
        peakMarkerList <- getMarkers(seMarker = markers, cutOff = cutOff, returnGR = TRUE)
        dx <- names(peakMarkerList)
        marker_rds_path <- file.path(proj_dir, str_glue("marker_feature_list_{dx}.rds"))
        writeMsg(str_glue("saveRDS({marker_rds_path})"))
        saveRDS(peakMarkerList[[1]], marker_rds_path)

        marker_csv_path <- file.path(proj_dir, str_glue("marker_feature_list_{dx}.csv"))
        writeMsg(str_glue("write_csv({marker_csv_path})"))
        write_csv(as.data.frame(peakMarkerList[[1]]), marker_csv_path)
        return(peakMarkerList[[1]])
    })
}

pwalk(list(archr_project, out_archr_project), function(proj_dir, out_dir) {
    writeMsg(str_glue("main({proj_dir}, {out_dir})"))
    main(proj_dir, out_dir)
})
