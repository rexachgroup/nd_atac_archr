# Call peaks and test differentially per subcluster.
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
cluster_args <- file.path(base_dir, "07_archr_harmony_subclustering", "cluster_args_tb.rds")
out_dir <- file.path(base_dir, "09_archr_subcluster_peak_calling")
batchtools_dir <- file.path(out_dir, "batchtools")
celltype_markers <- normalizePath("../ext/20210128_cell_markers_noependymial.csv")
orig_dir <- getwd()

setBold <- system("tput rev", intern = TRUE)
setNorm <- system("tput sgr0", intern = TRUE)

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(" ", msg, " "), side = "both", pad = "-"))
}

RESOURCES <- list(
    ncpus = 8,
    memory = 32,
    walltime = 86400
)

main <- function() {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    marker_tb <- read_csv(celltype_markers)
    
    if (!dir.exists(batchtools_dir)) {
        reg <- makeRegistry(batchtools_dir)
    } else {
        reg <- loadRegistry(batchtools_dir, writeable = TRUE)
    }

    cluster_args_tb <- readRDS(cluster_args) %>%
        dplyr::rename(proj_dir = out_path) %>%
        mutate(out_dir = file.path(out_dir, paste("peak_call", proj_name, sep = "_")))

    reg$packages <- liblist
    batchExport(mget(ls()))
    clearRegistry()
    ids <- batchMap(call_worker,
        args = list(proj_dir = cluster_args_tb$proj_dir, out_dir = cluster_args_tb$out_dir))
    submitJobs(ids, RESOURCES)
    waitForJobs()

    saveRDS(cluster_args_tb, file.path(out_dir, "cluster_args_tb.rds"))

    clearRegistry()
    ids <- batchMap(marker_feature_worker, args = cluster_args_tb, reg = reg)
    submitJobs(ids, resources = RESOURCES, reg = reg)
    waitForJobs()
    cluster_args_tb <- mutate(cluster_args_tb, cluster_peak_features = map(ids$job.id, loadResult))

    pwalk(cluster_args_tb, function(...) {
            cr <- list(...)
            plot_dir <- file.path(cr$out_dir, "Plots")
            project <- loadArchRProject(cr$out_dir)

            plot_markers(project, cr$cluster_peak_features, marker_tb)
            save_markers(project, cr$cluster_peak_features)
        })
    
    pmap(cluster_args_tb, function(...) {
        cr <- list(...)    
        filelist <- c(
            list.files(file.path(cr$out_dir, "Plots"), pattern = "*.pdf", full.names = TRUE),
            list.files(file.path(cr$out_dir), pattern = "marker_feature_list_.*csv", full.names = TRUE, recursive = TRUE)
        )
        print(filelist)
        unlink(file.path(out_dir, cr$proj_name), recursive = TRUE)
        dir.create(file.path(out_dir, cr$proj_name), showWarnings = FALSE)
        file.copy(filelist, file.path(out_dir, cr$proj_name), overwrite = TRUE)
    })

}

call_worker <- function(proj_dir, out_dir) {
    addArchRGenome("hg38")
    addArchRThreads(16)

    writeMsg(str_glue("copy {proj_dir} {out_dir}"))
    project <- loadArchRProject(path = proj_dir)
    project <- saveArchRProject(project, out_dir, load = TRUE)
    orig_dir <- getwd()
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

    setwd(orig_dir)
    saveArchRProject(project, load = FALSE)
}

### per-cluster peak regions
per_cluster_peak_marker_features <- function(project, maxCells) { 
    proj_dir <- project@projectMetadata$outputDirectory
    writeMsg(str_glue("output to {proj_dir}"))
    writeMsg(str_glue("per_cluster_peak_marker_features"))
    peakMarkers <- getMarkerFeatures(
        project,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters",
        maxCells = maxCells,
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )
    saveRDS(peakMarkers, file.path(proj_dir, "marker_features.rds"))
    return(peakMarkers)
}

marker_feature_worker <- function(...) {
    cr <- list(...)
    #unlink(file.path(cr$out_dir, "Plots"), recursive = TRUE)
    plot_dir <- file.path(cr$out_dir, "Plots")
    #dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    project <- loadArchRProject(cr$out_dir)

    writeMsg("dx marker features")
    cluster_counts <- as_tibble(project@cellColData, rownames = "cell_ids") %>%
        group_by(Clusters, Clinical.Dx) %>%
        summarize(n = n())
    maxCells = floor(max(cluster_counts$n) / 2)
    writeMsg(str_glue("maxCells: {maxCells}"))

    per_cluster_peak_marker_features(project, maxCells)
}

plot_markers <- function(project, marker_list, marker_tb) {
    writeMsg("plot_markers")
    proj_dir <- project@projectMetadata$outputDirectory
    orig_dir <- getwd()
    setwd(proj_dir)
    cutOff <- "FDR <= 0.1 & abs(Log2FC) >= 0.3"
    peakMarkerList <- getMarkers(marker_list, cutOff = cutOff, returnGR = TRUE)
    if (length(peakMarkerList) > 0) {
        tryCatch({
            volcanoMarkers <- plotMarkers(seMarker = marker_list, cutOff = cutOff, plotAs = "Volcano")
            plotPDF(volcanoMarkers, name = str_glue("Peak-Marker-Volcano_{basename(proj_dir)}_cluster"), width = 8, height = 6, ArchRProj = project, addDOC = FALSE)
        }, error = print)
        writeMsg(str_glue("heatmap"))
        tryCatch({
            heatmapPeaks <- plotMarkerHeatmap(seMarker = marker_list, cutOff = cutOff, transpose = TRUE, plotLog2FC = TRUE)
            plotPDF(heatmapPeaks, name = str_glue("Peak-Marker-Heatmap_cluster"), width = 8, height = 6, ArchRProj = project, addDOC = FALSE)

        }, error = print)
        writeMsg(str_glue("genome track peaks"))
        tryCatch({
            peakTracks <- plotBrowserTrack(
                project,
                groupBy = "Clusters",
                geneSymbol = marker_tb$gene_symbol,
                features = peakMarkerList,
                upstream = 50000,
                downstream = 50000
            )
            plotPDF(peakTracks, name = str_glue("peak_marker_track_cluster"), width = 8, height = 6, ArchRProj = project, addDOC = FALSE)
        }, error = print)
    }
    setwd(orig_dir)
}

save_markers <- function(project, marker_list) {
    writeMsg("save_markers")
    proj_dir <- project@projectMetadata$outputDirectory
    orig_dir <- getwd()
    setwd(proj_dir)
    cutOff <- "FDR <= 0.1 & abs(Log2FC) >= 0.3"
    peakMarkerList <- getMarkers(seMarker = marker_list, cutOff = cutOff, returnGR = TRUE)
    dx <- names(peakMarkerList)
    marker_rds_path <- file.path(proj_dir, str_glue("PeakCalls/marker_feature_list_cluster.rds"))
    writeMsg(str_glue("saveRDS({marker_rds_path})"))
    saveRDS(peakMarkerList[[1]], marker_rds_path)

    marker_csv_path <- file.path(proj_dir, str_glue("PeakCalls/marker_feature_list_cluster.csv"))
    writeMsg(str_glue("write_csv({marker_csv_path})"))
    write_csv(as.data.frame(peakMarkerList[[1]]), marker_csv_path)
    setwd(orig_dir)
    return(peakMarkerList[[1]])
}

### differential dx peak regions
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
    dx_features <- setNames(dx_features, c("AD", "bvFTD", "PSP-S"))
    saveRDS(dx_features, file.path(proj_dir, "marker_features_dx.rds"))
    return(dx_features)
}

plot_dx_markers <- function(project, dx_marker_list, marker_tb) {
    writeMsg("plot_dx_markers")
    proj_dir <- project@projectMetadata$outputDirectory
    orig_dir <- getwd()
    setwd(proj_dir)
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
                plotPDF(heatmapPeaks, name = str_glue("Peak-Marker-Heatmap_{dx}"), width = 8, height = 6, ArchRProj = project, addDOC = FALSE)

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
    setwd(orig_dir)
    
}

save_dx_features <- function(project, dx_marker_list) {
    writeMsg("save_dx_features")
    proj_dir <- project@projectMetadata$outputDirectory
    orig_dir <- getwd()
    setwd(proj_dir)
    cutOff <- "FDR <= 0.1 & abs(Log2FC) >= 0.3"
    map(dx_marker_list, function(markers) {
        peakMarkerList <- getMarkers(seMarker = markers, cutOff = cutOff, returnGR = TRUE)
        dx <- names(peakMarkerList)
        marker_rds_path <- file.path(proj_dir, str_glue("PeakCalls/marker_feature_list_{dx}.rds"))
        writeMsg(str_glue("saveRDS({marker_rds_path})"))
        saveRDS(peakMarkerList[[1]], marker_rds_path)

        marker_csv_path <- file.path(proj_dir, str_glue("PeakCalls/marker_feature_list_{dx}.csv"))
        writeMsg(str_glue("write_csv({marker_csv_path})"))
        write_csv(as.data.frame(peakMarkerList[[1]]), marker_csv_path)
        return(peakMarkerList[[1]])
    })
    setwd(orig_dir)
}

if (!interactive()) {
    main()
}
