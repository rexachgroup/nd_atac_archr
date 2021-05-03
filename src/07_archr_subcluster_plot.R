liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "patchwork", "ragg")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

SAMPLE_META <- normalizePath("../data/snATAC_metadata_summary_2021_e.xlsx")
base_dir <- normalizePath("../data/archr/atac-2020-all")
cluster_args <- file.path(base_dir, "07_archr_harmony_subclustering", "cluster_args_tb.rds")
celltype_markers <- normalizePath("../ext/20210128_cell_markers_noependymial.csv")
out_dir <- file.path(base_dir, "07_archr_subcluster_plot")
batchtools <- file.path(out_dir, "batchtools")

RESOURCES <- list(
    ncpus = 4,
    memory = 16,
    walltime = 86400
)

main <- function() {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(batchtools)
    } else {
        reg <- loadRegistry(batchtools, writeable = TRUE)
    }

    cluster_args_tb <- readRDS(cluster_args) %>%
        rename("out_path" = "proj_dir")
    
    reg$packages <- liblist
    batchExport(mget(ls()))
    clearRegistry()
    ids <- batchMap(plot_worker,
        args = list(proj_dir = cluster_args_tb$proj_dir))
    submitJobs(ids, RESOURCES)
    waitForJobs()

    pmap(cluster_args_tb, function(...) {
        cr <- list(...)    
        filelist <- list.files(str_glue("{cr$proj_dir}/Plots/"), pattern = "*.png", full.names = TRUE)
        dir.create(file.path(out_dir, cr$proj_name), showWarnings = FALSE)
        file.copy(filelist, file.path(out_dir, cr$proj_name), overwrite = TRUE)
    })

}

plot_worker <- function(proj_dir) { 
    addArchRGenome("hg38")
    writeLines(str_glue("load {proj_dir}"))
    project <- loadArchRProject(path = proj_dir)
    setwd(proj_dir)
    marker_genescores <- readRDS(file.path(proj_dir, "marker_genescores.rds"))
    marker_tb <- read_csv(celltype_markers)

    # sample / cluster distribution heatmap
    cmat <- confusionMatrix(project$Clusters, project$Sample)

    cprop <- cmat / rowSums(cmat)
    cm_heatmap <- pheatmap(
        as.matrix(cprop),
        color = paletteContinuous("whiteBlue"),
        border_color = "black"
    )
    agg_png(file.path(proj_dir, "Plots", "cluster_sample_heatmap.png"), height = 0.25 * nrow(cprop), width = 0.2 * ncol(cprop), units = "in", res = 300)
    print(cm_heatmap)
    graphics.off()

    agg_png(file.path(proj_dir, "Plots", "sample_umap.png"), height = 10, width = 10, units = "in", res = 300)
    print(plotEmbedding(project, colorBy = "cellColData", name = "Sample", embedding = "UMAP"))
    graphics.off()
    agg_png(file.path(proj_dir, "Plots", "cluster_umap.png"), height = 10, width = 10, units = "in", res = 300)
    print(plotEmbedding(project, colorBy = "cellColData", name = "Clusters", embedding = "UMAP"))
    graphics.off()

    # marker genescores
    heatmap_genescores <- plotMarkerHeatmap(
        seMarker = marker_genescores,
        labelMarkers = marker_tb$gene_symbol,
        transpose = TRUE
    )
    agg_png(file.path(proj_dir, "Plots", "marker_heatmap.png"), width = 14, height = 7, units = "in", res = 300)
    print(heatmap_genescores)
    graphics.off()

    # marker UMAP
    marker_plot_list <- plotEmbedding(project,
        colorBy = "GeneScoreMatrix",
        name = marker_tb$gene_symbol,
        embedding = "UMAP",
        quantCut = c(0.01, 0.95),
        imputeWeights = NULL
    )
    graphics.off()
    marker_plot_fmt <- map(marker_plot_list, function(x) {
        x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank()
        )
    })
    plotlist_png_square(marker_plot_fmt, 
        file.path(proj_dir, "Plots", "markerGene_UMAP.png"), width = 2, height = 2, units = "in", res = 300)
    graphics.off()

    # metadata UMAP
    meta_plot_list <- plotEmbedding(project,
        colorBy = "cellColData",
        name = c("Clusters", "FinalSite", "Clinical.Dx", "region", "PrepBatch", "SeqBatch"),
        embedding = "UMAP",
        quantCut = c(0.01, 0.95),
        imputeWeights = NULL
    )
    graphics.off()
    meta_plot_fmt <- map(meta_plot_list, function(x) {
        x +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    })
    plotlist_png_square(meta_plot_fmt,
        file.path(proj_dir, "Plots", "meta_plot_UMAP.png"), width = 4, height = 4, units = "in", res = 300)
    graphics.off()

    # Individual variables.

    dx_plotlist <- plot_individual_meta_lvls(project, "Clinical.Dx", "ggHex")
    plotlist_png_square(dx_plotlist, file.path(proj_dir, "Plots", "clinical_dx.png"), width = 4, res = 300)

    region_plotlist <- plot_individual_meta_lvls(project, "region", "ggHex")
    plotlist_png_square(region_plotlist, file.path(proj_dir, "Plots", "region.png"), width = 4, res = 300)

    library_plotlist <- plot_individual_meta_lvls(project, "Sample", "ggPoint")
    plotlist_png_square(library_plotlist, file.path(proj_dir, "Plots", "library_id.png"), width = 4, res = 200)

    # per-region metadata UMAP
    walk(unique(project$region), function(s) {
        writeLines(s)
        region_cells <- which(project$region %in% s)
        print(length(region_cells))
        tmp_proj <- project[region_cells, ]
        print(tmp_proj)
        plot_list <- plotEmbedding(tmp_proj,
            colorBy = "cellColData",
            name = c("Clusters", "FinalSite", "Clinical.Dx", "region", "PrepBatch", "SeqBatch"),
            embedding = "UMAP",
            quantCut = c(0.01, 0.95),
            imputeWeights = NULL
        )
        graphics.off()
        plotlist_png_square(plot_list,
            file.path(proj_dir, "Plots", str_glue("meta_{s}_UMAP.png")), width = 4, res = 300) 
        graphics.off()
    })

    project@cellColData %>%
        as_tibble(rownames = "cell_id") %>%
        group_by(Sample, Clinical.Dx, Clusters) %>%
        summarize(n = n()) %>%
        write_csv(file.path(proj_dir, "subcluster_cell_counts.csv"))

    # marker tracks
    marker_tracks <- plotBrowserTrack(project,
        groupBy = "Clusters",
        geneSymbol = marker_tb$gene_symbol,
        upstream = 50000,
        downstream = 50000
    )
    plotPDF(marker_tracks,
        name = "marker_tracks.pdf",
        project = project,
        addDOC = FALSE,
        width = 4, height = 4)
    graphics.off()
}

plot_individual_meta_lvls <- function(project, column, type = "ggPoint") {
    meta <- project@cellColData
    writeLines(str_glue("{column}: plotting {length(unique(meta[[column]]))} levels"))
    umap_x <- getEmbedding(project)[[1]]
    umap_y <- getEmbedding(project)[[2]]

    plotlist <- map(unique(meta[[column]]), function(lvl) {
        writeLines(lvl)
        is_lvl <- meta[[column]] == lvl
        print(table(as.numeric(is_lvl)))
        if (type == "ggPoint") {
            gg <- ggPoint(umap_x, umap_y, as.numeric(is_lvl), highlightPoints = which(is_lvl), 
                discrete = TRUE, labelMeans = FALSE, labelAsFactors = TRUE, title = lvl, alpha = 0.5, keepAxis = FALSE)
        } else if (type == "ggHex") { 
            gg <- ggHex(umap_x, umap_y, as.numeric(is_lvl), labelMeans = FALSE, labelAsFactors = TRUE, title = lvl, keepAxis = FALSE, FUN = "mean", bins = 200, hexCut = c(0, 1))
        }
        
        gg +
            theme_ArchR(baseSize = 6.5) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
            theme(
                axis.text.x=element_blank(), 
                axis.ticks.x=element_blank(), 
                axis.text.y=element_blank(), 
                axis.ticks.y=element_blank()
            )

    })
    # ggPoint opens a graphics device that is never used
    graphics.off()
    return(plotlist)
}

plotlist_png_square <- function(plotlist, path, width = 4, height = width, units = "in", res = 200) {
    plot_ncol <- plot_nrow <- ceiling(sqrt(length(plotlist)))
    if (exists("agg_png")) {
        agg_png(path, width = width * plot_ncol, height = height * plot_nrow, units = units, res = res)
    } else {
        png(path, width = width * plot_ncol, height = height * plot_nrow, units = units, res = res)
    }
    print(wrap_plots(plotlist, ncol = plot_ncol, nrow = plot_nrow))
    dev.off()
}

if (!interactive()) {
    main()
}
