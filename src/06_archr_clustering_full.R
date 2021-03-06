liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

SAMPLE_META <- normalizePath("../data/snATAC_metadata_summary_2021_e.xlsx")
base_dir <- normalizePath("../data/archr/atac-2020-all")
data_dir <- file.path(base_dir, "data")
plot_dir <- file.path(base_dir, "plot")
archr_project <- file.path(base_dir, "preprocess")
sample_table <- file.path(base_dir, "preprocess", "sample_tb.rds")
out_archr_project <- file.path(base_dir, "clustering_full")
celltype_markers <- normalizePath("../ext/20210128_cell_markers_noependymial.csv")

main <- function() {
    #addArchRThreads(threads = 6)
    addArchRGenome("hg38")
    writeLines(str_glue("load {archr_project}"))
    project <- loadArchRProject(path = archr_project)
    orig_dir <- getwd()
    setwd(out_archr_project)
    project@projectMetadata$outputDirectory <- out_archr_project
    dir.create(file.path(out_archr_project, "Plots"), showWarning = FALSE) 
    marker_tb <- read_csv(celltype_markers)

    # Copy metadata from subject to cellColData.
    sample_meta <- read_xlsx(SAMPLE_META)
    cell_meta <- as_tibble(project@cellColData, rownames = "cell_id")
    meta_join <- left_join(cell_meta, sample_meta, by = c("Sample" = "ATAC_fastq_name"))
    stopifnot(rownames(project) == meta_join$cell_id)
    #project@cellColData <- DataFrame(meta_join, row.names = meta_join$cell_id)
    new_cols <- setdiff(colnames(meta_join), colnames(cell_meta))
    project@cellColData[, new_cols] <- meta_join[, new_cols]

    # Copy project.
    #project <- saveArchRProject(project, out_archr_project)

    # Doublet removal before clustering.
    filterDoublets(project)
    
    # additional filtering -- TSSEnrichment > 4, BlacklistRatio < 0.1
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")
    filtered_cells <- cell_data %>%
        filter(TSSEnrichment > 4, BlacklistRatio < 0.1) %>%
        pluck("cell_id")
    writeLines(str_glue("{length(filtered_cells)} / {nrow(cell_data)} cells after filtering"))

    writeLines("subset")
    project <- project[filtered_cells, ]
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")

    # LSI
    project <- addIterativeLSI(project,
        useMatrix = "TileMatrix",
        name = "IterativeLSI",
        iterations = 4,
        clusterParams = list(
            resolution = c(0.1, 0.2, 0.4),
            sampleCells = 10000,
            n.start = 10
        ),
        varFeatures = 15000,
        dimsToUse = 1:30,
        outDir = out_archr_project,
        force = TRUE
    )
    #     project <- addHarmony(project,
    #         force = TRUE
    #     )
    #saveArchRProject(project, out_archr_project, load = FALSE)
    print(project@reducedDims)
    # Clusters
    project <- addClusters(project,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        resolution = 0.8,
        force = TRUE
    )
    #saveArchRProject(project, out_archr_project, load = FALSE)

    # sample / cluster distribution heatmap
    cmat <- confusionMatrix(project$Clusters, project$Sample)

    cprop <- cmat / rowSums(cmat)
    cm_heatmap <- pheatmap(
        as.matrix(cprop),
        color = paletteContinuous("whiteBlue"),
        border_color = "black"
    )
    pdf(file.path(out_archr_project, "Plots", "cluster_sample_heatmap_50k.pdf"), height = 0.25 * nrow(cprop), width = 0.2 * ncol(cprop))
    print(cm_heatmap)
    dev.off()

    # cell UMAP
    project <- addUMAP(project,
        reducedDims = "IterativeLSI",
        name = "UMAP",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine"
    )

    umap_samples <- plotEmbedding(project, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
    umap_clusters <- plotEmbedding(project, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
    plotPDF(
        umap_samples,
        umap_clusters,
        name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = project,
        addDOC = TRUE,
        width = 10,
        height = 10
    )
 
    # marker genes
    marker_genescores <- getMarkerFeatures(project,
        useMatrix = "GeneScoreMatrix",
        groupBy = "Clusters",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )
    saveRDS(marker_genescores, file.path(out_archr_project, "marker_genescores.rds"))
    marker_list <- getMarkers(marker_genescores, cutOff = "FDR <=0.01 & Log2FC >=1.25")

    heatmap_genescores <- plotMarkerHeatmap(
        seMarker = marker_genescores,
        labelMarkers = marker_tb$gene_symbol,
        transpose = TRUE
    )
    dev.off()
    png(file.path(out_archr_project, "Plots", "marker_heatmap.png"), width = 14, height = 7, units = "in", res = 100)
    print(heatmap_genescores)
    dev.off()

    # marker UMAP
    project <- addImputeWeights(project)
    marker_plot_list <- plotEmbedding(project,
        colorBy = "GeneScoreMatrix",
        name = marker_tb$gene_symbol,
        embedding = "UMAP",
        quantCut = c(0.01, 0.95),
        imputeWeights = NULL
    )
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
    marker_ncol <- marker_nrow <- ceiling(sqrt(length(marker_plot_fmt)))
    png(file.path(out_archr_project, "Plots", "markerGene_UMAP.png"), width = 2 * marker_ncol, height = 2 * marker_nrow, units = "in", res = 100)
    wrap_plots(marker_plot_fmt, ncol = marker_ncol, nrow = marker_nrow)
    dev.off()
    
    # metadata UMAP
    meta_plot_list <- plotEmbedding(project,
        colorBy = "cellColData",
        name = c("Clusters", "FinalSite", "Clinical.Dx", "region", "PrepBatch", "SeqBatch"),
        embedding = "UMAP",
        quantCut = c(0.01, 0.95),
        imputeWeights = NULL
    )
    dev.off()
    meta_plot_fmt <- map(meta_plot_list, function(x) {
        x +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    })
    meta_plot_ncol <- meta_plot_nrow <- ceiling(sqrt(length(meta_plot_fmt)))
    png(file.path(out_archr_project, "Plots", "meta_plot_UMAP.png"), width = 6 * meta_plot_ncol, height = 6 * meta_plot_nrow, units = "in", res = 100)
    wrap_plots(meta_plot_fmt, ncol = meta_plot_ncol, nrow = meta_plot_nrow)
    dev.off()

    # individual library id UMAP
    umap_x <- getEmbedding(project)[[1]]
    umap_y <- getEmbedding(project)[[2]]
    library_id_plotlist <- map(unique(project$Sample), function(s) {
        writeLines(s)
        cell_color <- as.numeric(project$Sample == s)
        print(table(cell_color))
        ggPoint(umap_x, umap_y, cell_color, highlightPoints = (which(project$Sample == s)), 
            discrete = TRUE, labelMeans = FALSE, labelAsFactors = TRUE, title = s, alpha = 0.8, keepAxis = FALSE) + 
            guides(color = FALSE, fill = FALSE) +
            theme_ArchR(baseSize = 6.5) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
            theme(
                axis.text.x=element_blank(), 
                axis.ticks.x=element_blank(), 
                axis.text.y=element_blank(), 
                axis.ticks.y=element_blank()
            )
    })
    dev.off()
    plotlist_png_square(library_id_plotlist, file.path(out_archr_project, "Plots", "library_id_UMAP.png"), width = 4, height = 4)

    dx_plotlist <- plot_individual_meta_lvls(project, "Clinical.Dx", "ggHex")
    plotlist_png_square(dx_plotlist, file.path(out_archr_project, "Plots", "clinical_dx.png"), width = 6)

    region_plotlist <- plot_individual_meta_lvls(project, "region", "ggHex")
    plotlist_png_square(region_plotlist, file.path(out_archr_project, "Plots", "region.png"), width = 6)

    saveArchRProject(project, out_archr_project, load = FALSE)
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
    png(path, width = width * plot_ncol, height = height * plot_nrow, units = units, res = res)
    print(wrap_plots(plotlist, ncol = plot_ncol, nrow = plot_nrow))
    dev.off()
}

if (!interactive()) {
    main()
}
