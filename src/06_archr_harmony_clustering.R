liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

SAMPLE_META <- normalizePath("../data/snATAC_metadata_summary_2021_e.xlsx")
base_dir <- normalizePath("../data/archr/atac-2020-all")
data_dir <- file.path(base_dir, "data")
plot_dir <- file.path(base_dir, "plot")
archr_project <- file.path(base_dir, "preprocess")
sample_table <- file.path(base_dir, "preprocess", "sample_tb.rds")
out_archr_project <- file.path(base_dir, "clustering_harmony")
celltype_markers <- normalizePath("../ext/20210128_cell_markers_noependymial.csv")
seurat_obj <- normalizePath("../../nucseq_nd_dpolioud/analysis/pci_import/pci_seurat.rds")
drop_samples <- c("P1_7_at1_7", "i3_6_at", "I1_7")

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(setBold, msg, setNorm), side = "both"))
}

main <- function() {
    addArchRThreads()
    addArchRGenome("hg38")
    writeLines(str_glue("load {archr_project}"))
    project <- loadArchRProject(path = archr_project)
    dir.create(file.path(out_archr_project, "Plots"), showWarning = FALSE, recursive = TRUE) 
    orig_dir <- getwd()
    setwd(out_archr_project)
    marker_tb <- read_csv(celltype_markers)

    writeMsg("copy metadata from sample_meta to cellColData")
    sample_meta <- read_xlsx(SAMPLE_META)
    cell_meta <- as_tibble(project@cellColData, rownames = "cell_id")
    meta_join <- left_join(cell_meta, sample_meta, by = c("Sample" = "ATAC_fastq_name"))
    stopifnot(rownames(project) == meta_join$cell_id)
    new_cols <- setdiff(colnames(meta_join), colnames(cell_meta))
    project@cellColData[, new_cols] <- meta_join[, new_cols]

    writeMsg("pre-cluster doublet removal")
    filterDoublets(project)
    
    # additional filtering on top of import filtering-- TSSEnrichment > 4, BlacklistRatio < 0.1
    # Drop samples P1_7at1_7. i3_6_at, I1_7
    writeMsg("pre-clustering cell / sample filtering")
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")
    filtered_cells <- cell_data %>%
        filter(TSSEnrichment > 4, BlacklistRatio < 0.1) %>%
        filter(!Sample %in% drop_samples) %>%
        pluck("cell_id")
    writeLines(str_glue("{length(filtered_cells)} / {nrow(cell_data)} cells after filtering"))

    writeLines("subset")
    project <- project[filtered_cells, ]
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")
    
    # Copy project.
    project <- saveArchRProject(project, out_archr_project, dropCells = TRUE)

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
    project <- addHarmony(project,
        groupBy = c("PrepBatch", "SeqBatch"),
        force = TRUE
    )
    saveArchRProject(project, out_archr_project, load = FALSE)
    print(project@reducedDims)

    # Clusters
    project <- addClusters(project,
        reducedDims = "Harmony",
        method = "Seurat",
        name = "Clusters",
        resolution = 0.8,
        force = TRUE
    )
    saveArchRProject(project, out_archr_project, load = FALSE)
    
    # Drop really small clusters.
    cluster_cts <- as.list(table(project$Clusters))
    small_clusters <- names(which(cluster_cts < 1000))
    project <- project[!project$Clusters %in% small_clusters, ]
    saveArchRProject(project, out_archr_project, load = FALSE)

    # cell UMAP
    project <- addUMAP(project,
        reducedDims = "Harmony",
        name = "UMAP",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine"
    )
    saveArchRProject(project, out_archr_project, load = FALSE)

 
    # marker genes
    project <- addImputeWeights(project)
    marker_genescores <- getMarkerFeatures(project,
        useMatrix = "GeneScoreMatrix",
        groupBy = "Clusters",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )
    saveRDS(marker_genescores, file.path(out_archr_project, "marker_genescores.rds"))

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
