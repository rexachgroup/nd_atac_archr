liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "ragg")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, "clustering_harmony")
out_archr_project <- file.path(base_dir, "liger_integration_precg")
plot_dir <- file.path(out_archr_project, "Plots")

seurat_object <- "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/pci_import/pci_seurat.rds"
liger_metadata <- "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/seurat_lchen/liger_subcluster_metadata.rds"
seurat_object_filter <- normalizePath("../data/subclusters_removed_byQC_final.xlsx")

main <- function() { 
    addArchRGenome("hg38")
    addArchRThreads(4)

    writeLines(str_glue("load {seurat_object}"))
    sobj <- readRDS(seurat_object)

    writeLines(str_glue("load {seurat_object_filter}"))
    sobj_badclusters <- read_xlsx(seurat_object_filter)

    writeLines(str_glue("load {archr_project}"))
    project <- loadArchRProject(path = archr_project)
    orig_dir <- getwd()
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(out_archr_project)
    
    writeLines(str_glue("load {liger_metadata}"))
    sobj_meta <- readRDS(liger_metadata) %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
        )
    
    writeLines("subset to precg")
    sobj_filtered_cells <- filter(sobj_meta, !ct_subcluster %in% sobj_badclusters$ct_subcluster) %>%
        filter(region == "preCG")
    writeLines(str_glue("sobj_filtered_cells: {nrow(sobj_filtered_cells)} / {nrow(sobj_meta)} {round(nrow(sobj_filtered_cells) / nrow(sobj_meta), 2)}"))
    sobj <- subset(sobj, cells = sobj_filtered_cells$cell_ids)
    sobj$ct_subcluster <- sobj_filtered_cells$ct_subcluster
    gc()

    archr_meta <- project@cellColData %>% as_tibble(rownames = "cell_ids")
    archr_filtered_meta <- filter(archr_meta, region == "PreCG")
    writeLines(str_glue("archr_filtered_meta: {nrow(archr_filtered_meta)} / {nrow(archr_meta)} {round(nrow(archr_filtered_meta) / nrow(archr_meta), 2)}"))
    project <- project[archr_filtered_meta$cell_ids, ]
    project <- saveArchRProject(project, out_archr_project, dropCells = TRUE, load = TRUE)
    
    nATAC <- 40000
    nRNA <- 40000

    project <- addGeneIntegrationMatrix(project,
        useMatrix = "GeneScoreMatrix",
        matrixName = "GeneIntegrationMatrix",
        reducedDims = "Harmony",
        seRNA = sobj,
        sampleCellsATAC = nATAC,
        sampleCellsRNA = nRNA,
        nGenes = 4000,
        addToArrow = FALSE,
        groupRNA = "ct_subcluster",
        nameCell = "predictedCell_Un",
        nameGroup = "predictedGroup_Un",
        nameScore = "predictedScore_Un"
    )
    project <- saveArchRProject(project, out_archr_project, load = TRUE)

    colorbar_tb <- colorbar_ct_tb(sobj_filtered_cells)
    clustering_mat <- clustering_data_mat(project, colorbar_tb)
    
    write_csv(as.data.frame(clustering_mat), file.path(out_archr_project, "clustering_data_mat.csv"))

    clustering_log_mat <- log10(clustering_mat + 1)
    clustering_heatmap <- pheatmap(
        clustering_log_mat, 
        main = str_glue("nATAC = {nATAC}, nRNA = {nRNA}, precg cells -- log10(n + 1) of overlap"),
        border_color = "black", 
        annotation_col = colorbar_tb, 
        annotation_colors = list(cluster_cell_type = paletteDiscrete(colorbar_tb$cluster_cell_type)),
        cluster_rows = TRUE, 
        cluster_cols = FALSE
    )

    pdf(file.path(out_archr_project, "Plots", "gene_integration_matrix.pdf"), width = 0.30 * ncol(clustering_log_mat), height = 0.30 * nrow(clustering_log_mat))
    print(clustering_heatmap)
    graphics.off()

}

colorbar_ct_tb <- function(seurat_meta) {
    seurat_meta %>%
        select(ct_subcluster, cluster_cell_type) %>%
        group_by(ct_subcluster) %>%
        slice_head(n = 1) %>%
        arrange(cluster_cell_type, ct_subcluster) %>%
        as.data.frame %>%
        column_to_rownames("ct_subcluster")
}

clustering_data_mat <- function(project, colorbar_ct_tb) {
    mat <- as.matrix(confusionMatrix(project$Clusters, project$predictedGroup_Un))
    clusters <- unique(project$Clusters)
    missing_predicted <- setdiff(rownames(colorbar_ct_tb), project$predictedGroup_Un)
    # add clusters that aren't present in predictedGroup_Un
    writeLines(str_glue("Missing subcluster in prediction: {missing_predicted}. Setting subcluster cols to 0."))
    blankmat <- matrix(
        0, 
        nrow = length(clusters), 
        ncol = length(missing_predicted),
        dimnames = list(clusters, missing_predicted)
    )
    mat <- cbind(mat, blankmat)
    # sort order for clustering: by natural ordering for atac clusters and cluster_cell_type for seurat clusters
    mat <- mat[str_sort(clusters, numeric = TRUE), as.character(rownames(colorbar_ct_tb))]
    return(mat)
}

if (!interactive()) {
    main()
}
