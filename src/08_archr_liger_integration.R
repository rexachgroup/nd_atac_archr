liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "ragg")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, "clustering_harmony")
out_archr_project <- file.path(base_dir, "liger_integration_insula")
plot_dir <- file.path(out_archr_project, "Plots")

seurat_object <- "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/pci_import/pci_seurat.rds"
liger_metadata <- "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/seurat_lchen/liger_subcluster_metadata.rds"
seurat_object_filter <- normalizePath("../data/subclusters_removed_byQC_final.xlsx")

main <- function() { 
    addArchRGenome("hg38")
    addArchRThreads(16)

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
    
    writeLines("subset to insula")
    sobj_filtered_cells <- filter(sobj_meta, !ct_subcluster %in% sobj_badclusters$ct_subcluster) %>%
        filter(region == "insula")
    writeLines(str_glue("sobj_filtered_cells: {nrow(sobj_filtered_cells)} / {nrow(sobj_meta)} {round(nrow(sobj_filtered_cells) / nrow(sobj_meta), 2)}"))
    sobj <- subset(sobj, cells = sobj_filtered_cells$cell_ids)
    sobj$ct_subcluster <- sobj_filtered_cells$ct_subcluster
    gc()

    archr_meta <- project@cellColData %>% as_tibble(rownames = "cell_ids")
    archr_filtered_meta <- filter(archr_meta, region == "midInsula")
    writeLines(str_glue("archr_filtered_meta: {nrow(archr_filtered_meta)} / {nrow(archr_meta)} {round(nrow(archr_filtered_meta) / nrow(archr_meta), 2)}"))
    project <- project[archr_filtered_meta$cell_ids, ]
    project <- saveArchRProject(project, out_archr_project, dropCells = TRUE, load = TRUE)
    
    nATAC <- 40000
    nRNA <- 25000

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
    saveArchRProject(project, out_archr_project, load = FALSE)
    
    clustering_col_ct_annotation <- sobj_filtered_cells %>% select(ct_subcluster, cluster_cell_type) %>%
        filter(!duplicated(ct_subcluster), ct_subcluster %in% project$predictedGroup_Un) %>%
        arrange(cluster_cell_type, ct_subcluster) %>%
        as.data.frame %>%
        column_to_rownames("ct_subcluster")

    unconstrained_clust <- as.matrix(confusionMatrix(project$Clusters, project$predictedGroup_Un))
    write_csv(as.data.frame(unconstrained_clust), file.path(out_archr_project, "unconstrained_clust.csv"))
    # sort order for clustering: by cluster_cell_type for seurat clusters and by natural ordering for atac clusters
    unconstrained_clust <- unconstrained_clust[str_sort(unique(project$Clusters), numeric = TRUE), as.character(rownames(clustering_col_ct_annotation))]
    unconstrained_clust <- log10(unconstrained_clust + 1)
    clustering_heatmap <- pheatmap(
        unconstrained_clust, 
        main = str_glue("nATAC = {nATAC}, nRNA = {nRNA}, mid-insula cells -- log10(n + 1) of overlap"),
        border_color = "black", 
        annotation_col = clustering_col_ct_annotation, 
        annotation_colors = list(cluster_cell_type = paletteDiscrete(clustering_col_ct_annotation$cluster_cell_type)),
        cluster_rows = TRUE, 
        cluster_cols = FALSE
    )
    pdf(file.path(out_archr_project, "Plots", "gene_integration_matrix.pdf"), width = 0.30 * ncol(unconstrained_clust), height = 0.30 * nrow(unconstrained_clust))
    print(clustering_heatmap)
    graphics.off()

}
