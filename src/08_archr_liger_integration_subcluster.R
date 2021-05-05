# Integration of the subclustered clusters with RNA data.
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "ragg")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(deparse.max.lines = 5)

liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR")
base_dir <- normalizePath("../data/archr/atac-2020-all")
cluster_args <- file.path(base_dir, "07_archr_harmony_subclustering", "cluster_args_tb.rds")
out_dir <- file.path(base_dir, "08_archr_liger_integration_subcluster")
batchtools <- file.path(out_dir, "batchtools")

seurat_object <- "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/pci_import/pci_seurat.rds"
liger_metadata <- "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/seurat_lchen/liger_subcluster_metadata.rds"
seurat_object_filter <- normalizePath("../data/subclusters_removed_byQC_final.xlsx")

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(" ", msg, " "), side = "both", pad = "|"))
}

RESOURCES <- list(
    ncpus = 16,
    memory = 128,
    walltime = 86400
)

nATAC <- 40000
nRNA <- 25000

main <- function() { 
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(batchtools)
    } else {
        reg <- loadRegistry(batchtools, writeable = TRUE)
    }
    
    cluster_args_tb <- readRDS(cluster_args) %>%
        rename("out_path" = "proj_path") %>%
        mutate(out_path = file.path(out_dir, paste("integration", proj_name, sep = "_")))
    
    writeMsg(str_glue("load {seurat_object_filter}"))
    sobj_badclusters <- read_xlsx(seurat_object_filter)

    # load and rename regions to match region names in cluster_args_tb.
    writeMsg(str_glue("load {liger_metadata}"))
    sobj_meta <- readRDS(liger_metadata) %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            region = fct_recode(region, "PreCG" = "preCG", "midInsula" = "insula")
        ) %>%
        filter(!ct_subcluster %in% sobj_badclusters$ct_subcluster)

    # Get sobj cellids for each region.
    cluster_args_tb <- cluster_args_tb %>%
        mutate(
            sobj_cellids = map(region, function(reg) { sobj_meta %>% filter(region == reg) %>% pluck("cell_ids") })
        )

    # batchtools run.
    reg$packages <- liblist
    batchExport(mget(ls()))
    clearRegistry()
    ids <- batchMap(integration_worker, args = select(cluster_args_tb, proj_path, out_path, sobj_cellids))
    ids <- getJobTable() %>% filter(is.na(done))
    submitJobs(ids, RESOURCES)
    waitForJobs()

    saveRDS(cluster_args_tb, file.path(out_dir, "cluster_args_tb.rds"))

    # plotting
    pwalk(cluster_args_tb, function(...) {
            cr <- list(...)
            writeMsg(str_glue("heatmap {cr$proj_name}"))
            project <- loadArchRProject(cr$out_path)
            clust_heatmap(project, cr$out_path, cr$sobj_cellids)
        })

    # copy outputs
    pwalk(cluster_args_tb, function(...) {
            cr <- list(...)
            file.copy(file.path(cr$out_path, "Plots", "gene_integration_matrix.pdf"), file.path(out_path, paste0(cr$proj_name, "-gene_integration_matrix.pdf")), overwrite = TRUE)
            file.copy(file.path(cr$out_path, "RNAIntegration", "unconstrained_clust.csv"), file.path(out_path, paste0(cr$proj_name, "-gene_integration_matrix.csv")), overwrite = TRUE)
        })
}

integration_worker <- function(proj_path, out_path, sobj_cellids) { 
    addArchRGenome("hg38")
    addArchRThreads(16)
    l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

    writeMsg(str_glue("load {seurat_object}"))
    sobj <- readRDS(seurat_object)

    writeMsg(str_glue("load {proj_path}"))
    project <- loadArchRProject(path = proj_path)
    orig_dir <- getwd()
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    setwd(out_path)
    
    writeMsg(str_glue("load {liger_metadata}"))
    sobj_meta <- readRDS(liger_metadata) %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
        )

    sobj_filtered_cells <- filter(sobj_meta, cell_ids %in% sobj_cellids)
    
    writeMsg(str_glue("sobj_filtered_cells: {nrow(sobj_filtered_cells)} / {nrow(sobj_meta)} {round(nrow(sobj_filtered_cells) / nrow(sobj_meta), 2)}"))
    sobj <- subset(sobj, cells = sobj_filtered_cells$cell_ids)
    sobj$ct_subcluster <- sobj_filtered_cells$ct_subcluster
    gc()

    archr_meta <- project@cellColData %>% as_tibble(rownames = "cell_ids")
    project <- saveArchRProject(project, out_path, dropCells = FALSE, load = TRUE)

    writeMsg(str_glue("drop plot dir: {file.path(out_path, 'Plots')}"))
    unlink(file.path(out_path, "Plots"), recursive = TRUE)
    dir.create(file.path(out_path, "Plots"))
    

    writeMsg("integration matrix")
    project <- addGeneIntegrationMatrix(project,
        useMatrix = "GeneScoreMatrix",
        matrixName = "GeneIntegrationMatrix",
        reducedDims = "IterativeLSI",
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
    unconstrained_clust <- as.matrix(confusionMatrix(project$Clusters, project$predictedGroup_Un))
    write_csv(as.data.frame(unconstrained_clust), file.path(out_path, "RNAIntegration", "unconstrained_clust.csv"))
    saveArchRProject(project, out_path, load = FALSE)

    setwd(orig_dir)
}

clust_heatmap <- function(project, out_path, sobj_cellids) {
    writeMsg("gene integration heatmap")
    writeMsg(str_glue("load {liger_metadata}"))
    sobj_meta <- readRDS(liger_metadata) %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
        )

    sobj_filtered_cells <- filter(sobj_meta, cell_ids %in% sobj_cellids)
    
    clustering_col_ct_annotation <- sobj_filtered_cells %>% select(ct_subcluster, cluster_cell_type) %>%
        filter(!duplicated(ct_subcluster), ct_subcluster %in% project$predictedGroup_Un) %>%
        arrange(cluster_cell_type, ct_subcluster) %>%
        as.data.frame %>%
        column_to_rownames("ct_subcluster")
    
    unconstrained_clust <- as.matrix(confusionMatrix(project$Clusters, project$predictedGroup_Un)) 
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
    pdf(file.path(out_path, "Plots", "gene_integration_matrix.pdf"), width = 0.30 * ncol(unconstrained_clust) + 4, height = 0.30 * nrow(unconstrained_clust) + 4)
    print(clustering_heatmap)
    graphics.off()
}

if (!interactive()) {
    main()
}
