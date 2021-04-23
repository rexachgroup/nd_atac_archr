# Subclustering of harmony clusters.
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, "clustering_harmony")
celltype_markers <- normalizePath("../ext/20210128_cell_markers_noependymial.csv")
seurat_obj <- normalizePath("../../nucseq_nd_dpolioud/analysis/pci_import/pci_seurat.rds")
out_dir <- file.path(base_dir, "07_archr_harmony_subclustering")
batchtools <- file.path(out_dir, "batchtools")
drop_samples <- c("P1_7_at1_7", "i3_6_at", "I1_7")

RESOURCES <- list(
    ncpus = 8,
    memory = 12 * 8,
    walltime = 86400
)

setBold <- system("tput rev", intern = TRUE)
setNorm <- system("tput sgr0", intern = TRUE)

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(setBold, msg, setNorm), side = "both"))
}

main <- function() {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(batchtools)
    } else {
        reg <- loadRegistry(batchtools, writeable = TRUE)
    }

    reg$packages <- liblist
    batchExport(mget(ls()))

    cluster_args_tb <- tribble(
        ~region,        ~clusters,                  ~proj_name,
        "PreCG",        c("C2", "C5"),              "precg-c25",
        "midInsula",    c("C2", "C5"),              "insula-c25",
        "PreCG",        c("C6", "C7", "C8", "C9"),  "precg-c6789",
        "midInsula",    c("C6", "C7", "C8", "C9"),  "insula-c6789"
    )
    cluster_args_tb <- mutate(cluster_args_tb, out_path = file.path(out_dir, paste("clustering", proj_name, sep = "_")))
    clearRegistry()
    ids <- batchMap(subcluster_worker, 
        args = list(
            f_region = cluster_args_tb$region,
            f_cluster = cluster_args_tb$clusters,
            out_path = cluster_args_tb$out_path
        ),
        more.args = list(
            proj_dir = archr_project,
            drop_samples = drop_samples
        )
    )
    submitJobs(ids, resources = RESOURCES)
    waitForJobs()
    saveRDS(cluster_args_tb, file.path(out_dir, "cluster_args_tb.rds"))
}

subcluster_worker <- function(proj_dir, f_region, f_cluster, drop_samples, out_path) {
    addArchRThreads(RESOURCES$ncpus)
    addArchRGenome("hg38")
    project <- loadArchRProject(proj_dir)
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")
    filtered_cells <- cell_data %>%
        filter(!Sample %in% drop_samples) %>%
        filter(region %in% f_region & Clusters %in% f_cluster)
    writeMsg(str_glue("{length(filtered_cells$cell_id)} / {nrow(cell_data)} cells after filtering"))
    writeMsg(str_glue("regions: {paste0(unique(filtered_cells$region), collapse = ' ')}"))
    writeMsg(str_glue("clusters: {paste0(unique(filtered_cells$Clusters), collapse = ' ')}"))

    writeMsg("subset")
    project <- subsetArchRProject(project, cells = filtered_cells$cell_id, outputDirectory = out_path, force = TRUE, dropCells = FALSE)
    setwd(out_path)
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")

    # LSI
    project <- addIterativeLSI(project,
        useMatrix = "TileMatrix",
        name = "IterativeLSI",
        iterations = 6,
        clusterParams = list(
            resolution = c(0.1, 0.2, 0.4, 0.8, 2),
            sampleCells = 10000,
            n.start = 10
        ),
        varFeatures = 100000,
        dimsToUse = 1:30,
        force = TRUE
    )
    
    print(project@reducedDims)

    # Clusters
    project <- addClusters(project,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        resolution = 0.8,
        force = TRUE
    )
    
    # Drop really small clusters.
    cluster_cts <- as.list(table(project$Clusters))
    small_clusters <- names(which(cluster_cts < 1000))
    project <- project[!project$Clusters %in% small_clusters, ]
    saveArchRProject(project, out_path, load = FALSE)

    # cell UMAP
    project <- addUMAP(project,
        reducedDims = "IterativeLSI",
        name = "UMAP",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine",
        force = TRUE
    )
 
    # marker genes
    project <- addImputeWeights(project)
    marker_genescores <- getMarkerFeatures(project,
        useMatrix = "GeneScoreMatrix",
        groupBy = "Clusters",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )
    saveRDS(marker_genescores, file.path(out_path, "marker_genescores.rds"))
    gene_names <- rowData(marker_genescores)$name
    marker_assay_matrices <- imap(as.list(assays(marker_genescores)), function(x, n) {
            assay_df <- as.data.frame(x)
            colnames(assay_df) <- paste(colnames(marker_genescores), n, sep = "_")
            return(assay_df)
        }) %>% bind_cols %>%
        mutate(gene = gene_names) %>%
        select(gene, everything())
    write_csv(marker_assay_matrices, file.path(cr$out_path, "marker_genescores.csv"))

    saveArchRProject(project, out_path, load = FALSE)
}

if (!interactive()) {
    main()
}
