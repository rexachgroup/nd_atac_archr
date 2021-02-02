liblist <- c("tidyverse", "Seurat", "Signac", "readxl", "GenomeInfoDb", "EnsDb.Hsapiens.v86")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

source("00_CONSTANTS.R")
round2_10x_aggr <- file.path(CELLRANGER_ATAC_AGGR_DIR, "precg-atac-2020-subsample")
round1_10x_aggr <- "/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/ATAC_aggr"
metadata <- "/geschwindlabshares/lchenprj01/nucseq_combined/data/snATAC_metadata_summary_2021_c.xlsx"
celltype_markers <- "../ext/20210128_cell_markers_noependymial.csv"
out_dir <- "../data/signac/round2_10x_aggr/"

data_dir <- file.path(out_dir, "data")
plot_dir <- file.path(out_dir, "plot")

main <- function() {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    meta <- read_xlsx(metadata)
    round2_meta <- cellranger_import_meta(round2_10x_aggr)
    round2_meta <- inner_join(round2_meta, meta, by = c("library_id" = "ATAC_sampleID")) %>%
        as.data.frame
    rownames(round2_meta) <- round2_meta$barcode
    round2_meta$cell_ids <- round2_meta$barcode
    round2_sobj <- signac_cellranger_import(round2_10x_aggr, round2_meta)
    print(dim(t(round2_meta)))
    print(dim(round2_sobj))
    saveRDS(round2_sobj, file.path(data_dir, "import.rds"), compress = FALSE)
    
    round2_sobj <- signac_annot_genes(round2_sobj)
    round2_sobj <- signac_annot_nucleosome(round2_sobj)
    round2_sobj <- signac_annot_tss(round2_sobj)

    pdf(file.path(plot_dir, "annot_plots.pdf"), height = 7, width = 14)
    tmp <- round2_sobj
    tmp$high.tss <- ifelse(tmp$TSS.enrichment > 2, 'High (TSS.enrichment > 2)', 'Low (TSS.Enrichment < 2)')
    tmp$nucleosome_group <- ifelse(tmp$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    TSSPlot(tmp, group.by = 'high.tss') + NoLegend()
    FragmentHistogram(object = tmp, group.by = 'nucleosome_group')
    dev.off()
    
    png(file.path(plot_dir, "qc_plots.png"), width = 20, height = 10, units = "in", res = 400)
    VlnPlot(
        object = round2_sobj,
        features = c('pct_reads_in_peaks', 'peak_region_fragments',
                   'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
        pt.size = 0.1,
        ncol = 5
    )
    dev.off()

    saveRDS(round2_sobj, file.path(data_dir, "prefilter.rds"))
    round2_sobj <- signac_preprocess_filter(round2_sobj)
    print(dim(round2_sobj))
    saveRDS(round2_sobj, file.path(data_dir, "preclustering.rds"))
    
    round2_sobj <- signac_nlme(round2_sobj)
    saveRDS(round2_sobj, file.path(data_dir, "postclustering.rds"))
    pdf(file.path(plot_dir, "clustering_umap.pdf"))
    DimPlot(round2_sobj, label = TRUE) + NoLegend()
    DimPlot(round2_sobj, label = FALSE, group.by = "library_id")
    dev.off()
    round2_sobj <- signac_gene_activity(round2_sobj)
    saveRDS(round2_sobj, file.path(data_dir, "postgeneact.rds"))

    pdf(file.path(plot_dir, "marker_umap.pdf"), width = 20, height = 20)
    DefaultAssay(round2_sobj) <- 'RNA'
    FeaturePlot(
        object = round2_sobj,
        features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
        pt.size = 0.1,
        max.cutoff = 'q95',
        ncol = 3
    )
    dev.off()

    ct_markers <- read_csv(celltype_markers)
    round2_sobj <- assign_celltype(round2_sobj, ct_markers)
    round2_sobj$cluster_ids <- round2_sobj[["peaks_snn_res.0.8"]]
    round2_sobj <- assign_cluster_celltype(round2_sobj)
    #     saveRDS(round2_sobj, file.path(data_dir, "postcelltype.rds"))
    saveRDS(round2_sobj, file.path(data_dir, "signac.rds"), compress = FALSE)

}


cellranger_import_meta <- function(cellranger_dir) {
    cellranger_meta <- read_csv(file.path(cellranger_dir, "outs", "singlecell.csv"))
    aggregation_meta <- read_csv(file.path(cellranger_dir, "outs", "aggregation_csv.csv"))

    cellranger_meta <- cellranger_meta %>%
        mutate(aggregation_no = str_extract(barcode, "\\d+$"))

    aggregation_meta <- aggregation_meta %>%
        mutate(aggregation_no = as.character(1:nrow(.)))

    cell_meta <- inner_join(cellranger_meta, aggregation_meta, by = "aggregation_no") 
    return(cell_meta)
}

signac_cellranger_import <- function(cellranger_dir, meta, downsample = NULL) {
    peak_path <- file.path(cellranger_dir, "outs", "filtered_peak_bc_matrix.h5")
    frag_path <- file.path(cellranger_dir, "outs", "fragments.tsv.gz")
    writeLines(str_glue("signac_cellranger_import: \n peak file:{peak_path} \nfrag file:{frag_path}"))
    counts <- Read10X_h5(peak_path)
    stopifnot(all(colnames(counts) %in% rownames(meta)))
    meta <- meta[colnames(counts), ]

    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        genome = "hg19",
        fragments = frag_path,
        min.cells = 1
    )
    sobj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = meta)

    if (!is.null(downsample)) {
        writeLines(str_glue("signac_cellranger_import: dowsampling to {downsample} cells"))
        cell_ids <- sample(colnames(sobj), downsample)
        sobj <- subset(sobj, cells = cell_ids)
    }
    return(sobj)
}

signac_annot_genes <- function(sobj) {
    writeLines("signac_annot_genes: add gene information")
    annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86, biotypes = c("protein_coding"))
    seqlevelsStyle(annot) <- "UCSC"
    genome(annot) <- "hg19"
    Annotation(sobj) <- annot
    return(sobj)
}

signac_annot_nucleosome <- function(sobj) {
    writeLines("signac_annot_nucleosome: NucleosomeSignal")
    sobj <- NucleosomeSignal(sobj)
    sobj$pct_reads_in_peaks <- sobj$peak_region_fragments / sobj$passed_filters * 100
    sobj$blacklist_ratio <- sobj$blacklist_region_fragments / sobj$peak_region_fragments
    return(sobj) 
}

signac_annot_tss <- function(sobj) {
    writeLines("signac_annot_tss: TSSEnrichment")
    return(TSSEnrichment(sobj, fast = FALSE))
}

signac_preprocess_filter <- function(sobj) {
    writeLines("signac_preprocess_filter")
    return(subset(sobj, subset = 
        peak_region_fragments > 1000 &
        peak_region_fragments < 20000 &
        pct_reads_in_peaks > 15 &
        blacklist_ratio < 0.05 &
        nucleosome_signal < 10 &
        TSS.enrichment > 2
    )) 
}

signac_nlme <- function(sobj, dims = 1:30) {
    writeLines("signac_nlme: normalization, linear dimensional reduction")
    sobj <- sobj %>%
        RunTFIDF %>%
        FindTopFeatures(min.cutoff = 'q0') %>%
        RunSVD(
            assay = "peaks",
            reduction.key = "LSI_",
            reduction.name = "lsi"
        )
    writeLines(str_glue("signac_nlme: non-linear dimensional reduction and clustering, dims = {min(dims)}:{max(dims)}"))
    sobj <- sobj %>%
        RunUMAP(reduction = "lsi", dims = dims) %>%
        FindNeighbors(reduction = "lsi", dims = dims) %>%
        FindClusters(verbose = FALSE)
    return(sobj)
}

signac_gene_activity <- function(sobj) { 
    writeLines("signac_gene_activity")
    gact <- GeneActivity(sobj, extend.upstream = 2000)
    sobj[["RNA"]] <- CreateAssayObject(counts = gact)
    sobj <- NormalizeData(sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = median(sobj$nCount_RNA))
    return(sobj)
}


assign_celltype <- function(sobj, markers) {
    writeLines("assign_celltype")
    markers <- markers %>%
        group_by(marker_for)
    ct_names <- group_keys(markers)$marker_for
    marker_split <- markers %>%
        group_split() %>%
        map(pluck("gene_symbol"))
    sobj@meta.data <- sobj@meta.data %>%
        dplyr::select(!contains(ct_names))
    sobj <- AddModuleScore(
        sobj,
        features = marker_split,
        assay = "RNA",
        name = ct_names
    )

    # Rename columns by dropping trailing digit.
    sobj@meta.data <- sobj@meta.data %>%
        rename_with(.cols = contains(ct_names), ~str_extract(., "^\\D+"))

    max_ct_score <- sobj@meta.data %>%
        as_tibble %>%
        dplyr::select(cell_ids, contains(ct_names)) %>%
        pivot_longer(cols = contains(ct_names), names_to = "cell_type", values_to = "celltype_score") %>%
        group_by(cell_ids) %>%
        dplyr::slice_max(celltype_score, with_ties = FALSE) %>%
        dplyr::select(cell_ids, cell_type)

    sobj@meta.data <- inner_join(sobj@meta.data, max_ct_score, by = "cell_ids")
    return(sobj)
}

assign_cluster_celltype <- function(sobj, cluster_ct_cutoff = 0.2) {
    writeLines("assign_cluster_celltype")

    cluster_counts <- sobj@meta.data %>%
        as_tibble %>%
        group_by(cluster_ids) %>%
        add_count(name = "cluster_total") %>%
        group_by(cluster_ids, cluster_total, cell_type) %>%
        count(name = "ct_total")

    cluster_assign <- cluster_counts %>%
        mutate(ct_frac = ct_total / cluster_total) %>%
        group_by(cluster_ids) %>%
        slice_max(ct_frac, with_ties = FALSE) %>%
        mutate(cluster_cell_type = ifelse(ct_frac > cluster_ct_cutoff, cell_type, "unknown")) %>%
        dplyr::select(cluster_ids, cluster_cell_type)

    sobj@meta.data <- inner_join(sobj@meta.data, cluster_assign, by = "cluster_ids")
    return(sobj)
}

if (!interactive()) {
    main()
}
