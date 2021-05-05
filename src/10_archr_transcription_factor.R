# TF / ChromVAR motif derivation

liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "circlize", "scales", "RColorBrewer", "ComplexHeatmap")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, c("peak_calling_precg_C2", "peak_calling_precg_C7", "peak_calling_insula_C2", "peak_calling_insula_C7"))
out_archr_project <- file.path(base_dir, c("chrom_var_precg_C2", "chrom_var_precg_C7", "chrom_var_insula_C2", "chrom_var_insula_C7"))
batchtools <- file.path(base_dir, "chromvar_batchtools")
setBold <- system("tput rev", intern = TRUE)
setNorm <- system("tput sgr0", intern = TRUE)

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(setBold, msg, setNorm), side = "both"))
}
alert <- function() { system("tput bel") }
RESOURCES <- list(
    ncpus = 8,
    memory = 4 * 8,
    walltime = 86400
)

main <- function() {
    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(batchtools)
    } else {
        reg <- loadRegistry(batchtools, writeable = TRUE)
    }
    clearRegistry()
    reg$packages <- liblist
    batchExport(mget(ls()))
    ids <- batchMap(tf_worker, args = list(proj_dir = archr_project, out_dir = out_archr_project))
    setJobNames(ids, basename(out_archr_project))
    submitJobs(ids, resources = RESOURCES)
    waitForJobs()

    dx_marker_features <- map(out_archr_project, tf_dx_marker_signif)
    pmap(list(out_archr_project, dx_marker_features), tf_dx_marker_signif_export)
    pmap(list(out_archr_project, dx_marker_features), chromvar_plot_worker)
    map(out_archr_project, chromvar_export)
}

tf_worker <- function(proj_dir, out_dir) {
    startTime <- Sys.time()
    addArchRGenome("hg38")
    addArchRThreads(8)

    writeMsg(str_glue("copy {proj_dir} {out_dir}"))
    project <- loadArchRProject(path = proj_dir)
    project <- saveArchRProject(project, out_dir, load = TRUE)
    plot_dir <- file.path(out_dir, "Plots")
    orig_dir <- getwd()
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(out_dir)

    writeMsg("TF motif annotations")
    project <- addMotifAnnotations(project, motifSet =  "JASPAR2018", species = "Homo sapiens", name = "Motif", force = TRUE)
    writeMsg("background peaks")
    project <- addBgdPeaks(project)
    writeMsg("derivations / chromVAR matrix")
    project <- addDeviationsMatrix(project, peakAnnotation = "Motif", force = TRUE)

    saveArchRProject(project, out_dir, load = FALSE)

    endTime <- Sys.time()
    writeMsg(str_glue("{format(endTime - startTime)}"))
}

tf_dx_marker_signif <- function(proj_dir) {
    project <- loadArchRProject(proj_dir)
    writeMsg(str_glue("get motif matrix {proj_dir}"))
    dx_wk <- function(dx) {
        writeMsg(str_glue("per_dx_marker_features {dx}"))
        peakMarkers <- getMarkerFeatures(
            project,
            useMatrix = "MotifMatrix",
            useSeqnames = "z",
            useGroups = dx,
            bgdGroups = "Control",
            groupBy = "Clinical.Dx",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = "wilcoxon"
        )
    }
    dx_marker_features <- map(c("AD", "bvFTD", "PSP-S"), ~tryCatch(dx_wk(.), error = print))
    setNames(dx_marker_features, c("AD", "bvFTD", "PSP-S"))
}

tf_dx_marker_signif_export <- function(proj_dir, dx_marker_signif_list) {
    writeMsg(str_glue("export motif dx summarizedExperiments to {proj_dir}/motifMatrix_dx.rds"))
    saveRDS(dx_marker_signif_list, file.path(proj_dir, "motifMatrix_dx.rds"))
    out_dir <- file.path(proj_dir, "TF_dx")
    dir.create(out_dir, showWarnings = FALSE)
    dx_wk <- function(marker_signif) {
        rd <- rowData(marker_signif)
        marker_assays <- setNames(Reduce(cbind, assays(marker_signif)), nm = names(assays(marker_signif)))
        marker_tb <- as.data.frame(cbind(rd, marker_assays))
        glimpse(marker_tb)
        out_path <- file.path(out_dir, str_glue("{basename(proj_dir)}_{colnames(marker_signif)}.csv"))
        writeMsg(str_glue("export motif dx test to {out_path}"))
        write_csv(marker_tb, out_path)
    }
    walk(dx_marker_signif_list, dx_wk)
}

chromvar_export <- function(proj_dir) {
    project <- loadArchRProject(proj_dir)
    motifMatrix <- getMatrixFromProject(project, useMatrix = "MotifMatrix")
    saveRDS(motifMatrix, file.path(proj_dir, "motifMatrix.rds"))
}

chromvar_plot_worker <- function(proj_dir, dx_marker_signif_list) {
    project <- loadArchRProject(path = proj_dir)
    proj_name <- basename(proj_dir)
    writeMsg(str_glue("get motif matrix {proj_dir}"))
    motifMatrix <- getMatrixFromProject(project, useMatrix = "MotifMatrix")

    deviation <- assays(motifMatrix)$deviations
    zscore <- assays(motifMatrix)$z
    anno <- colData(motifMatrix)

    writeMsg("deviations plot")
    var_dev_plot <- getVarDeviations(project, name = "MotifMatrix", plot = TRUE)
    pdf(file.path(proj_dir, "Plots", str_glue("TFvar_{proj_name}.pdf")))
    print(var_dev_plot)
    graphics.off() 


    writeMsg("motif enrichment significance heatmaps")
    dx_motif_heatmaps <- map(dx_marker_signif_list, function(marker_features) {
        plot_name <- colnames(marker_features)[[1]]
        marker_pval <- assay(marker_features, "Pval")
        motif_order <- order(marker_pval)[1:50]
        motif_label <- rowData(marker_features)$name[motif_order]
        
        dx_palette <- setNames(brewer.pal(length(unique(anno$Clinical.Dx)), name = "GnBu"), nm = unique(anno$Clinical.Dx))


        dx_col_annotation <- HeatmapAnnotation(dx = anno$Clinical.Dx, col = list(dx = dx_palette))
        motif_row_annotation <- HeatmapAnnotation(link = anno_mark(at = motif_order, labels = motif_label), which = "row")

        zscore_color <- colorRamp2(
            breaks = c(min(zscore), 0, max(zscore)), 
            colors = c("blue", "white", "red")
        )

        hmap <- Heatmap(
            mat = as.matrix(zscore),
            col = zscore_color,
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            right_annotation = motif_row_annotation,
            bottom_annotation = dx_col_annotation,
            use_raster = TRUE,
            name = plot_name,
            column_title = str_glue(basename(proj_dir), plot_name)
        )

        return(hmap)
    })
    writeMsg("motif enrichment significance volcano")
    pdf(file.path(proj_dir, "Plots", "zscore_dx.pdf"), height = 14, width = 8)
    print(dx_motif_heatmaps)
    graphics.off()
}

if (!interactive()) {
    main()
}
