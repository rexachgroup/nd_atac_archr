# TF / chromvar motif derivation per subcluster group.

liblist <- c("tidyverse", "batchtools", "TFBSTools","readxl", "Seurat", "ArchR", "circlize", "scales", "RColorBrewer", "ComplexHeatmap")
l <- lapply(liblist, library, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- list(
    "subcluster_peak_call_insula-c25" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_insula-c25")
    #"subcluster_peak_call_precg-c25" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_precg-c25"),
    #"subcluster_peak_call_insula-c6789" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_insula-c6789"),
    #"subcluster_peak_call_precg-c6789" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_precg-c6789")
)
combined_motif_rds <- file.path(base_dir, "00_import_scenic_transfac_db/scenic_transfac_db.rds")
out_dir <- file.path(base_dir, "10_archr_subcluster_transcription_factor")
batchtools_dir <- file.path(out_dir, "batchtools")

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(" ", msg, " "), side = "both", pad = "-"))
}

RESOURCES <- list(
    ncpus = 8,
    memory = 64,
    walltime = 86400,
    partition = "bigmem"
)

main <- function() {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    orig_dir <- getwd()
    setwd(out_dir)
    if (!dir.exists(batchtools_dir)) {
        reg <- makeRegistry(batchtools_dir)
    } else {
        reg <- loadRegistry(batchtools_dir, writeable = TRUE)
    }
    cluster_args_tb <- tibble(archr_project = unlist(archr_project), name = names(archr_project)) %>%
        left_join(., tibble(motif_type = c("jaspar", "scenic")), by = character()) %>%
        mutate(out_proj = file.path(out_dir, str_glue("{name}-{motif_type}"))) %>%
        print()
    combined_motifs <- readRDS(combined_motif_rds)
    
    writeMsg("batchtools tf derivation")
    reg$packages <- liblist
    batchExport(mget(ls()))
    clearRegistry()
    batchMap(
        function(...) {
            cr <- list(...)
            print(cr)
            tf_worker(cr$archr_project, cr$out_proj, cr$motif_type, combined_motifs)
        },
        args = cluster_args_tb
    )
    cluster_args_tb$job.id <- getJobTable()$job.id
    submitJobs(cluster_args_tb$job.id, RESOURCES)
    waitForJobs()
    saveRDS(cluster_args_tb, file.path(out_dir, "tb.rds"))
    writeMsg("done batchtools tf derivation")
    cluster_args_tb <- cluster_args_tb %>% mutate(proj_obj = map(out_proj, loadArchRProject))
    writeMsg("MotifMatrix dx signif")
    cluster_args_tb <- cluster_args_tb %>% mutate(motif_dx_signif = map(proj_obj, tf_signif_motif_dx))
    writeMsg("MotifMatrix plot")
    cluster_args_tb <- cluster_args_tb %>% mutate(chromvar_plots = pmap(., function(...) {
        cr <- list(...)
        chromvar_plots(cr$proj_obj, cr$motif_dx_signif)
    }))
    cluster_args_tb %>%
        pwalk(function(...) {
            cr <- list(...)
            motif_plot <- file.path(out_dir, str_glue("{cr$name}-{cr$motif_type}-vardev.pdf"))
            heat_plot <- file.path(out_dir, str_glue("{cr$name}-{cr$motif_type}-heatmap.pdf"))
            pdf(motif_plot)
            print(cr$chromvar_plots$dev_plot)
            graphics.off()
            pdf(heat_plot, height = 17, width = 10)
            print(cr$chromvar_plots$motif_heat)
            graphics.off()
        })

    cluster_args_tb %>%
        pwalk(function(...) {
            cr <- list(...)
            map(cr$motif_dx_signif, function(marker_signif) {
                rd <- rowData(marker_signif)
                marker_assays <- setNames(Reduce(cbind, assays(marker_signif)), nm = names(assays(marker_signif)))
                marker_tb <- as.data.frame(cbind(rd, marker_assays))
                glimpse(marker_tb)
                out_path <- file.path(out_dir, str_glue("{cr$name}_{colnames(marker_signif)}.csv"))
                writeMsg(str_glue("export motif dx test to {out_path}"))
                write_csv(marker_tb, out_path)
            })
        })
}

devmat_jaspar <- function(project) {
    writeMsg("JASPAR TF motif annotations")
    project <- addMotifAnnotations(project, motifSet =  "JASPAR2018", species = "Homo sapiens", name = "Motif", force = TRUE)
    writeMsg("background peaks")
    project <- addBgdPeaks(project, force = TRUE)
    writeMsg("derivations / chromVAR matrix")
    project <- addDeviationsMatrix(project, peakAnnotation = "Motif", force = TRUE) 
}

devmat_scenic <- function(project, motif_pwmlist) {
    writeMsg("SCENIC / TRANSFAC TF motif annotations")
    project <- addMotifAnnotations(project, motifSet = "custom", motifPWMs = motif_pwmlist, name = "Motif", force = TRUE)
    writeMsg("background peaks")
    project <- addBgdPeaks(project, force = TRUE)
    writeMsg("derivations / chromVAR matrix")
    project <- addDeviationsMatrix(project, peakAnnotation = "Motif", force = TRUE) 
}

tf_worker <- function(proj_dir, out_dir, motif_type, scenic_db) {
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

    if (motif_type == "jaspar") {
        project <- devmat_jaspar(project)
    } else if (motif_type == "scenic"){
        project <- devmat_scenic(project, scenic_db)
    }
    saveArchRProject(project, out_dir, load = FALSE)

    endTime <- Sys.time()
    writeMsg(str_glue("{format(endTime - startTime)}"))
}

tf_signif_motif_dx <- function(project) {
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

chromvar_plots <- function(project, motif_dx_signif) {
    writeMsg(str_glue("get motif matrix"))
    motifMatrix <- getMatrixFromProject(project, useMatrix = "MotifMatrix")

    deviation <- assays(motifMatrix)$deviations
    zscore <- assays(motifMatrix)$z
    anno <- colData(motifMatrix)

    writeMsg("deviations plot")
    var_dev_plot <- getVarDeviations(project, name = "MotifMatrix", plot = TRUE)

    writeMsg("motif enrichment significance heatmaps")
    dx_motif_heatmaps <- map(motif_dx_signif, function(marker_features) {
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
            column_title = str_glue(plot_name)
        )

        return(hmap)
    })
    return(list(dev_plot = var_dev_plot, motif_heat = dx_motif_heatmaps))
}

main()
