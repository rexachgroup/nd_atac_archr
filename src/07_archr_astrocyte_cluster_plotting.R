liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "patchwork", "ragg")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, "clustering_harmony")
out_dir <- file.path(base_dir, "07_archr_astrocyte_cluster_plotting")
out_archr_subset <- file.path(out_dir, "ArchRSubset")

main <- function() {
    dir.create(out_dir)
    orig_dir <- getwd()
    setwd(out_dir)
    if (!file.exists(out_archr_subset)) {
        writeLines(str_glue("load {archr_project}"))
        project <- loadArchRProject(path = archr_project)
            
        cells_f <- project@cellColData %>%
            as.data.frame %>%
            rownames_to_column("barcode") %>%
            filter(Clusters %in% paste("C", 18:25, sep = ""))

        proj_subset <- subsetArchRProject(project, cells = cells_f$barcode, force = T)
    } else {
        proj_subset <- loadArchRProject(out_archr_subset)
    }
    
    
    pdf(file.path(out_dir, "marker_tracks_cluster_astro.pdf"), width = 4, height = 8)
    marker_tracks <- plotBrowserTrack(proj_subset,
        groupBy = "Clusters",
        geneSymbol = c("REST", "CUX2")
    )
    walk(marker_tracks, plot)
    marker_tracks <- plotBrowserTrack(proj_subset,
        groupBy = "Clusters",
        geneSymbol = c("REST", "CUX2"),
        upstream = 100000,
        downstream = 100000
    )
    walk(marker_tracks, plot)
    graphics.off()
    
    proj_subset@cellColData$cluster_group <- proj_subset@cellColData %>%
        str_glue_data("{Clusters}-{Clinical.Dx}")
    
    dx_region <- GRanges(tribble(
        ~seqnames, ~start,      ~end,
        'chr7',     6652819,    6653319,
        'chr2',     120988358,  120988858,
        'chr5',     15936706,   15937206
    ))
    marker_tracks <- plotBrowserTrack(proj_subset,
        groupBy = "cluster_group",
        region = dx_region,
        tileSize = 10
    )
    
    pdf(file.path(out_dir, "marker_tracks_cluster_astro_dx.pdf"), width = 4, height = 20)
    walk(marker_tracks, plot)
    graphics.off()
    
    proj_subset@cellColData$cluster_group <- proj_subset@cellColData %>%
        str_glue_data("{Clusters}-{Clinical.Dx}-{region}")
    marker_tracks <- plotBrowserTrack(proj_subset,
        groupBy = "cluster_group",
        geneSymbol = c("CUX2", "ZMAT4", "REST", "NSF", "RORB", "SAT2B", "NRGN")
    )
    
    pdf(file.path(out_dir, "marker_tracks_cluster_astro_region.pdf"), width = 4, height = 20)
    walk(marker_tracks, plot)
    graphics.off()
    
    proj_subset@cellColData$cluster_group <- proj_subset@cellColData %>%
        str_glue_data("{Clinical.Dx}-{Sample}")
    marker_tracks <- plotBrowserTrack(proj_subset,
        groupBy = "cluster_group",
        geneSymbol = c("CUX2", "ZMAT4", "REST", "NSF", "RORB", "SAT2B", "NRGN")
    )
    
    pdf(file.path(out_dir, "marker_tracks_sample_astro.pdf"), width = 4, height = 20)
    walk(marker_tracks, plot)
    graphics.off()
}

main()
