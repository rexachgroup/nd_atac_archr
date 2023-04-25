# browser tracks centered around certain motifs / peaks.
liblist <- c("tidyverse", "readxl", "Seurat", "ArchR", "GenomicRanges", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
out_dir <- file.path(base_dir, "11_archr_browser_track/")
batchtools <- file.path(out_dir, "batchtools")
archr_project <- list(
    "cluster_insula-c2" = file.path(base_dir, "peak_calling_insula_C2"),
    "subcluster_peak_call_insula-c25" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_insula-c25")
)

peak_tb <- tribble(
    ~name,              ~seqnames,  ~start,     ~end,       ~offset,
    "VCP",              "chr9",     35056064,   35073249,   1150,
    "rs78011262",       "chr17",    45620000,   45835900,   0,
    "rs2045091",        "chr8",     130052000,  130063700,  0,
    "NPTX2",             "chr7",     98616428,   98629868,   0,
)

# plot_param <- tibble(
#     tileSize = c(10, 50, 100)
# )

main <- function() { 
    orig_dir <- getwd()
    setwd(out_dir)
    
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    tb <- tibble(
        project_names = names(archr_project),
        proj_dir = archr_project,
    )
    
    tb <- tb %>% mutate(proj = map(proj_dir, loadArchRProject))
    tb <- tb %>% mutate(peaks = map(proj, getPeakSet))

    # annot <- getGenomeAnnotation(tb$proj[[1]])
    # iterate over every row in peak_tb and plot_param.
    # If offset does not fall in end - start, extend in that direction, then extend 2000 kb.
    peak_coords <- peak_tb %>%
        pmap(function(...) {
            cr <- data.frame(...)
            start <- cr$start
            end <- cr$end
            if (cr$offset < 0) {
                start <- start + cr$offset
            } else if (cr$offset > (end - start)){
                end <- end + (cr$offset - (end - start))
            }
            cr$start <- start - 2000
            cr$end <- end + 2000
            return(cr)
        }) %>%
        bind_rows()

    plot_gr <- GRanges(peak_coords)

    tb <- tb %>%
        mutate(track_obj = pmap(., function(...) {
            cr <- list(...)
            plot_browser_tracks(cr$proj, cr$project_names, plot_gr)
        }))

    pwalk(tb, function(...) {
        cr <- list(...)
        out <- str_glue("{out_dir}/{cr$project_names}_track.pdf")
        pdf(out, width = 10, height = 10)
        tryCatch(print(cr$track_obj), error = print)
        dev.off()
    })
}

# For each archr project and peak loation, call plotBrowserTrack across 
# 1. the whole peak set
# 2. the whole peak set, split by dx
# 2. each cluster, split by dx.
plot_browser_tracks <- function(project, project_name, plot_grange, tileSize = 50) {
    cluster_group_meta <- project@cellColData %>% 
       as.data.frame %>%
       rownames_to_column("barcode") %>%
       group_by(Clusters) %>%
       group_nest(keep = T)
            
    map(1:length(plot_gr), function(range_i) {
        plot_range <- plot_gr[range_i, ]

        data_title <- str_glue("{project_name} {plot_range$name} all cluster atac signal sum")
        data_plot <- wrap_plots(plotBrowserTrack(project, region = plot_range, tileSize = tileSize)) + plot_annotation(title = data_title)

        dx_title <- str_glue("{project_name} {plot_range$name} all cluster atac signal sum, split dx") 
        dx_split_plot <- wrap_plots(
            plotBrowserTrack(project, region = plot_range, groupBy = "Clinical.Dx", tileSize = tileSize)
        ) + plot_annotation(title = dx_title)

        cluster_split_plots <- pmap(cluster_group_meta, function(...) {
            m <- list(...)
            subset_proj <- project[, m$barcode]
            title <- str_glue("{project_name} {plot_range$name} {m$Clusters} only atac signal sum, split dx") 
            track_obj <- plotBrowserTrack(subset_proj, region = plot_range, groupBy = "Clinical.Dx", tileSize = tileSize)
            return(wrap_plots(track_obj) + plot_annotation(title = title))
        })
        return(list(data_plot, dx_split_plot, cluster_split_plots))

    })

}

main()
