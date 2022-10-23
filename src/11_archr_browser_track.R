# browser tracks centered around certain motifs / peaks.
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "GenomicRanges", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
out_dir <- file.path(base_dir, "11_archr_browser_track/")
batchtools <- file.path(out_dir, "batchtools")
archr_project <- list(
    "precg_C2" = file.path(base_dir, "10_tf_transfac_c2", "chrom_var_precg_C2"),
    "insula_C2" = file.path(base_dir, "10_tf_transfac_c2", "chrom_var_insula_C2"),
    "subcluster_peak_call_insula-c25" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_insula-c25"),
    "subcluster_peak_call_precg-c25" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_precg-c25"),
    "subcluster_peak_call_insula-c6789" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_insula-c6789"),
    "subcluster_peak_call_precg-c6789" = file.path(base_dir, "09_archr_subcluster_peak_calling", "peak_call_precg-c6789")

    
)
peak_tb <- tibble(
    seqnames = c("chr7", "chr7"),
    start = c("98616427", "98616927"),
    end = c("98616428", "98616928")
)

plot_param <- tibble(
    tileSize = c(10, 50, 250)
)

main <- function() {
    
    orig_dir <- getwd()
    setwd(out_dir)
    
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(batchtools)
    } else {
        reg <- loadRegistry(batchtools, writeable = TRUE)
    }

    tb <- tibble(
        project_names = names(archr_project),
        proj_dir = archr_project,
    )
    sub_gr <- GRanges(peak_tb)
    #     peak_dx_tb <- map(dx_peaks_tb, function(tb) {
    #         unlist(tb) %>%
    #             as.data.frame %>%
    #             rownames_to_column("dxpeak") %>%
    #             as_tibble
    #     }) %>% bind_rows
    
    tb <- tb %>% mutate(proj = map(proj_dir, loadArchRProject))
    tb <- tb %>% mutate(peaks = map(proj, getPeakSet))
    tb <- tb %>%
        mutate(peak_sub = map(peaks, function(p) {
            peaks <- subsetByOverlaps(p, sub_gr, ignore.strand = T)
            flank(peaks, width = 5000, both = T)
        }))

    tb <- tb %>%
        mutate(track_obj = pmap(., function(...) {
            cr <- list(...)
            pmap(plot_param, function(...) {
                pr <- list(...)
                browser()
                title <- str_glue("{cr$project_names} {deparse(pr)} atac signal sum, split dx")
                print(title)
                track_obj <- plotBrowserTrack(cr$proj, region = cr$peak_sub, groupBy = "Clinical.Dx", ...)
                wrap_plots(track_obj) + plot_annotation(title = title)
            })
        }))

    pwalk(tb, function(...) {
        cr <- list(...)
        out <- str_glue("{out_dir}/{cr$project_names}_track.pdf")
        pdf(out)
        print(cr$track_obj)
        dev.off()
    })
}

