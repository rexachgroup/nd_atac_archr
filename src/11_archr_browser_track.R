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

peak_tb <- tribble(
    ~name,      ~seqnames,  ~start,     ~end,       ~offset,
    "CHMP2B",   "chr3",     87227271,   87255548,   398,
    "PSAP",     "chr10",    71816298,   71851375,   67,
    "PSAP",     "chr10",    71816298,   71851375,   -9263,
    "C9orf72",  "chr9",     27546545,   27573866,   -346,
    "VCP",      "chr9",     35056064,   35056064,  1150
)

plot_param <- tibble(
    tileSize = c(50)
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
    
    tb <- tb %>% mutate(proj = map(proj_dir, loadArchRProject))
    tb <- tb %>% mutate(peaks = map(proj, getPeakSet))

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
            pmap(plot_param, function(...) {
                pr <- list(...)
                map(1:length(plot_gr), function(range_i) {
                    gr <- plot_gr[range_i,]
                    title <- str_glue("{cr$project_names} {deparse(pr)} {gr$name} atac signal sum, split dx")
                    track_obj <- plotBrowserTrack(cr$proj, region = gr, groupBy = "Clinical.Dx", ...)
                    print(title)
                    return(wrap_plots(track_obj) + plot_annotation(title = title))
                })
            })
        }))

    pwalk(tb, function(...) {
        cr <- list(...)
        out <- str_glue("{out_dir}/{cr$project_names}_track.pdf")
        pdf(out, width = 10, height = 10)
        tryCatch(print(cr$track_obj), error = print)
        dev.off()
    })
}

main()
