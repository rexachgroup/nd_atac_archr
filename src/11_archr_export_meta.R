# export metadata / sample-level summary for external processing.
liblist <- c("readxl", "Seurat", "ArchR", "tidyverse")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))
base_dir <- normalizePath("../data/archr/atac-2020-all")
out_dir <- file.path(base_dir, "11_archr_meta_export")
CLUSTER_ARCHR <- file.path(base_dir,"clustering_harmony")
SUBCLUSTER_NAMES <- c("precg_C7", "insula_C7", "precg_C2", "insula_C2")
SUBCLUSTER_ARCHR <- setNames(file.path(base_dir, paste0("peak_calling_", SUBCLUSTER_NAMES)), SUBCLUSTER_NAMES)

main <- function() {
    dir.create(out_dir)
    oldwd <- getwd()
    setwd(out_dir)
    cluster_obj <- loadArchRProject(CLUSTER_ARCHR, showLogo = FALSE)

    subcluster_objs <- tibble(name = SUBCLUSTER_NAMES, path = SUBCLUSTER_ARCHR)

    cluster_obj@cellColData %>%
        as_tibble %>%
        group_by(Sample, region, Clinical.Dx) %>%
        slice_head(n = 1) %>%
        group_by(Clinical.Dx, region) %>%
        dplyr::summarize(m = sum(Sex == 'M'), f = sum(Sex == 'F'), pmi = mean(PMI)) %>%
        arrange(region, Clinical.Dx) %>%
        write_csv(file.path(out_dir, "sample_summary.csv"))

    # sample_meta
    cluster_obj@cellColData %>%
        as_tibble %>%
        group_by(Sample, region, Clinical.Dx) %>%
        select(Sample, Clinical.Dx, Autopsy.ID, region, Age, Sex, PMI) %>%
        slice_head(n = 1) %>%
        write_csv(file.path(out_dir, "sample_meta.csv"))
    # sample_ct
    cluster_obj@cellColData %>%
        as_tibble %>%
        group_by(Sample, region, Clinical.Dx) %>%
        summarize(n = n()) %>%
        write_csv(file.path(out_dir, "sample_ct.csv"))


    subcluster_objs <- subcluster_objs %>%
        mutate(obj = map(path, loadArchRProject, showLogo = F))

    subcluster_objs <- subcluster_objs %>%
        mutate(meta = map(obj, ~.x@cellColData %>% as_tibble))
    
    subcluster_objs <- subcluster_objs %>%
        mutate(dx_sum = map(meta, function(x) {
            x %>%
                group_by(Clinical.Dx) %>%
                dplyr::summarize(n = n())
        }))
}
