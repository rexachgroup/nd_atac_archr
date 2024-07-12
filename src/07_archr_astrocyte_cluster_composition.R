liblist <- c("tidyverse", "readxl", "Seurat", "ArchR")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, "clustering_harmony/")
out_dir <- file.path(base_dir, "07_archr_cluster_composition/")
main <- function() { 
    dir.create(out_dir)
    # astrocyte clusters composition
    project <- loadArchRProject(archr_project)
    astro_clusters <- paste("C", 18:25, sep = "")
    proj_meta <- project@cellColData %>%
        as.data.frame() %>%
        as.tibble() %>%
        mutate(Clinical.Dx = fct_relevel(Clinical.Dx, "Control"))
    astro_clusters <- intersect(paste("C", 18:25, sep = ""), proj_meta$Clusters)

    library_counts <- proj_meta %>%
        filter(Clusters %in% astro_clusters) %>%
        group_by(Clinical.Dx, Autopsy.ID, region) %>%
        dplyr::summarize(library_ct = n())
    
    subcluster_counts <- proj_meta %>%
        filter(Clusters %in% astro_clusters) %>%
        group_by(Clinical.Dx, Autopsy.ID, Clusters, region) %>%
        dplyr::summarize(subcluster_ct = n())

    prop_tb <- subcluster_counts %>%
        left_join(library_counts, by = c("Clinical.Dx", "Autopsy.ID", "region")) %>%
        mutate(subcluster_prop = subcluster_ct / library_ct)

    pdf(file.path(out_dir, "astro_prop.pdf"))
    ggplot(prop_tb, aes(x = Clusters, y = subcluster_prop, fill = Clinical.Dx)) +
        geom_boxplot() + 
        geom_point(position = position_jitterdodge())
    graphics.off()

    astro_clusters_prop_test <- prop_tb %>%
        group_by(Clusters) %>%
        group_nest() %>%
        mutate(prop_test = pmap(., function(...) {
            cr <- list(...)
            dat <- cr$data
            form <- "subcluster_prop ~ Clinical.Dx"
            lm_tidy <- broom::tidy(lm(form, data = dat))
            lm_tidy <- mutate(lm_tidy,
                p.value.adj = p.adjust(p.value, method = "fdr"),
                form = form
            )
        })) %>% unnest(prop_test) %>%
        filter(str_detect(term, "Clinical.Dx"))

    astro_clusters_prop_test %>%
        select(!data) %>%
        write_csv(file.path(out_dir, "astro_clusters_prop_test.csv"))
}

if (!interactive()) main()
