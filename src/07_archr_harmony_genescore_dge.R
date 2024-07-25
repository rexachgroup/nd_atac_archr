liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "pheatmap", "future.apply")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

base_dir <- normalizePath("../data/archr/atac-2020-all")
archr_project <- file.path(base_dir, "clustering_harmony/")
out_dir <- file.path(base_dir, "07_archr_harmony_genescore_dge/")

main <- function() {
    addArchRThreads(16)
    addArchRGenome("hg38")
    dir.create(out_dir)
    writeLines(str_glue("load {archr_project}"))
    project <- loadArchRProject(path = archr_project)
    
    cells <- project@cellColData %>%
        as.data.frame
    cells_f <- cells %>%
        rownames_to_column("barcode") %>%
        filter(Clusters == "C22")
    subsetArchRProject(project, cells = cells_f$barcode, force = T)
    proj_subset <- loadArchRProject("ArchRSubset")
    gscm <- getMatrixFromProject(proj_subset)
    gsc_a <- as.matrix(assays(gscm)$GeneScoreMatrix)

    proj_subset$Clinical.Dx <- fct_relevel(proj_subset$Clinical.Dx, "Control")

    rownames(gsc_a) <- gscm@elementMetadata$name

    plan(multicore, workers = 8)
    gene_fit <- future_apply(gsc_a, MARGIN = 1, function(vec) {
          broom::tidy(lm(vec ~ Clinical.Dx + TSSEnrichment, data = cells_f))
    })
    gene_fit_tb <- gene_fit %>% bind_rows(.id = "gene") %>% 
        as_tibble() %>%
        mutate(formula = "expr ~ Clinical.Dx + TSSEnrichment") %>%
        group_by(term) %>%
        dplyr::mutate(p.adj.n = n(), p.adj = p.adjust(p.value, n = p.adj.n))

    write_csv(gene_fit_tb, file.path(out_dir, "c22_genescore_lm.csv"))
    
}

if (!interactive()) main()
