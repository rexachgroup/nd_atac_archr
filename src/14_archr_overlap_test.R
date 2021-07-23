# Test overlap of scenic / regulon gene lists for each transcription factor's gene target.

base_dir <- normalizePath("../data/archr/atac-2020-all/")
liblist <- c("tidyverse", "progress", "batchtools", "readxl", "Seurat", "ArchR")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))
base_out_dir <- file.path(base_dir, "14_archr_motif_overlap_test")
cluster_args <- file.path(base_dir, "11_archr_peak_overlaps", "args_tb.rds")
scenic_excitatory <- normalizePath("../data/excitatory_preCGregulons_genes_GRNBoost_weights.csv")
scenic_microglia <- normalizePath("../data/microglia_preCGregulons_genes_GRNBoost_weights.csv")

main <- function() {
    args_tb <- readRDS(cluster_args)
    args_tb <- filter(args_tb, project_names == "precg_C7")
    args_tb <- mutate(args_tb,
        proj_dir = out_dir,
        out_dir = file.path(base_out_dir, project_names)
    )

    args_tb <- args_tb %>% mutate(fishers_tb = pmap(., function(...) {
        cr <- list(...)
        project <- loadArchRProject(cr$proj_dir)
        scenic_tb <- read_csv(scenic_microglia)
        archr_tb <- motif_match_nearest_gene(project)

        archr_tfs <- unique(archr_tb$tf)
        scenic_tfs <- unique(scenic_tb$TF)
        writeLines(str_glue("{length(archr_tfs)}, {length(scenic_tfs)}"))

        test_tb <- cross_df(list(archr = archr_tfs, scenic = scenic_tfs))
        writeLines(str_glue("{nrow(test_tb)}"))

        archr_gene_tb <- select(archr_tb, tf, nearestGene) %>% nest(archr_genes = nearestGene)
        scenic_gene_tb <- select(scenic_tb, TF, gene) %>% nest(scenic_genes = gene)
        background_size <- length(union(scenic_tb$gene, archr_tb$nearestGene))

        writeLines(str_glue("{background_size}"))
        test_gene_tb <- inner_join(test_tb, archr_gene_tb, by = c("archr" = "tf")) %>%
            inner_join(scenic_gene_tb, by = c("scenic" = "TF"))
        pb <- progress_bar$new(total = nrow(test_tb))
        test_gene_tb <- test_gene_tb %>% mutate(fishers = pmap(., function(...) {
            cr <- list(...)
            pb$tick()
            archr_genes <- cr$archr_genes %>% pluck("nearestGene") %>% unique
            scenic_genes <- cr$scenic_genes %>% pluck("gene") %>% unique
            gene_union <- union(archr_genes, scenic_genes)
            gene_intersect <- intersect(archr_genes, scenic_genes)
            cont_tbl <- matrix(
                    c(background_size - length(gene_union),
                    length(archr_genes) - length(gene_intersect),
                    length(scenic_genes) - length(gene_intersect),
                    length(gene_intersect)),
                    ncol = 2
                )
            try(fisher.test(cont_tbl, alternative = "greater"))
        }))

        return(test_gene_tb)
    }))
    saveRDS(args_tb, file.path(base_out_dir, "args_tb.rds"))
    args_tb <- unnest(args_tb, fishers_tb)
    args_tb <- mutate(args_tb, 
            fisher_estimate = map_dbl(fishers, ~ifelse(is.list(.), .$estimate, 0.0)),
            fisher_pval = map_dbl(fishers, ~ifelse(is.list(.), .$p.value, 1.0))
        )
    args_tb %>%
        select(project_names, archr, scenic, fisher_estimate, fisher_pval) %>%
        arrange(fisher_pval) %>%
        write_csv(file.path(base_out_dir, "archr_scenic_fisher_ovarlap.csv"))
}

motif_match_nearest_gene <- function(project) {
    motif_matches <- getMatches(project)
    peak_rowdata <- as_tibble(rowData(motif_matches))
    peak_rowdata$peak_dx <- rownames(rowData(motif_matches))
    motif_peak_idx <- apply(assay(motif_matches), 2, which)

    tf_peakinfo_tb <- bind_rows(map(motif_peak_idx, function(ids) {
        filter(peak_rowdata, idx %in% ids, peakType == "Promoter")
    }), .id = "tf")
    return(tf_peakinfo_tb)
}

