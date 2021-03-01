liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

SAMPLE_META <- "/geschwindlabshares/lchenprj01/nucseq_combined/data/snATAC_metadata_summary_2021_d.xlsx"
out_dir <- "../data/archr/atac-2020-all"
data_dir <- normalizePath(file.path(out_dir, "data"))
plot_dir <- normalizePath(file.path(out_dir, "plot"))
archr_project <- normalizePath(file.path(out_dir, "preprocess"))
sample_table <- normalizePath(file.path(out_dir, "preprocess", "sample_tb.rds"))

main <- function() {
    addArchRThreads(threads = 32)
    addArchRGenome("hg38")
    writeLines(str_glue("load {archr_project}"))
    project <- loadArchRProject(path = archr_project)
    #orig_dir <- getwd()
    #setwd(archr_project)

    sample_meta <- read_xlsx(SAMPLE_META)
    sample_tb <- readRDS(sample_table)
    #check size 
    paste0("Memory Size = ", round(object.size(project) / 10^6, 3), " MB")
    #available matrices
    getAvailableMatrices(project)
    #access cell names, sample names, tss score, 
    head(project$cellNames); head(project$Sample); quantile(project$TSSEnrichment)

    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")
    filtered_cells <- cell_data %>%
        filter(nFrags > 1000, nFrags < 100000, TSSEnrichment > 2, BlacklistRatio < 0.1) %>%
        pluck("cell_id")
    writeLines(str_glue("{length(filtered_cells)} / {nrow(cell_data)} cells after filtering"))

    # subset + reload cellColData
    writeLines("subset")
    project <- project[filtered_cells, ]
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")

    # TSS vs. Fragments
    gg_tss <- ggPoint(
        x = log10(cell_data$nFrags),
        y = cell_data$TSSEnrichment,
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = c(log10(500), quantile(log10(cell_data$nFrags), probs = 0.99)),
        ylim = c(0, quantile(cell_data$TSSEnrichment, probs = 0.99))
    ) + 
    geom_hline(yintercept = 4, lty = "dashed") + 
    geom_vline(xintercept = 3, lty = "dashed")

    #ggsave(file.path(plot_dir, "TSS-vs-Frags.pdf"), gg_tss)
    plotPDF(gg_tss, "TSS-vs-Frags.pdf", ArchRProj = project, addDOC = FALSE)

    # Ridge / violin plots for all. Computes fragment / TSS sizes
    writeLines("plotFragmentSizes / plotTSSEnrichment")
    #fragment_sizes <- plotFragmentSizes(project, returnDF = TRUE, threads = 128)
    #saveRDS(fragment_sizes, file.path(archr_project, "fragment_sizes.rds"))
    p1 <- plotFragmentSizes(ArchRProj = project)
    p2 <- plotTSSEnrichment(ArchRProj = project)

    plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)
    cell_data <- as_tibble(project@cellColData, rownames = "cell_id")
    write_csv(cell_data, file.path(archr_project, "cell_data.csv"))
    
    # Ridge / violin plots per sample
    writeLines("plotGroups")
    p1 <- plotGroups(
        ArchRProj = project,
        groupBy = "Sample", 
        colorBy = "cellColData", 
        name = "TSSEnrichment",
        plotAs = "ridges"
    )

    p2 <- plotGroups(
        ArchRProj = project, 
        groupBy = "Sample", 
        colorBy = "cellColData", 
        name = "TSSEnrichment",
        plotAs = "violin",
        alpha = 0.4,
        addBoxPlot = TRUE
    )

    p3 <- plotGroups(
        ArchRProj = project, 
        groupBy = "Sample", 
        colorBy = "cellColData", 
        name = "log10(nFrags)",
        plotAs = "ridges"
    )

    p4 <- plotGroups(
        ArchRProj = project, 
        groupBy = "Sample", 
        colorBy = "cellColData", 
        name = "log10(nFrags)",
        plotAs = "violin",
        alpha = 0.4,
        addBoxPlot = TRUE
    )
    p5 <- plotGroups(
        ArchRProj = project,
        groupBy = "Sample",
        colorBy = "cellColData",
        name = 
    )

    plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = project, addDOC = FALSE, width = 10, height = 10)
}

if (!interactive()) {
    main()
}
