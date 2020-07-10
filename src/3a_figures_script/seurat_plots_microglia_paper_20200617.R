# vim: set shiftwidth=2
# Damon Polioudakis
# 2020-02-06
# Plots of seurat analysis

# Must load modules:
#   module load R/3.6.0

# submit as job:
# qsub -N sp_mg_pr -l h_data=128G,h_rt=12:00:00,highp qsub_r_script.sh -p seurat_plots_microglia_paper_20200617.R
################################################################################

rm(list = ls())
set.seed(27)

require(methods)
require(Seurat)
require(Matrix)
require(cowplot)
require(viridis)
require(RColorBrewer)
require(ggcorrplot)
require(tidyverse)
require(gridExtra)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
sessionInfo()

## inputs
#in_seurat <- "/u/home/d/dpolioud/project-geschwind/nucseq_nd/analysis/seurat/20191116/percent_mt_8/reg_umi/p1-5_c1-5_filtered.rdat"
in_seurat <- "/u/project/geschwind/chenlo/nucseq_combined/ext/test/single_cell_seurat/2020_06_10/testing_pfc_processed.rds"
#in_seurat_raw <-
#  "/u/home/d/dpolioud/project-geschwind/nucseq_nd/analysis/seurat/20191116/percent_mt_8/reg_umi/p1-5_c1-5_raw.rdat"
resources_dir <- "/u/project/geschwind/dpolioud/nucseq_nd/resources"
marker_genes_tb <- read_csv(file.path(resources_dir,
  "cluster_markers_mixed_20181019.csv"))
marker_genes_refined_tb <- read_csv(file.path(resources_dir,
  "20191028_cell_markers_refined.csv"))
astro_markers_tb <- read_csv(file.path(resources_dir,
  "AstrocyteMarkers_Brie_20190110.csv"))
#in_exc_marker <- file.path(resources_dir, "excitatory_markers_20191023.csv")
polo_ad_gwas_genes_tb <- read_csv(file.path(resources_dir, "grubman_2019_st5_AD_GWAS_genes_ad.csv"))
#metadata_corrections_tb <- read_csv(file.path(resources_dir, "2020.06.13_corrected_metadata_jessica.csv"))

## Variables
date <- format(Sys.Date(), "%Y%m%d")
script_name <- paste0("seurat_plots_microglia_paper.R ", date)
graph_subtitle <- "testing_pfc_processed PFC1-3, 5% MT filter"

## outputs
out_graph <- file.path("../../data/05a-seurat-graphs/seurat_graph-")
out_table <- file.path("../../data/05a-seurat-graphs/seurat_table-")

# make directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
################################################################################

### functions

main_function <- function(){

  print("main_function")

  print("loading seurat object...")
  nd_so <- readRDS(in_seurat)
  print("done loading seurat object...")
  
  # alias variables
  nd_so[["number_umi"]] <- nd_so[["nCount_RNA"]]
  nd_so[["number_genes"]] <- nd_so[["nFeature_RNA"]]
  nd_so[["percent_mito"]] <- nd_so[["percent.mt"]]
  nd_so[["pmi_h"]] <- nd_so[["PMI"]] 
  nd_so[["rin"]] <- nd_so[["RIN"]]
  nd_so[["weight_g"]] <- nd_so[["weight"]]
  nd_so[["age"]] <- nd_so[["Age"]]
  nd_so[["sex"]] <- nd_so[["Sex"]]
  nd_so[["region"]] <- "PFC"
  nd_so@meta.data <- select(nd_so@meta.data, -c("nCount_RNA", "nFeature_RNA", "percent.mt", "PMI", "RIN", "weight", "Age", "Sex"))
  #nd_so[["finalsite"]] <- nd_so[["Site"]]

  # fix sample swaps
  metadata_to_fix <- c("clinical_dx", "region", "age", "sex", "pmi_h", "rin", "weight_g")
  # metadata_corrections_tb <- read_csv("../metadata/2020.06.13_corrected_metadata_jessica.csv")
  # check before fixing
  stopifnot(all(metadata_to_fix %in% colnames(nd_so@meta.data)))
  nd_so[[c("library_id", metadata_to_fix)]] %>%
      distinct() %>% data.frame() %>% print()

  # fix
#   for(x in metadata_to_fix){
#     corrected_tb <- tibble(
#       library_id = metadata_corrections_tb$library_id,
#       corrected_variable = metadata_corrections_tb[[x]])
#     nd_so[[x]] <- nd_so[[c("library_id", x)]] %>%
#       left_join(., corrected_tb, by = "library_id") %>%
#       mutate(x = ifelse(
#         ! is.na(corrected_variable),
#         corrected_variable,
#         .[[x]])) %>%
#       pull(x)
#   }
# check after fixing
#   nd_so[[c("library_id", metadata_to_fix)]] %>%
#     filter(grepl("P5", library_id)) %>%
#     distinct() %>% data.frame() %>% print()
# 
#   nd_so$cluster_ids <- nd_so[["RNA_snn_res.0.6"]]
#   nd_so$cell_ids <- nd_so[[]] %>% rownames
# 
  # subset to preCG microglia
  #   subset_so <- subset(nd_so, subset = region == "preCG")
  # subset_so <- subset(nd_so, subset = prep %in% c(1,2,3,4))
  #   subset_so <- subset(subset_so, subset = clinical_dx %in% c("bvFTD", "Control"))
  subset_so <- nd_so
  print(table(subset_so[["clinical_dx"]]))


  # remove P1_7
  #   subset_so <- subset(
  #     subset_so, subset = library_id == c("P1_7"), invert = TRUE)
  #   subset_so[[c("library_id")]] %>% table() %>% print()

  # add metadata variable that labels each cell TRUE / FALSE == microglia
  clinical_dx_microglia <- subset_so[["clinical_dx"]]
  clinical_dx_microglia[subset_so[["cluster_cell_type"]] != "microglia"] <- "Other cell type"
  subset_so[["clinical_dx_microglia"]] <- clinical_dx_microglia

  # add metadata variable that labels cluster TRUE / FALSE == microglia
  microglia_clusters <- subset_so[["cluster_cell_type"]]
  microglia_clusters[subset_so[["cluster_cell_type"]] != "microglia"] <- "other cell type"
  subset_so[["microglia_clusters"]] <- microglia_clusters
  subset_so[["microglia_clusters"]] %>% table() %>% print()

  # subset to microglia for some plots
  microglia_so <- subset(subset_so, subset = cluster_cell_type == "microglia")
  microglia_so[[c("library_id")]] %>% table() %>% print()
  # output table of cell number by library id
  microglia_so[[c("library_id")]] %>% tibble() %>% count(library_id) %>%
    write_csv(., path = paste0(out_table, "FTD_control_preCG_microglia_cell_number_by_library_id.csv"))
  # run PCA
  microglia_so <- RunPCA(microglia_so, online.pca = TRUE,
    npcs = 100, do.print = FALSE, verbose = FALSE, approx = FALSE)

  ## plots

  plot_qc_summary_stats_violin_boxplots(seurat_obj = subset_so)

  plot_qc_summary_stats_violin_boxplots(
    seurat_obj = microglia_so, out_graph_suffix = "_microglia")

  ##### I would need to filter the raw data down to just FTD and control too
  # plot_qc_metrics(seurat_obj = subset_so, in_seurat_raw = in_seurat_raw)

  plot_r_squared_matrix_of_pcs_and_metadata(seurat_obj = subset_so, pcs = 1:10)
  plot_variance_explained_by_expression_pc(
    seurat_obj = subset_so, pcs = 1:10)

  plot_r_squared_matrix_of_pcs_and_metadata(
    seurat_obj = microglia_so, pcs = 1:10, out_file_suffix = "_microglia")
  plot_variance_explained_by_expression_pc(
    seurat_obj = microglia_so, pcs = 1:10, out_file_suffix = "_microglia")

  plot_metadata_by_clinical_dx_numerical_boxplot_wrapper(
    seurat_obj = subset_so)

  plot_metadata_by_clinical_dx_numerical_boxplot_wrapper(
    seurat_obj = microglia_so, out_graph_suffix = "_microglia")

  plot_metadata_by_clinical_dx_factor_bar_plot_wrapper(
    seurat_obj = subset_so)

  plot_metadata_by_clinical_dx_factor_bar_plot_wrapper(
    seurat_obj = microglia_so, out_graph_suffix = "_microglia")

  # plot_cell_type_by_cluster_stacked_bar_plot(
    # seurat_obj = subset_so, cluster_col_name = "cluster_ids")

  plot_number_of_cells_per_cluster_bar_plot(
    seurat_obj = subset_so,
    cluster_col_name = "cluster_ids",
    out_graph_file_type = ".pdf"
  )

  plot_metadata_by_cluster_stacked_bar_plot(
    seurat_obj = subset_so,
    cluster_col_name = "cluster_ids",
    out_graph_file_type = ".pdf"
  )

  plot_metadata_by_cluster_percent_stacked_bar_plot(
    seurat_obj = subset_so,
    cluster_col_name = "cluster_ids",
    out_graph_file_type = ".pdf"
  )

  # plot_metadata_by_variable_stacked_barplot
  {
    var_string <- "library_id"
    
    metadata_vars <- c(
      "cell_type", 
      "RNA_snn_res.0.4",
      "RNA_snn_res.0.5",
      "RNA_snn_res.0.6",
      "RNA_snn_res.0.7",
      "RNA_snn_res.0.8",
      "RNA_snn_res.0.9",
      "RNA_snn_res.1"
    )
    seurat_obj <- subset_so
    
    # Select and pivot_longer all non-numeric variables other than var_string.
    metadata_tb <- seurat_obj@meta.data %>%
      as_tibble %>%
      rownames_to_column("UMI") %>%
      mutate(group = .data[[var_string]]) %>%
      select(group, matches(metadata_vars)) %>%
      select_if(negate(is.numeric)) %>%
      pivot_longer(cols = -group) %>%
      glimpse()

    legend_nrow = ceiling(length(unique(metadata_tb$name)) / 2)
    

    tb_pct <- metadata_tb %>%
      group_nest(name) %>%
      pmap(function(name, data, ...){
        ggtb <- data %>%
          group_by(group) %>%
          count(value) %>%
          mutate(
                 percent = n / sum(n),
                 variable = name,
                 group = factor(group)
                 ) %>%
          glimpse()
        ggplot(ggtb, aes(x = group, y = percent, fill = value)) +
          geom_bar(stat = "identity") +
          {
            valuect <- length(unique(ggtb$value))
            if (valuect > 11){
              scale_fill_manual(name = name, values = 
                                colorRampPalette(brewer.pal(n = 11, name = "Set3"))(valuect)
              )
            } else {
              scale_fill_brewer(name = name, palette = "Set3")
            }
          } +
         guides(fill = guide_legend(ncol = ceiling(valuect/legend_nrow))) + 
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
      })
      meta_grid <- plot_grid_wrapper(tb_pct, ncol = 2, title = make_plot_title("Metadata by library_id"))
      ggsave(paste0(out_graph, "metadata_by_library_percent_stackedbar.png"), meta_grid,
             width = 24, height = 4 * legend_nrow, limitsize = FALSE)
  }

  # T-cell population plots.
  {
    genes <- rownames(subset_so@assays$RNA)

    immune_genes <- unique(c("CD3E", "LCK", "TRAC", "TRAJ16", "CD8A", "CCR7", "CD8A", "CD44", "CD8A", "IFNG", "CD8A", "CD69", "ITGAE", "CD4", "IL7R", "FOXP3", "CTLA4", "CD19", "IGHM", "MS4A1", "CD40", "CD38", "SDC1",  "CD27"))
    immune_genes %in% genes
    plot_genes_of_interest_expression_dim_reduction(
      seurat_obj = subset_so, reduction = "umap", out_suffix = "immune_cell_markers", genes = immune_genes)
    plot_number_of_cells_expressing_genes_box_jitter_plot(
      seurat_obj = subset_so, genes = immune_genes)
  }

  {
    gene_tb <- FetchData(subset_so, vars = c(immune_genes, "clinical_dx", "library_id")) %>%
      as_tibble() %>%
      pivot_longer(cols = matches(immune_genes), names_to = "gene", values_to = "expression")
    gene_sum <- gene_tb %>%
      group_by(library_id, gene) %>%
      summarize(clinical_dx = unique(clinical_dx), cells_expressed = sum(expression > 0), total_cells = n(),  pct = cells_expressed / total_cells, .groups = "keep") %>%
      arrange(library_id, gene)
  }

  # plot UMAPs with microglia clusters colored gray
  reduction <- "umap"
  var_name <- "microglia_clusters"
  plot_dim_reduction_colored_by_variable(
    dim_1 = Embeddings(object = subset_so, reduction = reduction)[ ,"UMAP_1"],
    dim_2 = Embeddings(object = subset_so, reduction = reduction)[ ,"UMAP_2"],
    variable_value = subset_so[[var_name]][ ,var_name],
    title = paste0(var_name),
    legend_title = paste0(var_name),
    size = 0.1,
    alpha = 0.5,
    guide_size = 2
    ) +
    scale_color_manual(values = c("#00BA38", "lightgrey")) +
    ggtitle(make_plot_title(paste0(
      "Cells colored by microglia clusters",
      "\nDimensionality reduction: ", reduction)))
  ggsave(paste0(out_graph, var_name, "_", reduction, ".png")
    , width = 6, height = 4.5, dpi = 200, limitsize = FALSE)

  # plot UMAPs with microglia colored by dx and all other cell types gray
  reduction <- "umap"
  var_name <- "clinical_dx_microglia"
  plot_dim_reduction_colored_by_variable(
    dim_1 = Embeddings(object = subset_so, reduction = reduction)[ ,"UMAP_1"],
    dim_2 = Embeddings(object = subset_so, reduction = reduction)[ ,"UMAP_2"],
    variable_value = subset_so[[var_name]][ ,var_name],
    title = paste0(var_name),
    legend_title = paste0(var_name),
    size = 0.1,
    alpha = 0.5,
    guide_size = 2
    ) +
    scale_color_manual(values = c("#e41a1c", "#337eb8", "#4daf4a", "lightgrey", "#984ea3", "#ff7f00")) +
    ggtitle(make_plot_title(paste0(
      "Cells colored by clinical_dx for microglia",
      "\nDimensionality reduction: ", reduction)))
  ggsave(paste0(out_graph, var_name, "_", reduction, ".png")
    , width = 6, height = 4.5, dpi = 200, limitsize = FALSE)

  # plot UMAPs colored by mean expression of groups of cell type markers
  plot_marker_refined_expression_dim_reduction(
    seurat_obj = subset_so, reduction = "umap")


  ## output data tables

  # make table of bvFTD and and control preCG microglia normalized expression matrix
  # save as rdat
  microglia_so <- subset(subset_so, subset = cluster_cell_type == "microglia")
  normalized_expression_matrix <- GetAssayData(microglia_so, slot = "data")
  print("dimensions of FTD preCG microglia expression matrix:")
  normalized_expression_matrix %>% dim() %>% print()
  save(normalized_expression_matrix, file = paste0(out_table, "FTD_control_preCG_microglia_expression_matrix.rdat"))
  # save as csv
  GetAssayData(microglia_so, slot = "data") %>%
    as_tibble(rownames = "gene") %>%
    write_csv(., path = paste0(out_table, "FTD_control_preCG_microglia_expression_matrix.csv"))
  # make table of bvFTD and and control preCG microglia normalized expression for genes of interest
  FetchData(microglia_so, slot = "data",
      vars = c(
        "B2M",
        "IFNGR1",
        "IRF8",
        "JAK2",
        "RUNX1",
        "SPP1",
        "TLR7",
        "CD74",
        "TRIM69",
        "IFIT3",
        "IFNGR2",
        "IL18",
        "CD276",
        "CTSS",
        "CSTB",
        "PRDX1",
        "ALCAM",
        "ADARB1",
        "DNASE1",
        "TRIM5",
        "HLA-B")
      ) %>%
    t() %>%
    as_tibble(rownames = "gene") %>%
    write_csv(., path = paste0(out_table, "FTD_control_preCG_microglia_genes_of_interest_expression_matrix.csv"))
  # make table of bvFTD and and control preCG microglia cell_ids and clinical_dx
  FetchData(microglia_so, vars = c("clinical_dx"))[] %>%
    as_tibble(rownames = "cell_id") %>%
    write_csv(., path = paste0(out_table, "FTD_control_preCG_microglia_cell_id_clinical_dx.csv"))
  # make table of bvFTD and and control preCG microglia metadata
  FetchData(microglia_so, vars = c(colnames(microglia_so[[]])))[] %>%
    as_tibble(rownames = "cell_id") %>%
    write_csv(., path = paste0(out_table, "FTD_control_preCG_microglia_metadata.csv"))
  # write table of immune cell markers
  write_csv(gene_sum, path = paste0(out_table, "immune_cell_expression_count.csv"))

}
################################################################################

plot_qc_summary_stats_violin_boxplots <- function(
  seurat_obj, out_graph_suffix = ""){

  print("plot_qc_summary_stats_violin_boxplots()")

  metadata <- c("cell_ids", "number_genes", "number_umi", "percent_mito"
    , "library_id")

  # violin boxplots
  summary_stats_tb <- seurat_obj[[metadata]] %>%
    gather(variable, value, -cell_ids, -library_id) %>%
      group_by(variable) %>%
      summarize(
        value_mean = round(mean(value, na.rm = TRUE), 1)
        , value_median = round(median(value, na.rm = TRUE),1)
      ) %>%
      mutate(label = paste0(
        "mean:\n", value_mean, "\nmedian:\n", value_median))
  seurat_obj[[metadata]] %>%
    gather(variable, value,  -cell_ids, -library_id) %>%
    group_by(variable) %>%
    ggplot(aes(x = variable, y = value, fill = variable)) +
      facet_wrap(~variable, scales = "free", ncol = 5) +
      geom_violin(fill = "lightgrey", color = "lightgrey") +
      geom_boxplot(width = 0.2, outlier.size = 0.1) +
      theme(
        legend.position = "none"
        , axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data = summary_stats_tb, aes(
        x = variable, y = Inf, label = label), vjust = 1, hjust = 2) +
      ggtitle(make_plot_title("QC metrics"))
    ggsave(paste0(
        out_graph, "qc_metrics_violin_boxplots", out_graph_suffix, ".pdf")
      , width = 9, height = 5)

  print("end of... plot_qc_summary_stats_violin_boxplots()")
}


plot_metadata_by_clinical_dx_numerical_boxplot <- function(
  seurat_obj, variable, ylim = NULL){

  print("plot_metadata_by_clinical_dx_numerical_boxplot()")

  variable_enquo <- enquo(variable)

  metadata_tb <- seurat_obj[[]] %>%
    select(clinical_dx, !!variable_enquo)

  if(is.null(ylim)){
    ylim <- max(metadata_tb[,variable], na.rm = TRUE) * 1.1
  }

  gg <- metadata_tb %>%
    ggplot(aes_string(x = "clinical_dx", y = variable, fill = "clinical_dx")) +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(y = c(0, ylim)) +
      theme(legend.position = "none") +
      # expand_limits(y = c(0, 5000)) +
      ggtitle(variable)

  print("end of... plot_metadata_by_clinical_dx_numerical_boxplot()")

  print(gg)
}

plot_metadata_by_clinical_dx_numerical_boxplot_wrapper <- function(
  seurat_obj = subset_so, out_graph_suffix = ""
  ){

  print("plot_metadata_by_clinical_dx_numerical_boxplot_wrapper()")

  plot_grid_wrapper(plotlist = list(
    plot_metadata_by_clinical_dx_numerical_boxplot(
      seurat_obj = seurat_obj, variable = "pmi_h"),
    plot_metadata_by_clinical_dx_numerical_boxplot(
      seurat_obj = seurat_obj, variable = "rin"),
    plot_metadata_by_clinical_dx_numerical_boxplot(
      seurat_obj = seurat_obj, variable = "weight_g"),
    plot_metadata_by_clinical_dx_numerical_boxplot(
      seurat_obj = seurat_obj, variable = "age"),
    plot_metadata_by_clinical_dx_numerical_boxplot(
      seurat_obj = seurat_obj, variable = "number_umi", ylim = 20000),
    plot_metadata_by_clinical_dx_numerical_boxplot(
      seurat_obj = seurat_obj, variable = "number_genes", ylim = 5000),
    plot_metadata_by_clinical_dx_numerical_boxplot(
      seurat_obj = seurat_obj, variable = "percent_mito")),
    title = make_plot_title("Numerical metadata by clinical dx"))
  ggsave(paste0(
      out_graph, "metadata_by_clinical_dx_boxplot", out_graph_suffix, ".pdf"),
    height = 10, width = 6)

  print("end of... plot_metadata_by_clinical_dx_numerical_boxplot_wrapper()")
}

plot_metadata_by_clinical_dx_factor_bar_plot <- function(
  seurat_obj, variable){

  print("plot_metadata_by_clinical_dx_factor_bar_plot()")

  variable_enquo <- enquo(variable)

  metadata_tb <- seurat_obj[[]] %>%
    select(clinical_dx, !!variable_enquo)

  gg <- metadata_tb %>%
    ggplot(aes_string(x = variable, fill = "clinical_dx")) +
      geom_bar(position = "dodge") +
      ylab("number of cells") +
      ggtitle(variable)

  print("end of... plot_metadata_by_clinical_dx_factor_bar_plot")

  print(gg)
}

plot_metadata_by_clinical_dx_factor_bar_plot_wrapper <- function(
  seurat_obj, out_graph_suffix = ""
  ){
print("plot_metadata_by_clinical_dx_factor_bar_plot_wrapper()")
   plot_grid_wrapper(plotlist = list(plot_metadata_by_clinical_dx_factor_bar_plot(seurat_obj = seurat_obj, variable = "sex")),
     title = make_plot_title("Number of cells for indicated metadata by clinical dx"), rel_height = 0.3)
   ggsave(paste0(
       out_graph, "metadata_by_clinical_dx_barplot", out_graph_suffix, ".pdf"),
     height = 7, width = 7)

  print("end of... plot_metadata_by_clinical_dx_factor_bar_plot_wrapper()")
}
################################################################################

### main function

main_function()
################################################################################
