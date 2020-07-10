# Damon Polioudakis
# 2019-05-23
# Function library

################################################################################

### asign cell type labels to cells

assign_cell_type_to_cell <- function(seurat_obj){

  print("assign_cell_type_to_cell")

  # marker_genes <- marker_genes_tb %>%
  #   filter(source == "tsai") %>%
  #   filter(gene_symbol %in% rownames(seurat_obj)) %>%
  #   split(x = ., f = .$marker_for) %>%
  #   map(., "gene_symbol")

  marker_genes <- marker_genes_refined_tb %>%
    filter(gene_symbol %in% rownames(seurat_obj)) %>%
    filter(! marker_for %in% c("ependymal")) %>%
    split(x = ., f = .$marker_for) %>%
    map(., "gene_symbol")

  # remove cell type module score columns if exist from previous use of function
  for (i in 1:length(names(marker_genes))){
    cell_type <- names(marker_genes)[i]
    seurat_obj_metadata_col_names <- colnames(seurat_obj[[]])
    idx <- grep(cell_type, seurat_obj_metadata_col_names)
    if(length(idx) > 0){
      seurat_obj[[seurat_obj_metadata_col_names[idx]]] <- NULL
    }
  }

  # AddModuleScore() appends a sequential digit to the end of the label
  # e.g. astrocyte to astrocyte1, endo to endo2
  seurat_obj <- AddModuleScore(object = seurat_obj, features = marker_genes,
    name = names(marker_genes))

  # remove digit added to end of labels by AddModuleScore()
  for (i in 1:length(names(marker_genes))){
    cell_type <- names(marker_genes)[i]
    seurat_obj_metadata_col_names <- colnames(seurat_obj[[]])
    idx <- grep(cell_type, seurat_obj_metadata_col_names)
    seurat_obj[[cell_type]] <- seurat_obj[[]][idx]
    seurat_obj[[seurat_obj_metadata_col_names[idx]]] <- NULL
  }

  seurat_obj$cell_type <- seurat_obj[[c("cell_ids", names(marker_genes))]] %>%
    as_tibble() %>%
    gather(key = "cell_type", value = "score", -cell_ids) %>%
    group_by(cell_ids) %>%
    dplyr::slice(which.max(score)) %>%
    dplyr::select(cell_ids, cell_type) %>%
    right_join(., seurat_obj[["cell_ids"]]) %>%
    pull(cell_type)

  return(seurat_obj)

}

assign_cell_type_cluster <- function(
  seurat_obj, cluster_col_name = "RNA_snn_res.0.6", percent_cutoff = 50){

  print("assign_cell_type_cluster")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]] %>% pull

  seurat_obj$cluster_cell_type <-
    seurat_obj[[c("cluster_ids", "cell_ids", "cell_type")]] %>%
      as_tibble() %>%
      group_by(cluster_ids) %>%
      count(cell_type) %>%
      mutate(percent = n/sum(n)*100) %>%
      dplyr::slice(which.max(percent)) %>%
      mutate(cluster_cell_type = ifelse(percent > percent_cutoff, cell_type, "mixed")) %>%
      right_join(., seurat_obj[[c("cluster_ids", "cell_ids", "cell_type")]],
        by = "cluster_ids") %>%
      pull(cluster_cell_type)

  return(seurat_obj)

}


assign_excitatory_layer_to_cell <- function(seurat_obj){

  print("assign_excitatory_layer_to_cell")

  marker_genes_exc_tb <- read_csv(in_exc_marker) %>%
    filter(! is.na(gene_symbol)) %>%
    mutate(groupings = marker_for) %>%
    filter(source == "lake")

  marker_genes <- bind_rows(
    marker_genes_exc_tb %>%
      slice(grep("L2.L3", .$marker_for)) %>%
      mutate(groupings = "L2.L3"),
    # marker_genes_exc_tb %>%
    #   slice(grep("L4", .$marker_for)) %>%
    #   mutate(groupings = "L4"),
    marker_genes_exc_tb %>%
      slice(grep("L4|L5|L6", .$marker_for)) %>%
      mutate(groupings = "L4.L5.L6")
      # ,
    # marker_genes_exc_tb %>%
    #   slice(grep("Gran", .$marker_for)) %>%
    #   mutate(groupings = "granule")
    ) %>%
    dplyr::select(gene_symbol, groupings) %>%
    split(x = ., f = .$groupings) %>%
    map(., "gene_symbol")

  # AddModuleScore() appends a sequential digit to the end of the label
  # e.g. astrocyte to astrocyte1, endo to endo2
  seurat_obj <- AddModuleScore(object = seurat_obj, features = marker_genes,
    name = names(marker_genes))

  # remove digit added to end of labels by AddModuleScore()
  for (i in 1:length(names(marker_genes))){
    cell_type <- names(marker_genes)[i]
    seurat_obj_metadata_col_names <- colnames(seurat_obj[[]])
    idx <- grep(cell_type, seurat_obj_metadata_col_names)
    seurat_obj[[cell_type]] <- seurat_obj[[]][idx]
    seurat_obj[[seurat_obj_metadata_col_names[idx]]] <- NULL
  }

  # assign layer to excitatory neurons
  layer_metadata_df <-
      seurat_obj[[c("cell_ids", "cell_type", names(marker_genes))]] %>%
        as_tibble() %>%
        gather(key = "layer", value = "score", -cell_ids, -cell_type) %>%
        group_by(cell_ids) %>%
        slice(which.max(score)) %>%
        # only assign layer to cells labeled as excitatory
        mutate(layer = ifelse(cell_type == "excitatory", layer, NA)) %>%
        # change L2.L3 to L2/L3, etc
        mutate(layer = gsub("\\.", "/", layer)) %>%
        mutate(cell_type_and_layer = ifelse(cell_type == "excitatory", layer, cell_type)) %>%
        dplyr::select(cell_ids, layer, cell_type_and_layer) %>%
        right_join(., seurat_obj[["cell_ids"]]) %>%
        as.data.frame()
  rownames(layer_metadata_df) <- layer_metadata_df$cell_ids
  seurat_obj <- AddMetaData(seurat_obj, layer_metadata_df)

  return(seurat_obj)

}

assign_cell_type_and_layer_to_cluster <- function(
  seurat_obj, cluster_col_name = "RNA_snn_res.0.6"){

  print("assign_cell_type_and_layer_to_cluster")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]] %>% pull

  if(! "cluster_cell_type" %in% colnames(seurat_obj[[]])){
    seurat_obj$cluster_cell_type <-
      seurat_obj[[c("cluster_ids", "cell_ids", "cell_type")]] %>%
        as_tibble() %>%
        group_by(cluster_ids) %>%
        count(cell_type) %>%
        mutate(percent = n/sum(n)*100) %>%
        slice(which.max(percent)) %>%
        mutate(cluster_cell_type = ifelse(percent > 50, cell_type, "mixed")) %>%
        right_join(., seurat_obj[[c("cluster_ids", "cell_ids", "cell_type")]],
          by = "cluster_ids") %>%
        pull(cluster_cell_type)
  }

  seurat_obj$cluster_cell_type_and_layer <-
    seurat_obj[[c("cluster_ids", "cell_ids", "cluster_cell_type", "layer")]] %>%
      as_tibble() %>%
      group_by(cluster_ids) %>%
      count(layer) %>%
      mutate(percent = n/sum(n)*100) %>%
      slice(which.max(percent)) %>%
      right_join(., seurat_obj[[c("cluster_ids", "cell_ids", "cluster_cell_type")]],
        by = "cluster_ids") %>%
      mutate(cluster_cell_type_and_layer = ifelse(percent > 30 && cluster_cell_type == "excitatory", layer, cluster_cell_type)) %>%
      pull(cluster_cell_type_and_layer)

  return(seurat_obj)

}
################################################################################

plot_dim_reduction_colored_by_variable <- function(
  dim_1, dim_2, variable_value, facet_variable = NULL
  , limit_high = NULL, limit_low = NULL
  , title = NULL, guide_size = 4, legend_title = NULL
  , alpha = 0.75, size = 0.3
  , expression_color_gradient = FALSE, zscore = FALSE){

    # browser()

  print("plot_dim_reduction_colored_by_variable")

  gg_tb <- tibble(dim_1 = dim_1, dim_2 = dim_2
    , value = variable_value)
  if (! is.null(facet_variable)){
    gg_tb <- mutate(gg_tb, facet_variable = facet_variable)
  }
  # Plot
  if (class(gg_tb$value) %in% c("character", "factor", "logical")){
    gg <- ggplot(gg_tb, aes(
        x = dim_1, y = dim_2, shape = value, col = value, group = value)) +
      geom_point(size = size, alpha = alpha, fill = NA) +
      scale_shape_manual(name = legend_title, values = rep(0:6,100)) +
      # option to facet
      { if (! is.null(facet_variable)){
        facet_wrap(~facet_variable, scales = "free")
      }}
    { if (length(unique(gg_tb$value)) > 18){
      gg <- gg + guides(color = guide_legend(title = legend_title
        , override.aes = list(size = guide_size, alpha = 1)
        , ncol = ceiling(length(unique(gg_tb$value))/18)))
      gg <- gg + guides(shape = guide_legend(title = legend_title
        , ncol = ceiling(length(unique(gg_tb$value))/18)))
    } else {
      gg <- gg + guides(color = guide_legend(title = legend_title
      , override.aes = list(size = guide_size, alpha = 1)))
    }}
  } else {
    if (is.null(limit_high)) {
      limit_high <- max(gg_tb$value)
    }
    if (is.null(limit_low)) {
      limit_low <- min(gg_tb$value)
    }
    gg <- ggplot(gg_tb, aes(x = dim_1, y = dim_2)) +
      geom_point(size = size, alpha = alpha, aes(col = value)) +
      # scale_shape_manual(values = rep(0:6,100)) +
      # option to facet
      { if (! is.null(facet_variable)){
        facet_wrap(~facet_variable, scales = "free")
      }} +
      # color scale options
      { if (expression_color_gradient == TRUE){
        scale_color_gradientn(name = legend_title
          , colours = c(
            "lightgrey", "#fee090", "#fdae61", "#f46d43", "#ca0020")
          )
      } else if (zscore == TRUE){
        scale_color_distiller(name = legend_title
            , type = "div", palette = 5, direction = -1
            , limits = c(limit_low, limit_high)
            , na.value = "grey90")
      } else {
        scale_color_viridis(
          name = legend_title, limits = c(limit_low, limit_high))
      }}
  }
  gg <- gg +
    ggplot_set_theme_publication +
    ggtitle(title)

  return(gg)
}

plot_dim_reduction_colored_by_expression <- function(
  genes
  , groupings = NULL
  , seurat_obj
  , assay_slot = "RNA"
  , expression_slot = "scale.data"
  , seurat_cluster_col
  , reduction = "tsne"
  , title = NULL
  , legend_title = "expression"
  , limit_high = 1.5
  , limit_low = -1.5
  , ncol = 4
  , zscore = FALSE
  , expression_color_gradient = FALSE){

  print("plot_dim_reduction_colored_by_expression")

  # function to collect tsne and expression values
  collect_dim_reduction_and_expression_values_from_seurat_obj <- function(genes){
    gg_tb <-
      # collect tsne values
      Embeddings(object = seurat_obj, reduction = reduction) %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble() %>%
      rename(dim_1 = 2, dim_2 = 3) %>%
      # add mean expression values to tibble
      inner_join(.
        # collect gene expression values and calculate mean
        , GetAssayData(object = seurat_obj[[assay_slot]], slot = expression_slot)[
          rownames(GetAssayData(object = seurat_obj[[assay_slot]], slot = expression_slot)) %in% genes, ,drop = FALSE] %>%
          t() %>%
          as.data.frame() %>%
          rownames_to_column("cell_ids") %>%
          as_tibble() %>%
          mutate(expression = rowMeans(.[,-1])) %>%
          mutate(expression = replace(
            expression, expression > limit_high, limit_high)) %>%
          mutate(expression = replace(
            expression, expression < limit_low, limit_low))
          # mutate(cell_ids = pull(seurat_obj[["cell_ids"]]))
        )
    return(gg_tb)
  }

  # collect tsne and expression values
  # if groupings argument exists then calculate mean expression for each group
  if (! is.null(groupings)){
    tsne_expression_ltb <- lapply(split(genes, groupings), function(genes){
      collect_dim_reduction_and_expression_values_from_seurat_obj(genes = genes)
      })
    # order by groupings argument, with any missing groupings removed
    list_order <- unique(groupings)[
      unique(groupings) %in% names(tsne_expression_ltb)]
    tsne_expression_ltb <- tsne_expression_ltb[list_order]
  # else collect expression for each gene individually
  } else {
    tsne_expression_ltb <- lapply(genes, function(gene){
    collect_dim_reduction_and_expression_values_from_seurat_obj(genes = gene)
    })
    names(tsne_expression_ltb) <- genes
  }

  # plot
  expression_gg_l <- lapply(names(tsne_expression_ltb), function(grouping){
    tsne_expression_tb <- tsne_expression_ltb[[grouping]]
    plot_dim_reduction_colored_by_variable(
      dim_1 = tsne_expression_tb$dim_1
      , dim_2 = tsne_expression_tb$dim_2
      , variable_value = tsne_expression_tb$expression
      , title = grouping
      , legend_title = legend_title
      , alpha = 0.25
      , expression_color_gradient = expression_color_gradient
      , zscore = zscore
      , limit_high = limit_high
      , limit_low = limit_low
    )
  })

  # tsne colored by cluster
  cluster_tb <- # collect tsne values
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    mutate(cluster = seurat_obj[[seurat_cluster_col]] %>% pull)
  cluster_gg <- plot_dim_reduction_colored_by_variable(
    dim_1 = cluster_tb$dim_1
    , dim_2 = cluster_tb$dim_2
    , variable_value = cluster_tb$cluster
    # , title = grouping
    , legend_title = "cluster"
    , alpha = 0.25
    , guide_size = 2
  )

  # combine graphs and plot
  # extract legend - make sure it exists (NA genes plot with no legend)
  plot_legend <- lapply(expression_gg_l, function(x) tryCatch(
    get_legend(x), error = function(e) NA)) %>%
      .[! is.na(.)] %>% .[[1]]
  # remove legends
  expression_gg_l <- lapply(expression_gg_l, function(gg){
    gg + theme(legend.position = "none")})
  # remove cluster legends if there are more than 10
  if (length(unique(cluster_tb$cluster)) > 10){
    cluster_gg <- cluster_gg + theme(legend.position = "none")}
  # plot
  plot_grid_wrapper(append(list(cluster_gg), expression_gg_l)
  , align = 'v', axis = 'r', ncol = ncol, rel_height = 0.2
  , title = title) %>%
    # add legend
    plot_grid(., plot_legend, rel_widths = c(1, .1))
}

plot_mean_expression_genes_heatmap <- function(
  seurat_obj, genes_to_plot, genes_groupings = NULL
  , cluster_col, subclust_keep = NULL, expression_slot = "data"
  , plot_title = NULL, limit_low = 0, limit_high = 5, zscore = FALSE){

    # browser()

  print("plot_mean_expression_marker_genes_heatmap")

  # set groups of genes to plot
  if(! is.null(genes_groupings)){
    gene_group_tb <- tibble(gene = genes_to_plot, group = genes_groupings) %>%
      mutate(gene_group = paste0(gene, "_", group))
  } else {
    gene_group_tb <- tibble(gene = genes_to_plot, group = genes_to_plot) %>%
      mutate(gene_group = paste0(gene, "_", group))
  }

  # cell id cluster id key
  cluster_ids <- factor(seurat_obj@meta.data[[cluster_col]])
  names(cluster_ids) <- rownames(seurat_obj@meta.data)
  Idents(seurat_obj) <- cluster_ids
  cellid_clusterid_tb <- Idents(seurat_obj) %>%
    enframe(name = "cell_id", value = "cluster")

  # expression z-scores
  mean_expr_zscores_tb <-
    # get expression data and subset to genes of interest
    GetAssayData(seurat_obj, slot = expression_slot)[
      rownames(GetAssayData(seurat_obj, slot = "data")) %in%
        gene_group_tb$gene, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    gather("cell_id", "expression", -gene) %>%
    # add gene group info
    right_join(., gene_group_tb, by = c("gene" = "gene")) %>%
    # add sub-clusters
    left_join(., cellid_clusterid_tb) %>%
    # filter(subcluster %in% subclust_keep) %>%
    # mean expression by subcluster
    group_by(cluster, group, gene, gene_group) %>%
    summarise(mean_expression = mean(expression)) %>%
    ungroup()

  # plot
  gg <- mean_expr_zscores_tb %>%
    # set limits
    mutate(mean_expression = if_else(mean_expression < limit_low, limit_low
      , if_else(mean_expression > limit_high, limit_high, mean_expression)
    )) %>%
    # select genes to plot
    # filter(group %in% genes_to_plot) %>%
    # set order of genes to plot
    mutate(gene_group = factor(gene_group
      , levels = rev(gene_group_tb$gene_group))) %>%
    # plot
    ggplot(aes(x = cluster, y = gene_group, fill = mean_expression)) +
      {if(! is.null(genes_groupings)) {
        facet_grid(rows = vars(group), scales = "free") }} +
      # facet_grid(rows = vars(group), scales = "free") +
      geom_tile() +
      # scale_fill_distiller(name = "Normalized\nexpression\nzscore"
      #     , type = "div", palette = 5, direction = -1
      #     , limits = c(limit_low, limit_high)
      #     , na.value = "grey90") +
          # color scale options
          { if (zscore == TRUE){
            scale_fill_distiller(name = "Normalized\nexpression\nzscore"
                , type = "div", palette = 5, direction = -1
                , limits = c(limit_low, limit_high)
                , na.value = "grey90")
          } else {
            scale_fill_viridis(name = "Normalized\nexpression"
                , limits = c(limit_low, limit_high)
                , na.value = "grey90")
          }} +

      ggplot_set_theme_publication +
      ggtitle(plot_title)
    return(gg)
}
################################################################################

### plot QC metrics (histograms, density plots)

plot_qc_metrics <- function(seurat_obj, in_seurat_raw){

  print("plot_qc_metrics")

  load(in_seurat_raw)
  nd_raw_so$cell_ids <- nd_raw_so[[]] %>% rownames


  # Plot gene detection curves

  nd_raw_so[[c("number_genes", "library_id")]] %>%
    as_tibble() %>%
    group_by(library_id) %>%
    mutate(rank = order(order(number_genes, decreasing = TRUE))) %>%
    ggplot(aes(x = rank, y = number_genes)) +
      geom_line() +
      facet_wrap(~library_id, scales = "free", ncol = 8) +
      geom_hline(yintercept = 200, color = "red") +
      coord_cartesian(ylim = c(0,4000)) +
      ggtitle(make_plot_title(paste0("Number of genes detected for each cell",
        "\nRed line indicates 200 genes detected")))
      ggsave(paste0(out_graph, "genes_detection_curve.pdf"),
        height = 22, width = 16)

  # Plot percent MT curves

  nd_raw_so[[c("percent_mito", "library_id")]] %>%
    as_tibble() %>%
    group_by(library_id) %>%
    mutate(rank = order(order(percent_mito, decreasing = FALSE))) %>%
    ggplot(aes(x = rank, y = percent_mito)) +
      geom_line() +
      facet_wrap(~library_id, scales = "free", ncol = 8) +
      geom_hline(yintercept = 5, color = "red") +
      coord_cartesian(ylim = c(0,100)) +
      ggtitle(make_plot_title(paste0("Percent MT for each cell",
        "\nRed line indicates 5% MT")))
      ggsave(paste0(out_graph, "percent_mt_curve.pdf"),
        height = 22, width = 16)

  # Plot percent MT by clinical_dx and region
  # raw data
  nd_raw_so[[c("percent_mito", "library_id", "clinical_dx", "region")]] %>%
    as_tibble() %>%
    mutate(clinical_dx = factor(clinical_dx,
      levels = c("Control", "AD", "PSP-S", "bvFTD"))) %>%
    ggplot(aes(x = clinical_dx, y = percent_mito, fill = clinical_dx)) +
      geom_boxplot() +
      facet_wrap(~region, scales = "free", ncol = 2) +
      coord_cartesian(ylim = c(0,25)) +
      ggtitle(make_plot_title(paste0(
        "Percent MT for each cell by clinical_dx and region",
        "\nRaw data")))
      ggsave(paste0(out_graph, "percent_mt_by_clinical_dx_raw_data_boxplot.pdf"),
        height = 6, width = 8)
  # QC filtered data
  seurat_obj[[c("percent_mito", "library_id", "clinical_dx", "region")]] %>%
    as_tibble() %>%
    mutate(clinical_dx = factor(clinical_dx,
      levels = c("Control", "AD", "PSP-S", "bvFTD"))) %>%
    ggplot(aes(x = clinical_dx, y = percent_mito, fill = clinical_dx)) +
      geom_boxplot() +
      facet_wrap(~region, scales = "free", ncol = 2) +
      coord_cartesian(ylim = c(0,25)) +
      ggtitle(make_plot_title(paste0(
        "Percent MT for each cell by clinical_dx and region",
        "\nQC filtered data")))
      ggsave(paste0(out_graph, "percent_mt_by_clinical_dx_boxplot.pdf"),
        height = 6, width = 8)

  # QC filtered data
  seurat_obj <- assign_cell_type_to_cell(seurat_obj)
  seurat_obj[[c("percent_mito", "library_id", "clinical_dx", "region", "cell_type")]] %>%
    as_tibble() %>%
    mutate(clinical_dx = factor(clinical_dx,
      levels = c("Control", "AD", "PSP-S", "bvFTD"))) %>%
    ggplot(aes(x = cell_type, y = percent_mito, fill = clinical_dx)) +
      geom_boxplot() +
      facet_wrap(~region, scales = "free", ncol = 2) +
      coord_cartesian(ylim = c(0,25)) +
      ggtitle(make_plot_title(paste0(
        "Percent MT for each cell by clinical_dx and region",
        "\nQC filtered data")))
      ggsave(paste0(out_graph, "percent_mt_by_clinical_dx_cell_type_boxplot.pdf"),
        height = 6, width = 14)

  # raw data
  nd_raw_so <- assign_cell_type_to_cell(nd_raw_so)
  nd_raw_so[[c("percent_mito", "library_id", "clinical_dx", "region", "cell_type")]] %>%
    as_tibble() %>%
    mutate(clinical_dx = factor(clinical_dx,
      levels = c("Control", "AD", "PSP-S", "bvFTD"))) %>%
    ggplot(aes(x = cell_type, y = percent_mito, fill = clinical_dx)) +
      geom_boxplot() +
      facet_wrap(~region, scales = "free", ncol = 2) +
      coord_cartesian(ylim = c(0,25)) +
      ggtitle(make_plot_title(paste0(
        "Percent MT for each cell by clinical_dx and region",
        "\nRaw data")))
      ggsave(paste0(out_graph, "percent_mt_by_clinical_dx_cell_type_raw_data_boxplot.pdf"),
        height = 6, width = 14)


  # Table of QC stats

  metadata <- c("cell_ids", "library_id", "number_genes", "number_umi", "percent_mito"
    )
    raw_tb <- nd_raw_so[[metadata]] %>%
      # rename(. ,replace = c("sample_number" = "sample_number")) %>%
      mutate(data_type = "raw data") %>%
      gather(variable, value, -data_type, -cell_ids, -library_id) %>%
      group_by(., library_id, variable) %>%
      summarize(.,
        mean = round(mean(value, na.rm = TRUE), 1),
        stdev = round(sd(value, na.rm = TRUE), 1)
      ) %>%
      ungroup %>%
      gather(key = "stat", value = "value",  -library_id, -variable) %>%
      unite("statistic", c("stat", "variable")) %>%
      spread(key = statistic, value = c("value"))

  filtered_tb <- seurat_obj[[c(metadata, "raw_cell_counts_per_lib", "cell_counts_per_lib")]] %>%
      mutate(data_type = "filtered data") %>%
      gather(variable, value, -data_type, -cell_ids, -library_id) %>%
      group_by(., library_id, variable) %>%
      summarize(.,
        mean = round(mean(value, na.rm = TRUE), 1),
        stdev = round(sd(value, na.rm = TRUE), 1)
      ) %>%
      ungroup %>%
      gather(key = "stat", value = "value",  -library_id, -variable) %>%
      unite("statistic", c("stat", "variable")) %>%
      spread(key = statistic, value = c("value"))

  statistics_tb <- left_join(raw_tb, filtered_tb, by = "library_id", suffix = c("_raw", "_filtered")) %>%
  left_join(., seurat_obj[[c("library_id", "clinical_dx", "rin", "pmi_h", "sscdnaquant_ng/ul", "concentration_qbit_ng/ul",  "Mean Reads per Cell", "Reads Mapped Confidently to Genome", "Reads Mapped Confidently to Transcriptome", "Reads Mapped Antisense to Gene")]] %>% distinct)

  write_csv(statistics_tb, path = paste0(out_table, "library_stats.csv"))


  ## Plots

  metadata <- c("cell_ids", "number_genes", "number_umi", "percent_mito"
    , "library_id")
    # , "mean_gc_content", "mean_cds_length")
  dat_tb <- nd_raw_so[[metadata]] %>%
    # rename(. ,replace = c("sample_number" = "sample_number")) %>%
    mutate(data_type = "raw data") %>%
    bind_rows(.
      , seurat_obj[[metadata]] %>%
        mutate(data_type = "filtered data")
      ) %>%
    mutate(data_type = factor(
      data_type, levels = c("raw data", "filtered data"))) %>%
    as_tibble() %>%
    mutate(mean_gc_content = 1) %>%
    mutate(mean_cds_length = 1)
    # gather(key = metric, value = value, -data_type, -cell_ids)

  # violin boxplots
  means_tb <- dat_tb %>%
  gather(variable, value, -data_type, -cell_ids, -library_id) %>%
    group_by(data_type, variable) %>%
    summarize(
      value_mean = round(mean(value, na.rm = TRUE), 1)
      , value_median = round(median(value, na.rm = TRUE),1)
    ) %>%
    mutate(label = paste0(
      "mean:\n", value_mean, "\nmedian:\n", value_median))
  dat_tb %>%
    gather(variable, value, -data_type, -cell_ids, -library_id) %>%
    group_by(data_type, variable) %>%
    ggplot(aes(x = data_type, y = value, fill = data_type)) +
      facet_wrap(~variable, scales = "free", ncol = 5) +
      geom_violin(fill = "lightgrey", color = "lightgrey") +
      geom_boxplot(width = 0.2, outlier.size = 0.1) +
      theme(
        legend.position = "none"
        , axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data = means_tb, aes(
        x = data_type, y = Inf, label = label), vjust = 1) +
      ggtitle(make_plot_title("QC metrics"))
    ggsave(paste0(out_graph, "qc_metrics_violin_boxplots.pdf")
      , width = 9, height = 5)

  plot_grid_wrapper(
    list(
      # number genes
      # histograms
      ggplot(dat_tb, aes(x = number_genes)) +
        geom_histogram(binwidth = 50) +
        facet_wrap(~data_type, ncol = 2) +
        ggtitle("Number genes (>0 counts). Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = number_genes, color = data_type)) +
        geom_density()
      # number umi
      # histograms
      , ggplot(dat_tb, aes(x = number_umi)) +
        geom_histogram(binwidth = 50) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = number_umi, color = data_type)) +
        geom_density() +
        ggtitle("Number UMI (>0 counts). Binwidth: 50")
      # percent MT
      # histograms
      , ggplot(dat_tb, aes(x = percent_mito)) +
        geom_histogram(binwidth = 1) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = percent_mito, color = data_type)) +
        geom_density() +
        ggtitle("Percentage MT. Binwidth: 1")

      # mean CDS length
      # histograms
      , ggplot(dat_tb, aes(x = mean_cds_length)) +
        geom_histogram(binwidth = 50) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = mean_cds_length, color = data_type)) +
        geom_density() +
        ggtitle("Mean CDS length. Binwidth: 50")
      # mean gc content
      # histograms
      , ggplot(dat_tb, aes(x = mean_gc_content)) +
        geom_histogram(binwidth = 1) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = mean_gc_content, color = data_type)) +
        geom_density() +
        ggtitle("Percentage GC. Binwidth: 1")
      )
    , rel_height = 0.1
    , title = (make_plot_title("Histograms and density plots of QC metrics"
      ))
  )
  ggsave(paste0(out_graph, "qc_metrics_histogram_density.pdf")
    , width = 13, height = 16)

  # number_genes
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = number_genes, fill = library_id)) +
        geom_histogram(binwidth = 50) +
        facet_grid(data_type~library_id) +
        coord_cartesian(xlim = c(0,6000)) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Number genes (>0 counts). Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = number_genes, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        coord_cartesian(xlim = c(0,6000)) +
          facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(paste0(
      "Histograms and density plots: number of genes detected per cell by sample"
      , "\nx axis limits 0-6,000"
    )))
  )
  ggsave(paste0(out_graph, "qc_number_genes_by_sample_histogram_density.pdf")
    , width = 22, height = 6)

  # number_umi
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = number_umi, fill = library_id)) +
        geom_histogram(binwidth = 50) +
        facet_grid(data_type~library_id) +
        coord_cartesian(xlim = c(0,10000)) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Number genes (>0 counts). Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = number_umi, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        coord_cartesian(xlim = c(0,10000)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(paste0(
      "Histograms and density plots: number of UMI per cell by sample"
      , "\nx axis limits 0-10,000"
    )))
  )
  ggsave(paste0(out_graph, "qc_number_umi_by_sample_histogram_density.pdf")
    , width = 22, height = 5)

  # percent_mito
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = percent_mito, fill = library_id)) +
        geom_histogram(binwidth = 1) +
        facet_grid(data_type~library_id) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Percentage MT. Binwidth: 1")
      # density plots
      , ggplot(dat_tb, aes(x = percent_mito, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        coord_cartesian(xlim = c(0,25)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(
      "Histograms and density plots: percent MT per cell by sample"
      ))
  )
  ggsave(paste0(out_graph, "qc_percent_mt_by_sample_histogram_density.pdf")
    , width = 22, height = 6)

  # mean cds length
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = mean_cds_length, fill = library_id)) +
        geom_histogram(binwidth = 50) +
        facet_grid(data_type~library_id) +
        # scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Mean CDS length. Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = mean_cds_length, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(
      "Histograms and density plots: mean CDS length per cell by sample"
      ))
  )
  ggsave(paste0(out_graph, "qc_mean_cds_length_by_sample_histogram_density.pdf")
    , width = 22, height = 5)

  # mean gc content
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = mean_gc_content, fill = library_id)) +
        geom_histogram(binwidth = 1) +
        facet_grid(data_type~library_id) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Mean GC content. Binwidth: 1")
      # density plots
      , ggplot(dat_tb, aes(x = mean_gc_content, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(
      "Histograms and density plots: mean gc content per cell by sample"
      ))
  )
  ggsave(paste0(out_graph, "qc_mean_gc_by_sample_histogram_density.pdf")
    , width = 22, height = 6)

}
################################################################################

plot_number_of_cells_filtered <- function(seurat_obj, in_seurat_raw){

  load(in_seurat_raw)

  plot_number_of_cells_filtered_stacked_bar_graph <- function(variable){

    print("plot_number_of_cells_filtered_stacked_bar_graph")

    # Cell counts per library id before filtering
    raw_cell_counts <- nd_raw_so[[variable]] %>%
      table %>%
      as_tibble %>%
      rename(raw_cell_counts = "n", variable = ".")

    # Cell counts per library id after filtering
    cells_remaining <- seurat_obj[[variable]] %>%
      table %>%
      as_tibble %>%
      rename(cells_remaining = "n", variable = ".")

    left_join(raw_cell_counts, cells_remaining, by = c("variable")) %>%
      mutate(cells_filtered = raw_cell_counts - cells_remaining) %>%
      distinct(variable, .keep_all = TRUE) %>%
      gather(key, value, -variable, -raw_cell_counts) %>%
      ggplot(aes(x = variable, y = value, fill = key)) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab(variable) +
        ylab("number of cells") +
        ggtitle(make_plot_title(
          paste0("Number of cells per ", variable, " before and after filtering")))
  }

  plot_number_of_cells_filtered_stacked_bar_graph("library_id")
  ggsave(paste0(out_graph, "number_cells_filtered_library_bar_graph.png")
    , width = 24, height = 7)
  plot_number_of_cells_filtered_stacked_bar_graph("clinical_dx")
  ggsave(paste0(out_graph, "number_cells_filtered_clinical_dx_bar_graph.png")
    , width = 7, height = 7)
  plot_number_of_cells_filtered_stacked_bar_graph("region")
  ggsave(paste0(out_graph, "number_cells_filtered_region_bar_graph.png")
    , width = 7, height = 7)

  variable <- c("clinical_dx", "region")

  # Cell counts per library id before filtering
  raw_cell_counts <- nd_raw_so[[variable]] %>%
    table %>%
    as_tibble %>%
    rename(raw_cell_counts = "n")
  # Cell counts per library id after filtering
  cells_remaining <- seurat_obj[[variable]] %>%
    table %>%
    as_tibble %>%
    rename(cells_remaining = "n")

  left_join(raw_cell_counts, cells_remaining, by = variable) %>%
    mutate(cells_filtered = raw_cell_counts - cells_remaining) %>%
    gather(key, value, -raw_cell_counts, -region, -clinical_dx) %>%
    ggplot(aes(x = region, y = value, fill = key)) +
      geom_bar(stat = "identity") +
      facet_grid(~clinical_dx) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab(variable) +
      ylab("number of cells") +
      ggtitle(make_plot_title(
        paste0("Number of cells per ", paste0(variable, collapse = " "), " before and after filtering")))
      ggsave(paste0(out_graph, "number_cells_filtered_region_clinical_dx_bar_graph.png")
        , width = 7, height = 7)
}
################################################################################

### PCA plots

plot_pca <- function(seurat_obj){

  print("plot_pca")

  # Plot genes with highest PC loadings
  pc_df <- Loadings(object = seurat_obj, reduction = "pca")[,1:8] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble %>%
    gather(., key = "pc", value = "pc_loading", -gene) %>%
    mutate(pc = gsub("V", "", pc)) %>%
    group_by(pc) %>%
    filter(! is.na(pc_loading)) %>%
    arrange(pc_loading) %>%
    filter(row_number() %in% c(1:15,(n() - 14):n())) %>%
    ungroup() %>%
    # filter(pc == "PC_1")
    as.data.frame()
  gg_l <- lapply(split(pc_df, pc_df$pc), function(df){
    df$gene <- factor(df$gene, levels = df$gene)
    ggplot(df, aes(x = pc_loading, y = gene)) +
      geom_point() +
      xlab("Loading") +
      ylab("Gene") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
      ggtitle(df$pc[1])
  })
  plot_grid_wrapper(gg_l, ncol = 8, align = 'v', axis = 'r', rel_height = 0.2
    , title = make_plot_title("Genes with highest PC loadings")
  )
  ggsave(paste0(out_graph, "pc_loading_plots.pdf"), width = 20, height = 7)

  # PC variance plot
  stdev <- seurat_obj@reductions$pca@stdev
  tibble(stdev = stdev, pc = c(1:length(stdev))) %>%
    ggplot(aes(x = pc, y = stdev)) +
      geom_point() +
      ggtitle(
        make_plot_title("PC variance")
        )
  ggsave(paste0(out_graph, "pc_variance_plot.pdf"), width = 5, height = 5)

}
################################################################################

### r squared matrix of PCs and metadata

plot_r_squared_matrix_of_pcs_and_metadata <- function(
  seurat_obj, pcs = 1:10, out_file_suffix = "",
  metadata_to_plot = c(
    "cluster_cell_type", "cell_ids",  "clinical_dx",
    "library_id", "pmi_h", "rin", "region", "sex", "weight_g", "age",
    "number_umi","number_genes", "percent_mito",
    # "concentration_qbit_ng/ul",
    "cell_counts_per_lib", "Fraction Reads in Cells", "Mean Reads per Cell",
    "Reads Mapped Confidently to Genome",
    "Reads Mapped Confidently to Transcriptome",
    "Reads Mapped Antisense to Gene")
  ){

  print("plot_r_squared_matrix_of_pcs_and_metadata")

  # make tibble of metadata and pc scores
  metadata_tb <-
    Embeddings(object = seurat_obj, reduction = "pca")[,pcs] %>%
      as.data.frame() %>%
      rownames_to_column(var = "cell_ids") %>%
      as_tibble() %>%
      # clean names
      set_names(., { gsub("V", "PC", colnames(.)) }) %>%
      # add metadata
      inner_join(.
        , seurat_obj[[metadata_to_plot]]
        ) %>%
      select(-cell_ids) %>%
      # clean and convert Cell Ranger percentages to numeric
      mutate_if(grepl("%", .[1,]), .funs = function(x){
        gsub("%", "", x) %>% as.numeric}
        ) %>%
      mutate_at(vars(contains("concentration_qbit_ng/ul")), as.numeric) %>%
      mutate_if(is.character, as.factor)
  # only keep factors with >= 2 levels
  metadata_tb <- metadata_tb %>%
    mutate_if(is.factor, droplevels)
  metadata_tb <- metadata_tb %>%
    select_if(function(col){
      length(levels(col)) > 1 |
      is.numeric(col)})
  #
  # metadata_tb %>%
  #   select_if(function(col){
  #     length(unique(col)) > 1
  #   })

  # calculate correlation of pcs and metadata
  r_squared_m <- matrix(nrow = ncol(metadata_tb), ncol = ncol(metadata_tb))
  colnames(r_squared_m) <- rownames(r_squared_m) <- colnames(metadata_tb)
  for(i in 1:nrow(r_squared_m)){
    for(j in 1:ncol(r_squared_m)){
      # print(rownames(r_squared_m)[i])
      # print(rownames(r_squared_m)[j])
      x <- metadata_tb %>% pull(rownames(r_squared_m)[i])
      y <- metadata_tb %>% pull(colnames(r_squared_m)[j])
      # tryCatch({
          if (class(x) == "numeric" & class(y) == "numeric") {
            x = scale(x, center = TRUE, scale = TRUE)
            y = scale(y, center = TRUE, scale = TRUE)
            lmout <- lm(y~x)
            r2 <- abs(summary(lmout)$adj.r.squared)
            # r2 <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
          } else {
            lmout <- lm(y~x)
            r2 <- abs(summary(lmout)$adj.r.squared)
          }
        # },
        #   error = function(cond){
        #     return(NA)
      # })
      r_squared_m[i,j] <- r2
    }
  }
  # fill in NA values from dependent variable being factor
  for (j in 1:ncol(r_squared_m)) {
    if (all(is.na(r_squared_m[ ,j])) == TRUE) {
      r_squared_m[ ,j] <- r_squared_m[j, ]
    }
  }

  # plot
  # reverse column and order for plotting
  r_squared_m <- r_squared_m[ ,ncol(r_squared_m):1]
  ggcorrplot(r_squared_m, type = "full", lab = TRUE, show.diag = TRUE
    , title = make_plot_title("|r^2| values for expression PCs and metadata")
  )
  ggsave(paste0(out_graph, "r_squared_matrix", out_file_suffix, ".pdf"),
    width = 16, height = 16)

  print("end of... plot_r_squared_matrix_of_pcs_and_metadata()")

}

plot_variance_explained_by_expression_pc <- function(
  seurat_obj, pcs, out_file_suffix = ""){

  print("plot_variance_explained_by_expression_pc()")

  # variance explained for expression PCs
  variance_by_pc <- Embeddings(object = seurat_obj, reduction = "pca") %>%
    apply(., 2, var)
  variance_explained_tb <- (variance_by_pc / sum(variance_by_pc) * 100) %>%
    as_tibble(rownames = "pc") %>%
    rename(percent_variance_explained = value) %>%
    mutate(pc = factor(pc, levels = pc)) %>%
    mutate(percent_variance_explained = round(percent_variance_explained, 2))

  ggplot(variance_explained_tb[pcs, ],
    aes(x = pc, y = percent_variance_explained, label = percent_variance_explained)) +
    geom_bar(stat = "identity") +
    geom_text() +
    ggtitle(make_plot_title(
      "Expression PCs percent variance explained
      PCA run on seurat variable genes (centered and scaled)"))
      ggsave(paste0(
          out_graph, "pc_variance_explained", out_file_suffix, ".pdf"),
        width = 6, height = 6)

  print("end of... plot_variance_explained_by_expression_pc()")
}
################################################################################

### Correlation of matrix of metadata

plot_correlation_matrix_of_metadata <- function(seurat_obj, in_seurat_raw){

  print("plot_correlation_matrix_of_metadata")

  load(in_seurat_raw)

  # make tibble of metadata
  metadata_tb <- nd_raw_so[[]] %>%
    as_tibble() %>%
    left_join(.,
      seurat_obj[[c("library_id", "cell_counts_per_lib")]] %>% distinct(),
      by = "library_id") %>%
    mutate(percent_cells_pass_qc = cell_counts_per_lib/raw_cell_counts_per_lib*100) %>%
    # clean and convert Cell Ranger percentages to numeric
    mutate_if(grepl("%", .[1,]), .funs = function(x){
      gsub("%", "", x) %>% as.numeric}
      ) %>%
    mutate_at("concentration_qbit_ng/ul", as.numeric) %>%
    select(-cell_ids, -orig.ident, -note, -notes, -cell_counts_per_lib, -volume_ul, -sampleposition, -side) %>%
    mutate_if(is.character, as.factor)

  # calculate correlation of metadata
  cor_m <- matrix(nrow = ncol(metadata_tb), ncol = ncol(metadata_tb))
  colnames(cor_m) <- rownames(cor_m) <- colnames(metadata_tb)
  for(i in 1:nrow(cor_m)){
    for(j in 1:ncol(cor_m)){
      x <- metadata_tb %>% pull(rownames(cor_m)[i])
      y <- metadata_tb %>% pull(colnames(cor_m)[j])
      if (class(x) == "numeric" & class(y) == "numeric") {
        r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
      } else {
        lmout <- lm(y~x)
        r <- sqrt(abs(summary(lmout)$adj.r.squared))
      }
      cor_m[i,j] <- r
    }
  }
  # fill in NA values from dependent variable being factor
  for (j in 1:ncol(cor_m)) {
    if (all(is.na(cor_m[ ,j])) == TRUE) {
      cor_m[ ,j] <- cor_m[j, ]
    }
  }

  # plot
  # reverse column and order for plotting
  cor_m <- cor_m[ ,ncol(cor_m):1]
  ggcorrplot(cor_m, type = "full", lab = TRUE, show.diag = TRUE
    , title = make_plot_title("|Spearman's rho| correlation values")
  )
  ggsave(paste0(out_graph, "correlation_matrix_metadata.png"), width = 20, height = 20)

}
################################################################################

### tSNE and clustering using different PCs

plot_tsne_clustering_pcs_test <- function(){

  print("plot_tsne_clustering_pcs_test")

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  plot_tsne_clustering_loom <- function(tsne_loom_path, seurat_cluster_col){
    tsne_gg <- filtered_loom[[tsne_loom_path]][,] %>% t() %>% as_tibble %>% rename(tsne1 = V1, tsne2 = V2) %>%
      mutate(cluster = filtered_loom[[seurat_cluster_col]][] %>% as.character) %>%
      ggplot(aes(x = tsne1, y = tsne2, color = cluster)) +
        geom_point(size = 0.05, alpha = 0.5) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        ggplot_set_theme_publication +
        ggtitle(paste0(tsne_loom_path, "\n", seurat_cluster_col))
    return(tsne_gg)
  }

  plot_tsne_clustering_loom("col_attrs/tsne_cell_embeddings_pc1to25", "col_attrs/cluster_pc1to25_res08")

  # organize plot_tsne_clustering_loom arguments with tibble
  # note cannot include () after plot_tsne_clustering_loom in pmap call
  tibble(
    tsne_loom_path = c(
      "col_attrs/tsne_cell_embeddings_pc1to25"
        , "col_attrs/tsne_cell_embeddings_pc1to50"
        , "col_attrs/tsne_cell_embeddings_pc1to75"
        , "col_attrs/tsne_cell_embeddings_pc1to100"
      )
    , seurat_cluster_col = c(
      "col_attrs/cluster_pc1to25_res08"
      , "col_attrs/cluster_pc1to50_res08"
      , "col_attrs/cluster_pc1to75_res08"
      , "col_attrs/cluster_ids"
    )) %>%
    # plot
    pmap(., plot_tsne_clustering_loom) %>%
    plot_grid_wrapper(plotlist = ., ncol = 2, rel_height = 0.2
      , align = 'v', axis = 'r'
      , title = paste0(script_name
        , "\n\nSeurat cluster and tSNE with different PCs used")
      )
    ggsave(paste0(out_graph, "pcs_test_tsne.png")
      , width = 12, height = 9)


  filtered_loom$close_all()
}
################################################################################

### Plot dim reduction colored by clustering with different Seurat resolutions

plot_dim_reduction_clustering_resolution_test <- function(
  seurat_obj, reduction = "tsne"){

  print("plot_dim_reduction_clustering_resolution_test")

  cluster_col_names <- c("RNA_snn_res.0.4", "RNA_snn_res.0.5"
    , "RNA_snn_res.0.6" , "RNA_snn_res.0.7", "RNA_snn_res.0.8")

  gg_l <- map(cluster_col_names, function(variable){
    gg_tb <-
      Embeddings(object = seurat_obj, reduction = reduction) %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble %>%
      rename(dim_1 = 2, dim_2 = 3) %>%
      inner_join(.
        , seurat_obj[[variable]] %>% rownames_to_column("cell_ids"))
    # gg_tb <- gg_l[[3]]
    plot_dim_reduction_colored_by_variable(
      dim_1 = gg_tb$dim_1
      , dim_2 = gg_tb$dim_2
      , variable_value = gg_tb[[variable]]
      , title = variable
      , legend_title = variable
      , size = 0.5
      , alpha = 0.5
      , guide_size = 2
    )
  })
  # set list names as metadata to use for file names
  names(gg_l) <- sapply(gg_l, function(gg){gg$labels$title})
  dir.create(paste0(out_graph, "resolution_test_", reduction))
  lapply(names(gg_l), function(name){
    gg_l[[name]] +
      # geom_point(size = 0.75, alpha = 0.25) +
      ggtitle(make_plot_title(paste0(
        "Cells colored by ", name
        , "\nDimensionality reduction: ", reduction)))
    ggsave(paste0(out_graph, "resolution_test_", reduction, "/", name, ".png")
      , width = 9, dpi = 400)
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.1
    , align = 'v', axis = 'r'
    , title = make_plot_title(
      paste0("Cells colored by cluster"
      , "\nDimensionality reduction: ", reduction))
    )
  ggsave(paste0(out_graph, "resolution_test_", reduction, ".png")
    , width = 28, height = 12, dpi = 200)

}
################################################################################

### dim red plot colored by metadata

plot_dim_reduction_colored_by_metadata <- function(
  seurat_obj, reduction = "tsne", metadata_to_plot = NULL, ncol = 4,
  plot_height = NULL, plot_width = 40){

  print("plot_dim_reduction_colored_by_metadata")

  if(is.null(metadata_to_plot)){
    metadata_to_plot <- c("cluster_ids"
      , "subclustering", "library_id"
      , "age", "clinical_dx",
      , "pmi_h",  "rin", "region", "sex", "weight_g",
      , "number_umi", "number_genes", "percent_mito", "concentration_qbit_ng/ul"
      , "cell_counts_per_lib", "Fraction Reads in Cells", "Mean Reads per Cell"
      , "Reads Mapped Confidently to Genome"
      , "Reads Mapped Confidently to Transcriptome", "Reads Mapped Antisense to Gene")
  }
    # , "mean_gc_content", "mean_cds_length")
  # make sure metadata to plot is in seurat object
  metadata_to_plot <- metadata_to_plot[
    metadata_to_plot %in% colnames(seurat_obj[[]])]

  gg_l <- map(metadata_to_plot, function(variable){
    gg_tb <-
      Embeddings(object = seurat_obj, reduction = reduction) %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble %>%
      rename(dim_1 = 2, dim_2 = 3) %>%
      inner_join(.
        , seurat_obj[[variable]] %>% rownames_to_column("cell_ids")) %>%
      # clean and convert Cell Ranger percentages to numeric
      mutate_if(grepl("%", .[1,]), .funs = function(x){
        gsub("%", "", x) %>% as.numeric}
        ) %>%
      mutate_at(vars(contains("concentration_qbit_ng/ul")), as.numeric)
    plot_dim_reduction_colored_by_variable(
      dim_1 = gg_tb$dim_1
      , dim_2 = gg_tb$dim_2
      , variable_value = gg_tb[[variable]]
      , title = variable
      , legend_title = variable
      , size = 0.5
      , alpha = 0.5
      , guide_size = 2
    )
  })
  # set list names as metadata to use for file names
  names(gg_l) <- sapply(gg_l, function(gg){gg$labels$title})
  dir.create(paste0(out_graph, "metadata_", reduction))
  lapply(names(gg_l), function(name){
    gg_l[[name]] +
      # geom_point(size = 0.75, alpha = 0.25) +
      ggtitle(make_plot_title(paste0(
        "Cells colored by ", name
        , "\nDimensionality reduction: ", reduction)))
    # remove "/"s from metadata labels for use in file names
    name <- gsub("/", "", name)
    ggsave(paste0(out_graph, "metadata_", reduction, "/", name, ".png")
      , width = 9, dpi = 400)
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = ncol, rel_height = 0.1
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0(
      "Cells colored by metadata"
      , "\nDimensionality reduction: ", reduction))
    )
  if (is.null(plot_height)){
    plot_height <- length(metadata_to_plot)/6+4
  }
  ggsave(paste0(out_graph, "metadata_", reduction, ".png"),
    width = plot_width, height = plot_height,
    dpi = 200, limitsize = FALSE)
    # , width = 6, height = 6, dpi = 200)
}
################################################################################

### dimensionality reduction colored by cluster

plot_dim_reduction_colored_by_individual_variable <- function(
  seurat_obj, var_name = "cluster_ids", reduction = "tsne"
  , factor_to_numeric = FALSE){

  print("plot_dim_reduction_colored_by_individual_variable")

  seurat_obj$variable_id <- seurat_obj[[var_name]]

  var_dim_reduction_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[["variable_id"]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble() %>%
      { if(factor_to_numeric == TRUE){
        mutate(., variable_id = variable_id %>%
          as.character() %>%
          as.numeric())
      } else .
      }
    )
  variables <- var_dim_reduction_tb$variable_id %>%
    unique() %>% sort()
  gg_l <- map(variables, function(variable){
    var_dim_reduction_tb %>%
      mutate(variable_membership = if_else(
        variable_id == variable, TRUE, FALSE)) %>%
      # sort so that TRUE points are printed on top of FALSE
      arrange(variable_membership) %>%
        # pull(variable_membership) %>% table
      {plot_dim_reduction_colored_by_variable(
        dim_1 = .$dim_1
        , dim_2 = .$dim_2
        , variable_value = .$variable_membership
        , title = paste0(var_name, ": ", variable)
        , legend_title = paste0(var_name, ": ", variable)
        , size = 0.1
        , alpha = 0.5
        , guide_size = 2
      ) + scale_color_manual(values = c("grey", "#00BFC4"))}
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.1
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0("Colored by ", var_name, "\n", reduction))
    )
  ggsave(paste0(out_graph, "variable_facet_", var_name, "_", reduction, ".png")
    , width = 26, height = 4.25+length(variables), dpi = 200, limitsize = FALSE)
    # , width = 6, height = 6, dpi = 200)
}
################################################################################

### Stacked bar charts of medadata by cluster

plot_cluster_metrics <- function(
  seurat_obj, cluster_col_name = "cluster_ids", out_graph_file_type = ".png"){

  print("plot_cluster_metrics()")

  plot_number_of_cells_per_cluster_bar_plot(
    seurat_obj = seurat_obj,
    cluster_col_name = cluster_col_name,
    out_graph_file_type = out_graph_file_type
  )

  plot_metadata_by_cluster_stacked_bar_plot(
    seurat_obj = seurat_obj,
    cluster_col_name = cluster_col_name,
    out_graph_file_type = out_graph_file_type
  )

  plot_metadata_by_cluster_percent_stacked_bar_plot(
    seurat_obj = seurat_obj,
    cluster_col_name = cluster_col_name,
    out_graph_file_type = out_graph_file_type
  )

  plot_metadata_by_cluster_percent_heatmap(
    seurat_obj = seurat_obj,
    cluster_col_name = cluster_col_name,
    out_graph_file_type = out_graph_file_type
  )

  output_percent_library_by_dx_cluster_table(
    seurat_obj = seurat_obj,
    cluster_col_name = cluster_col_name,
    out_graph_file_type = out_graph_file_type
  )

  print("end of... plot_cluster_metrics()")
}

gather_seurat_obj_metadata <- function(seurat_obj){

  print("gather_seurat_obj_metadata()")

  metadata_tb <- seurat_obj[[]] %>%
    as_tibble %>%
    select_if(negate(is.numeric)) %>%
    # make sure clusters plot 0 to highest number
    mutate(cluster_ids = factor(
      cluster_ids,
        levels = cluster_ids %>%
        as.character() %>%
        as.numeric() %>%
        unique() %>%
        sort
      ))

  print("end of... gather_seurat_obj_metadata()")
  return(metadata_tb)
}

plot_number_of_cells_per_cluster_bar_plot <- function(
  seurat_obj, cluster_col_name = "cluster_ids", out_graph_file_type = ".png"
  ){

  print("plot_number_of_cells_per_cluster_bar_plot()")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  metadata_tb <- gather_seurat_obj_metadata(seurat_obj = seurat_obj) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  print("plot_number_of_cells_per_cluster_bar_plot()")
  metadata_tb %>%
    # mutate(cluster_ids = as.numeric(cluster_ids)) %>%
    filter(variable == "library_id") %>%
    group_by(cluster_ids) %>%
    summarise(count = n()) %>%
    ggplot(aes(x = cluster_ids, y = count)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = count, x = cluster_ids, y = count)
        , vjust = 0) +
      ggtitle(make_plot_title("Number of cells per cluster"))
      ggsave(paste0(
        out_graph, "number_cells_per_cluster_bargraph", out_graph_file_type)
        , width = 11)

  print("end of... plot_number_of_cells_per_cluster_bar_plot")
}

plot_metadata_by_cluster_stacked_bar_plot <- function(
  seurat_obj, cluster_col_name = "cluster_ids", out_graph_file_type = ".png"
  ){

  print("plot_metadata_by_cluster_stacked_bar_plot()")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  metadata_tb <- gather_seurat_obj_metadata(seurat_obj = seurat_obj) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  # Count stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){
      gg <- .[[name]]
      ggplot(gg, aes(x = cluster_ids, fill = value)) +
        facet_wrap(~variable) +
        geom_bar() +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }}
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "r", ncol = 2, rel_height = 0.1
      , title = make_plot_title("Metadata by cluster"))
     # +
        # guide = guide_legend(nrow=2)) +
      ggsave(paste0(
        out_graph, "metadata_by_cluster_stackedbar", out_graph_file_type)
        , width = 24, height = 70, limitsize = FALSE)
}

plot_metadata_by_cluster_percent_stacked_bar_plot <- function(
  seurat_obj, cluster_col_name = "cluster_ids", out_graph_file_type = ".png"
  ){

  print("plot_metadata_by_cluster_percent_stacked_bar_plot()")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  metadata_tb <- gather_seurat_obj_metadata(seurat_obj = seurat_obj) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  # Percent stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){
      gg <- .[[name]]
      gg <- gg %>%
        group_by(cluster_ids) %>%
        count(value) %>%
        mutate(percent = n/sum(n)) %>%
        mutate(variable = name) %>%
        ungroup() %>%
        mutate(cluster_ids = as.numeric(as.character(cluster_ids))) %>%
        mutate(cluster_ids = factor(cluster_ids, levels = sort(unique(cluster_ids))))
      ggplot(gg, aes(x = cluster_ids, y = percent, fill = value)) +
        facet_wrap(~variable) +
        geom_bar(stat = "identity") +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }}
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "r", ncol = 2, rel_height = 0.1
      , title = make_plot_title("Metadata by cluster"))
      ggsave(paste0(
        out_graph, "metadata_by_cluster_percent_stackedbar", out_graph_file_type)
        , width = 24, height = 70, limitsize = FALSE)

  print("end of... plot_metadata_by_cluster_percent_stacked_bar_plot()")
}

plot_metadata_by_cluster_percent_heatmap <- function(
  seurat_obj, cluster_col_name = "cluster_ids", out_graph_file_type = ".png"
  ){

  print("plot_metadata_by_cluster_percent_heatmap()")

  metadata_tb <- gather_seurat_obj_metadata(seurat_obj = seurat_obj) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]
  # heatmap of percent of cells from each library for each cluster
  name <- "library_id"
  metadata_tb %>%
    filter(variable %in% name) %>%
    group_by(cluster_ids) %>%
    count(value) %>%
    mutate(percent = n/sum(n)*100) %>%
    mutate(variable = name) %>%
    ungroup() %>%
    complete(cluster_ids, value, variable
      , fill = list(n = 0, percent = 0)) %>%
    ggplot(aes(x = cluster_ids, y = value, fill = percent)) +
      geom_tile() +
      geom_text(aes(label = paste0(round(percent, 1), "%  (", n, ")"))) +
      scale_fill_gradient(
        low = "white", high = "red", na.value = "grey", limits = c(0,25)) +
      ylab(name) +
      ggplot_set_theme_publication +
      ggtitle(make_plot_title(
        paste0("Percent of cells from each library for each cluster"
          , "\nEach column sums to 100%, () indicates number of cells")))
      ggsave(paste0(out_graph, "metadata_by_cluster_percent_heatmap.png")
        , width = 30, height = 17)

  print("end of... plot_metadata_by_cluster_percent_heatmap()")
}

output_percent_library_by_dx_cluster_table <- function(
  seurat_obj, cluster_col_name = "cluster_ids"
  ){

  print("output_percent_library_by_dx_cluster_table()")

  metadata_tb <- gather_seurat_obj_metadata(seurat_obj = seurat_obj) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]
    # table of variance in percentage of cells for each cluster by disease by library
    metadata_tb %>%
      filter(variable %in% c("library_id", "clinical_dx")) %>%
      spread(key = variable, value = value) %>%
      group_by(cluster_ids) %>%
      count(library_id, clinical_dx) %>%
      mutate(percent = n/sum(n)*100) %>%
      group_by(cluster_ids, clinical_dx) %>%
      summarize(min(percent), max(percent), var(percent), sd(percent)) %>%
      write_csv(path = paste0(out_table, "percent_library_by_dx_cluster.csv"))

  print("end of... output_percent_library_by_dx_cluster_table()")

}
################################################################################

### stacked bar plot of metadata by percent of cells for each cluster

plot_metadata_by_cluster_percent_stacked_barplot <- function(
  seurat_obj, cluster_col_name = "cluster_ids", metadata_to_plot = NULL,
  legend_nrow = 18, plot_height = 70, plot_width = 24){

  print("plot_metadata_by_cluster_percent_stacked_barplot")

  if(is.null(metadata_to_plot)){
    metadata_to_plot <- colnames(seurat_obj[[]])
  }

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  metadata_tb <- seurat_obj[[c("cell_ids", "cluster_ids", "orig.ident", metadata_to_plot)]] %>%
    as_tibble %>%
    select_if(negate(is.numeric)) %>%
    # make sure clusters plot 0 to highest number
    mutate(cluster_ids = factor(
      cluster_ids
      , levels = cluster_ids %>%
        as.character() %>%
        as.numeric() %>%
        unique() %>%
        sort
      )) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  # Percent stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){

      gg <- .[[name]]
      gg <- gg %>%
        group_by(cluster_ids) %>%
        count(value) %>%
        mutate(percent = n/sum(n)) %>%
        mutate(variable = name) %>%
        ungroup() %>%
        mutate(cluster_ids = as.numeric(as.character(cluster_ids))) %>%
        mutate(cluster_ids = factor(cluster_ids, levels = sort(unique(cluster_ids))))

      ggplot(gg, aes(x = cluster_ids, y = percent, fill = value)) +
        facet_wrap(~variable) +
        geom_bar(stat = "identity") +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }} +
        guides(fill = guide_legend(ncol = ceiling(length(unique(gg$value))/legend_nrow)))
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "r", ncol = 2, rel_height = 0.1
      , title = make_plot_title("Metadata by cluster"))
      ggsave(paste0(out_graph, "metadata_by_cluster_percent_stackedbar.png")
        , width = plot_width, height = plot_height, limitsize = FALSE)

}

### stacked bar plot of metadata per generalized faceting variable

plot_metadata_by_variable_stacked_barplot <- function(
  seurat_obj, var_name = "cluster_ids", metadata_to_plot = NULL,
  legend_nrow = 18, plot_height = 70, plot_width = 24){

  print("plot_metadata_by_cluster_percent_stacked_barplot")

  if(is.null(metadata_to_plot)){
    metadata_to_plot <- colnames(seurat_obj[[]])
  }

  metadata_tb <- seurat_obj@meta.data

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  metadata_tb <- seurat_obj[[c("cell_ids", "cluster_ids", "orig.ident", metadata_to_plot)]] %>%
    as_tibble %>%
    select_if(negate(is.numeric)) %>%
    # make sure clusters plot 0 to highest number
    mutate(cluster_ids = factor(
      cluster_ids
      , levels = cluster_ids %>%
        as.character() %>%
        as.numeric() %>%
        unique() %>%
        sort
      )) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  # Percent stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){

      gg <- .[[name]]
      gg <- gg %>%
        group_by(cluster_ids) %>%
        count(value) %>%
        mutate(percent = n/sum(n)) %>%
        mutate(variable = name) %>%
        ungroup() %>%
        mutate(cluster_ids = as.numeric(as.character(cluster_ids))) %>%
        mutate(cluster_ids = factor(cluster_ids, levels = sort(unique(cluster_ids))))

      ggplot(gg, aes(x = cluster_ids, y = percent, fill = value)) +
        facet_wrap(~variable) +
        geom_bar(stat = "identity") +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }} +
        guides(fill = guide_legend(ncol = ceiling(length(unique(gg$value))/legend_nrow)))
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "r", ncol = 2, rel_height = 0.1
      , title = make_plot_title("Metadata by cluster"))
      ggsave(paste0(out_graph, "metadata_by_cluster_percent_stackedbar.png")
        , width = plot_width, height = plot_height, limitsize = FALSE)

}

################################################################################

### tSNE colored by PC score

plot_tsne_colored_by_pc_score <- function(seurat_obj){

  print("plot_tsne_colored_by_pc_score")

  gg_tb <- map(1:15, function(pc){
    pc_variable_name <- paste0("pc_", pc)
    gg_tb <-
      # collect tsne values
      Embeddings(object = seurat_obj, reduction = "tsne") %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble %>%
      rename(tsne1 = tSNE_1, tsne2 = tSNE_2) %>%
      # collect pca scores
      inner_join(.
        , Embeddings(object = seurat_obj, reduction = "pca")[,pc] %>%
          enframe() %>%
          as_tibble() %>%
          rename(cell_ids = name, pc_scores = value)
        )
    plot_dim_reduction_colored_by_variable(
      dim_1 = gg_tb$tsne1
      , dim_2 = gg_tb$tsne2
      , variable_value = gg_tb$pc_scores
      , title = pc_variable_name
      , legend_title = pc_variable_name
    )
  }) %>%
  plot_grid_wrapper(plotlist = ., ncol = 4, rel_height = 0.1
    , title = make_plot_title(
      "Aggregated P1 and P2 samples colored by PC score")
    )
  ggsave(paste0(out_graph, "pc_score_tsne.png")
    , width = 14, height = 12, dpi = 200)

}
################################################################################

### dim plot colored by expression of known marker genes

plot_marker_expression_dim_reduction <- function(
  seurat_obj, cluster_col_name = "cluster_ids", reduction = "tsne",
  plot_individually = TRUE){

  print("plot_marker_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  marker_genes_tb <- marker_genes_tb %>% filter(! is.na(gene_symbol))

  # append source to cell type names to distinguish groupings
  groupings <- with(marker_genes_tb
    , paste0(marker_for, "_", clean_strings(source)))

  marker_genes_tb$groupings <- groupings

  # plot mean of groups
  plot_dim_reduction_colored_by_expression(
    genes = marker_genes_tb$gene_symbol
    , groupings = groupings
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , reduction = reduction
    , expression_color_gradient = TRUE
    , ncol = 4
    , limit_high = 3
    , limit_low = 0
    , title = make_plot_title(
      "Mean normalized expression of cell type markers")
  )
  ggsave(paste0(
    out_graph, "markers_gene_expression_", reduction, ".png")
    , width = 18, height = 18, dpi = 200)

  # plot individually and output in sub directory
  if (plot_individually == TRUE){
    dir.create(paste0(out_graph, "markers_individually_gene_expression_", reduction))
    # set list names as metadata to use for file names
    marker_genes_ltb <- split(marker_genes_tb, marker_genes_tb$groupings)
    lapply(names(marker_genes_ltb), function(name){
      marker_genes_tb <- marker_genes_ltb[[name]]
      tryCatch(
        {
          plot_dim_reduction_colored_by_expression(
            genes = marker_genes_tb$gene_symbol
            , seurat_obj = seurat_obj
            , seurat_cluster_col = "cluster_ids"
            , expression_slot = "data"
            , reduction = reduction
            , expression_color_gradient = TRUE
            , ncol = 4
            , limit_high = 3
            , limit_low = 0
            , title = make_plot_title(paste0(name
              , "\nNormalized expression of cell type markers")))
        },
          error = function(cond){
            return(NA)
        })
      ggsave(paste0(
        out_graph, "markers_individually_gene_expression_", reduction, "/", name, ".png")
        , width = 18, height = 2+0.75*length(marker_genes_tb$gene_symbol)
        , dpi = 200, limitsize = FALSE)
    })
  }

}
################################################################################

### dim plot colored by expression of known marker genes

plot_marker_refined_expression_dim_reduction <- function(
  seurat_obj, cluster_col_name = "cluster_ids"
  , reduction = "tsne"){

  print("plot_marker_refined_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  marker_genes_tb <- marker_genes_refined_tb %>% filter(! is.na(gene_symbol))

  marker_genes_tb$groupings <- marker_genes_tb$marker_for

  # plot mean of groups
  plot_dim_reduction_colored_by_expression(
    genes = marker_genes_tb$gene_symbol
    , groupings = marker_genes_tb$groupings
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , reduction = reduction
    , expression_color_gradient = TRUE
    , ncol = 4
    , limit_high = 3
    , limit_low = 0
    , title = make_plot_title(
      "Mean normalized expression of cell type markers")
  )
  ggsave(paste0(
    out_graph, "markers_refined_gene_expression_", reduction, ".png")
    , width = 18, height = 14, dpi = 200)

  # plot individually and output in sub directory
  dir.create(paste0(out_graph, "markers_refined_individually_gene_expression_", reduction))
  # set list names as metadata to use for file names
  marker_genes_ltb <- split(marker_genes_tb, marker_genes_tb$groupings)
  lapply(names(marker_genes_ltb), function(name){
    marker_genes_tb <- marker_genes_ltb[[name]]
    tryCatch(
      {
        plot_dim_reduction_colored_by_expression(
          genes = marker_genes_tb$gene_symbol
          , seurat_obj = seurat_obj
          , seurat_cluster_col = "cluster_ids"
          , expression_slot = "data"
          , reduction = reduction
          , expression_color_gradient = TRUE
          , ncol = 4
          , limit_high = 3
          , limit_low = 0
          , title = make_plot_title(paste0(name
            , "\nNormalized expression of cell type markers")))
      },
        error = function(cond){
          return(NA)
      })
    ggsave(paste0(
      out_graph, "markers_refined_individually_gene_expression_", reduction, "/", name, ".png")
      , width = 18, height = 2+length(marker_genes_tb$gene_symbol)
      , dpi = 200, limitsize = FALSE)
  })

}
################################################################################

### dim plot colored by expression of known marker genes

plot_excitatory_marker_expression_dim_reduction <- function(
  seurat_obj, cluster_col_name = "cluster_ids", reduction = "umap",
  in_exc_marker = in_exc_marker){

  print("plot_excitatory_marker_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  marker_genes_tb <- read_csv(in_exc_marker)
  marker_genes_tb <- marker_genes_tb %>% filter(! is.na(gene_symbol))
  marker_genes_tb$groupings <- marker_genes_tb$marker_for
  marker_genes_tb <- marker_genes_tb %>% filter(source == "hodge")

  # plot mean of groups
  plot_dim_reduction_colored_by_expression(
    genes = marker_genes_tb$gene_symbol
    # , groupings = marker_genes_tb$groupings
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , reduction = reduction
    , expression_color_gradient = TRUE
    , ncol = 4
    , limit_high = 3
    , limit_low = 0
    , title = make_plot_title(
      "Mean normalized expression of cell type markers")
  )
  ggsave(paste0(
    out_graph, "markers_excitatory_hodge_gene_expression_", reduction, ".png")
    , width = 18, height = 30, dpi = 200)

  marker_genes_tb <- read_csv(in_exc_marker)
  marker_genes_tb <- marker_genes_tb %>% filter(! is.na(gene_symbol))
  marker_genes_tb$groupings <- marker_genes_tb$marker_for
  marker_genes_tb <- marker_genes_tb %>% filter(source == "lake")

  # plot mean of groups
  plot_dim_reduction_colored_by_expression(
    genes = marker_genes_tb$gene_symbol
    # , groupings = marker_genes_tb$groupings
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , reduction = reduction
    , expression_color_gradient = TRUE
    , ncol = 4
    , limit_high = 3
    , limit_low = 0
    , title = make_plot_title(
      "Mean normalized expression of cell type markers")
  )
  ggsave(paste0(
    out_graph, "markers_excitatory_lake_gene_expression_", reduction, ".png")
    , width = 18, height = 30, dpi = 200)

}
################################################################################

### dimensionality reduction plot colored by expression of astrocyte markers

plot_marker_astrocyte_expression_dim_reduction <- function(
  seurat_obj, cluster_col_name = "cluster_ids", reduction = "tsne"){

  print("plot_marker_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # clean
  astro_markers_tb <- astro_markers_tb %>%
    filter(! is.na(gene_symbol)) %>%
    mutate(gene = gsub("\\*", "", .$gene_symbol)) %>% as.data.frame()

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = astro_markers_tb$gene
    , groupings = astro_markers_tb$marker_for
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , expression_color_gradient = TRUE
    , reduction = reduction
    , ncol = 4
    , limit_high = 2
    , limit_low = 0
    , title = make_plot_title("Mean normalized expression of astrocyte markers")
  )
  ggsave(paste0(out_graph, "markers_astrocytes_expression_", reduction, ".png")
    , width = 18, height = 8, dpi = 200)

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = astro_markers_tb$gene
    , groupings = astro_markers_tb$marker_for
    , seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , seurat_cluster_col = "cluster_ids"
    , reduction = reduction
    , legend_title = "expression\nz-score"
    , zscore = TRUE
    , ncol = 4
    , limit_high = 1
    , limit_low = -1
    , title = make_plot_title("Mean normalized expression of astrocyte markers")
  )
  ggsave(paste0(
    out_graph, "markers_astrocytes_expression_zscore_", reduction, ".png")
    , width = 18, height = 8, dpi = 200)

  # plot individually and output in sub directory
  dir.create(paste0(
    out_graph, "markers_astrocytes_individually_expression_", reduction))
  # set list names as metadata to use for file names
  marker_genes_ltb <- split(astro_markers_tb, astro_markers_tb$marker_for)
  lapply(names(marker_genes_ltb), function(name){
    astro_markers_tb <- marker_genes_ltb[[name]]
    tryCatch(
      {
    plot_dim_reduction_colored_by_expression(
      genes = astro_markers_tb$gene
      , seurat_obj = seurat_obj
      , seurat_cluster_col = "cluster_ids"
      , expression_slot = "data"
      , reduction = reduction
      , expression_color_gradient = TRUE
      , ncol = 4
      , limit_high = 2
      , limit_low = 0
      , title = make_plot_title(paste0(name
        , "\nNormalized expression of astrocyte markers")))
      },
        error = function(cond){
          return(NA)
      })
    ggsave(paste0(
      out_graph, "markers_astrocytes_individually_expression_", reduction, "/"
      , name, ".png")
      , width = 18, height = 2+0.75*length(astro_markers_tb$gene_symbol)
      , dpi = 200, limitsize = FALSE)
  })

  # plot individually and output in sub directory
  dir.create(paste0(
    out_graph, "markers_astrocytes_individually_expression_zscore_", reduction))
  # set list names as metadata to use for file names
  marker_genes_ltb <- split(astro_markers_tb, astro_markers_tb$marker_for)
  lapply(names(marker_genes_ltb), function(name){
    astro_markers_tb <- marker_genes_ltb[[name]]
    tryCatch(
      {
    plot_dim_reduction_colored_by_expression(
      genes = astro_markers_tb$gene
      , seurat_obj = seurat_obj
      , expression_slot = "scale.data"
      , seurat_cluster_col = "cluster_ids"
      , legend_title = "expression\nz-score"
      , reduction = reduction
      , zscore = TRUE
      , ncol = 4
      , limit_high = 1
      , limit_low = -1
      , title = make_plot_title(paste0(name
        , "\nNormalized expression of astrocyte markers")))
      },
        error = function(cond){
          return(NA)
      })
    ggsave(paste0(
      out_graph, "markers_astrocytes_individually_expression_zscore_", reduction, "/"
      , name, ".png")
      , width = 18, height = 2+0.75*length(astro_markers_tb$gene_symbol)
      , dpi = 200, limitsize = FALSE)
  })

}
################################################################################

plot_dim_reduction_colored_by_module_score <- function(
  seurat_obj, reduction = "umap"){

  print("plot_dim_reduction_colored_by_module_score")

  cell_marker_scores_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble() %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[[c("RNA_snn_res.0.6", "astrocyte1", "endo2", "excitatory3", "inhibitory4", "micro5", "oligo6", "OPC7")]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble()
    )

    gg_l <- map(c("RNA_snn_res.0.6", "astrocyte1", "endo2", "excitatory3", "inhibitory4", "micro5", "oligo6", "OPC7"), function(variable){
      gg_tb <- cell_marker_scores_tb
      plot_dim_reduction_colored_by_variable(
        dim_1 = gg_tb$dim_1
        , dim_2 = gg_tb$dim_2
        , variable_value = gg_tb[[variable]]
        , title = variable
        , legend_title = variable
        , expression_color_gradient = TRUE
        , size = 0.5
        , alpha = 0.5
        , guide_size = 2
      )
    })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.3
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0("Colored by marker expression score", "\n", reduction))
    )
  ggsave(paste0(out_graph, "marker_score_", reduction, ".png")
    , width = 20, height = 8, dpi = 200)

}

plot_dim_reduction_colored_by_individual_variable <- function(
  seurat_obj, var_name = "cell_type_and_layer", reduction = "umap",
  factor_to_numeric = FALSE){

  print("plot_dim_reduction_colored_by_individual_variable")

  seurat_obj$variable_id <- seurat_obj[[var_name]]

  var_dim_reduction_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[["variable_id"]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble() %>%
      { if(factor_to_numeric == TRUE){
        mutate(., variable_id = variable_id %>%
          as.character() %>%
          as.numeric())
      } else .
      }
    )
  variables <- var_dim_reduction_tb$variable_id %>%
    unique() %>% sort()
  gg_l <- map(variables, function(variable){
    var_dim_reduction_tb %>%
      mutate(variable_membership = if_else(
        variable_id == variable, TRUE, FALSE)) %>%
      # sort so that TRUE points are printed on top of FALSE
      arrange(variable_membership) %>%
        # pull(variable_membership) %>% table
      {plot_dim_reduction_colored_by_variable(
        dim_1 = .$dim_1
        , dim_2 = .$dim_2
        , variable_value = .$variable_membership
        , title = paste0(var_name, ": ", variable)
        , legend_title = paste0(var_name, ": ", variable)
        , size = 0.1
        , alpha = 0.5
        , guide_size = 2
      ) + scale_color_manual(values = c("grey", "#00BFC4"))}
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.3
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0("Colored by ", var_name, "\n", reduction))
    )
  ggsave(paste0(out_graph, "variable_facet_", var_name, "_", reduction, ".png")
    , width = 18, height = 3.25+0.6*length(variables), dpi = 200, limitsize = FALSE)
}

plot_cell_type_by_cluster_stacked_bar_plot <- function(
  seurat_obj, cluster_col_name = "RNA_snn_res.0.6"){

  print("plot_cell_type_by_cluster_stacked_bar_plot")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # collect metadata
  metadata_tb <- seurat_obj[[]] %>%
    as_tibble %>%
    select_if(negate(is.numeric)) %>%
    # make sure clusters plot 0 to highest number
    mutate(cluster_ids = factor(
      cluster_ids
      , levels = cluster_ids %>%
        as.character() %>%
        as.numeric() %>%
        unique() %>%
        sort
      )) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

    # percent stacked bar chart
    # calculate percent
    gg_tb <- metadata_tb %>%
      filter(variable == "cell_type") %>%
      group_by(cluster_ids) %>%
      count(value) %>%
      mutate(percent = n/sum(n)) %>%
      mutate(variable = "cell_type") %>%
      ungroup() %>%
      mutate(cluster_ids = as.numeric(as.character(cluster_ids))) %>%
      mutate(cluster_ids = factor(cluster_ids, levels = sort(unique(cluster_ids))))
    # plot
    ggplot(gg_tb, aes(x = cluster_ids, y = percent, fill = value)) +
      geom_bar(stat = "identity") +
      { if (length(unique(gg_tb$value)) > 11){
        scale_fill_manual(name = "cell_type", values = colorRampPalette(
          brewer.pal(n = 11, name = "Set3"))(length(unique(gg_tb$value))))
      } else {
        scale_fill_brewer(name = "cell_type", palette = "Set3")
      }} +
      ggtitle(make_plot_title("Percent assigned cell type by cluster"))
    ggsave(paste0(out_graph, "cell_type_by_cluster_percent_stackedbar.png"),
      width = 12, height = 6, limitsize = FALSE)

}

plot_cell_type_and_layer_by_cluster_stacked_bar_plot <- function(
  seurat_obj, cluster_col_name = "RNA_snn_res.0.6"){

  print("plot_cell_type_and_layer_by_cluster_stacked_bar_plot")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # collect metadata
  metadata_tb <- seurat_obj[[]] %>%
    as_tibble %>%
    select_if(negate(is.numeric)) %>%
    # make sure clusters plot 0 to highest number
    mutate(cluster_ids = factor(
      cluster_ids
      , levels = cluster_ids %>%
        as.character() %>%
        as.numeric() %>%
        unique() %>%
        sort
      )) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

    # percent stacked bar chart
    # calculate percent
    gg_tb <- metadata_tb %>%
      filter(variable == "cell_type_and_layer") %>%
      group_by(cluster_ids) %>%
      count(value) %>%
      mutate(percent = n/sum(n)) %>%
      mutate(variable = "cell_type_and_layer") %>%
      ungroup() %>%
      mutate(cluster_ids = as.numeric(as.character(cluster_ids))) %>%
      mutate(cluster_ids = factor(cluster_ids, levels = sort(unique(cluster_ids))))
    # plot
    ggplot(gg_tb, aes(x = cluster_ids, y = percent, fill = value)) +
      geom_bar(stat = "identity") +
      { if (length(unique(gg_tb$value)) > 11){
        scale_fill_manual(name = "cell_type", values = colorRampPalette(
          brewer.pal(n = 11, name = "Set3"))(length(unique(gg_tb$value))))
      } else {
        scale_fill_brewer(name = "cell_type", palette = "Set3")
      }} +
      ggtitle(make_plot_title("Percent assigned cell type by cluster"))
    ggsave(paste0(out_graph, "cell_type_and_layer_by_cluster_percent_stackedbar.png"),
      width = 12, height = 6, limitsize = FALSE)

}

plot_dim_reduction_colored_by_cell_type_assignment <- function(
  seurat_obj, reduction = "umap"){

  print("plot_dim_reduction_colored_by_cell_type_assignment")

  cell_type_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble() %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[[c("RNA_snn_res.0.6", "cell_type")]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble()
    )

  plot_dim_reduction_colored_by_variable(
    dim_1 = cell_type_tb$dim_1
    , dim_2 = cell_type_tb$dim_2
    , variable_value = cell_type_tb[["cell_type"]]
    , title = make_plot_title(paste0("Cell type assigned by max cell type marker score", "\n", reduction))
    , legend_title = "cell_type"
    , size = 0.5
    , alpha = 0.5
    , guide_size = 2
  )

  ggsave(paste0(out_graph, "cell_type_", reduction, ".png")
    , width = 7, height = 6, dpi = 200)

}

plot_dim_reduction_colored_by_cell_type_and_layer_assignment <- function(
  seurat_obj, reduction = "umap"){

  print("plot_dim_reduction_colored_by_cell_type_assignment")

  cell_type_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble() %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[[c("RNA_snn_res.0.6", "cell_type_and_layer")]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble()
    )

  plot_dim_reduction_colored_by_variable(
    dim_1 = cell_type_tb$dim_1
    , dim_2 = cell_type_tb$dim_2
    , variable_value = cell_type_tb[["cell_type_and_layer"]]
    , title = make_plot_title(paste0("Cell type assigned by max cell type marker score", "\n", reduction))
    , legend_title = "cell_type"
    , size = 0.5
    , alpha = 0.5
    , guide_size = 2
  )

  ggsave(paste0(out_graph, "cell_type_and_layer_", reduction, ".png")
    , width = 7, height = 6, dpi = 200)

}

plot_dim_reduction_colored_by_cluster_cell_type_assignment <- function(
  seurat_obj, reduction = "umap"){

  print("plot_dim_reduction_colored_by_cluster_cell_type_assignment")

  cell_type_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble() %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[[c("RNA_snn_res.0.6", "cluster_cell_type_and_layer")]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble()
    )

  plot_dim_reduction_colored_by_variable(
    dim_1 = cell_type_tb$dim_1
    , dim_2 = cell_type_tb$dim_2
    , variable_value = cell_type_tb[["cluster_cell_type_and_layer"]]
    , title = make_plot_title(paste0("Cell type assigned to cluster by cell type marker score", "\n", reduction))
    , legend_title = "cluster_cell_type_and_layer"
    , size = 0.5
    , alpha = 0.5
    , guide_size = 2
  )

  ggsave(paste0(out_graph, "cluster_cell_type_", reduction, ".png")
    , width = 7, height = 6, dpi = 200)

}
################################################################################

### dimensionality reduction plot colored by expression of genes of interest

plot_genes_of_interest_expression_dim_reduction <- function(
  seurat_obj,
  genes,
  cluster_col_name = "cluster_ids",
  reduction = "tsne",
  out_suffix = "genes_of_interest"){

  print("plot_genes_of_interest_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = genes
    , seurat_obj = seurat_obj
    , expression_slot = "data"
    , seurat_cluster_col = "cluster_ids"
    , expression_color_gradient = TRUE
    , reduction = reduction
    , ncol = 4
    , limit_high = 2
    , limit_low = 0
    , title = make_plot_title("Mean normalized expression of genes of interest")
  )
  ggsave(paste0(out_graph, "expression_", reduction, "_", out_suffix, ".png")
    , width = 18, height = 10, dpi = 200)

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = genes
    , seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , seurat_cluster_col = "cluster_ids"
    , legend_title = "expression\nz-score"
    , reduction = reduction
    , zscore = TRUE
    , ncol = 4
    , limit_high = 1
    , limit_low = -1
    , title = make_plot_title("Mean normalized expression of genes of interest")
  )
  ggsave(
    paste0(out_graph, "expression_zscore_", reduction, "_", out_suffix, ".png")
    , width = 18, height = 10, dpi = 200)

}
################################################################################

plot_ad_gwas_genes_expression_dim_reduction <- function(
  seurat_obj, cluster_col_name = "cluster_ids", reduction = "umap"){

  print("plot_ad_gwas_genes_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # clean
  polo_ad_gwas_genes_tb <- polo_ad_gwas_genes_tb %>%
    clean_variable_names()

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = polo_ad_gwas_genes_tb$gwas_gene
    , groupings = polo_ad_gwas_genes_tb$`disease/trait`
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , expression_color_gradient = TRUE
    , reduction = reduction
    , ncol = 4
    , limit_high = 2
    , limit_low = 0
    , title = make_plot_title("Mean normalized expression of AD GWAS genes")
  )
  ggsave(paste0(out_graph, "ad_gwas_genes_expression_", reduction, ".png")
    , width = 18, height = 28, dpi = 200)

}
################################################################################

### top cluster enriched genes expression heatmap

plot_cluster_enriched_genes_heatmap <- function(
  seurat_obj, cluster_enriched_de_tb){

  print("plot_cluster_enriched_genes_heatmap")

  # gene_group_tb <- cluster_enriched_de_tb %>%
  #   filter(avg_logFC > 0) %>%
  #   group_by(cluster) %>%
  #   top_n(n = 20) %>%
  #   select(gene, cluster)
  gene_group_tb <- cluster_enriched_de_tb %>%
    filter(log2_fold_change > 0) %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = log2_fold_change) %>%
    select(gene, cluster) %>%
    rename(group = cluster)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$group
    , cluster_col = "cluster_ids"
    , plot_title = make_plot_title("Top enriched genes for each cluster")
    , limit_high = 4
  )
  ggsave(paste0(out_graph, "cluster_top_enriched_heatmap.png")
    , height = 30)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$group
    , cluster_col = "cluster_ids"
    , plot_title = make_plot_title("Top enriched genes for each cluster")
    , limit_low = -1
    , limit_high = 1
    , zscore = TRUE
  )
  ggsave(paste0(out_graph, "cluster_top_enriched_heatmap_zscore.png")
    , height = 30)

}
################################################################################

### top cluster expressed genes expression heatmap

plot_cluster_top_expressed_genes_heatmap <- function(
  seurat_obj, top_expressed_genes_tb){

  print("plot_cluster_enriched_genes_heatmap")

  gene_group_tb <- top_expressed_genes_tb %>%
    group_by(cluster) %>%
    top_n(n = 40, wt = mean_expression)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$cluster
    , cluster_col = "cluster_ids"
    , plot_title = make_plot_title("Top expressed genes for each cluster")
    , limit_high = 4
  )
  ggsave(paste0(out_graph, "cluster_top_expressed_heatmap.png")
    , height = 6*length(unique(top_expressed_genes_tb$cluster))
    , limitsize = FALSE)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$cluster
    , cluster_col = "cluster_ids"
    , plot_title = make_plot_title("Top expressed genes for each cluster")
    , limit_low = -1
    , limit_high = 1
    , zscore = TRUE
  )
  ggsave(paste0(out_graph, "cluster_top_expressed_heatmap_zscore.png")
    , height = 6*length(unique(top_expressed_genes_tb$cluster))
    , limitsize = FALSE)

}
################################################################################

### cluster enriched genes violin plots

plot_cluster_enriched_genes_violin_plots <- function(
  seurat_obj, cluster_enriched_de_tb){

  print("plot_cluster_enriched_genes_violin_plots")

  gene_group_tb <- cluster_enriched_de_tb %>%
    filter(log2_fold_change > 0) %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = log2_fold_change) %>%
    select(gene, cluster) %>%
    rename(group = cluster)

  # cell id cluster id key
  cellid_clusterid_tb <- Idents(seurat_obj) %>%
    enframe(name = "cell_id", value = "cluster")

  # expression z-scores
  mean_expr_zscores_tb <-
    # get expression data and subset to genes of interest
    GetAssayData(seurat_obj, slot = "data")[
      rownames(GetAssayData(seurat_obj, slot = "data")) %in%
        gene_group_tb$gene, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    gather("cell_id", "expression", -gene) %>%
    # add gene group info
    right_join(., gene_group_tb, by = c("gene" = "gene")) %>%
    # add sub-clusters
    left_join(., cellid_clusterid_tb)

  mean_expr_zscores_tb %>%
    mutate(gene = factor(gene, levels = gene_group_tb$gene)) %>%
    # group_by(group) %>%
    ggplot(aes(x = cluster, y = expression, fill = cluster)) +
      facet_wrap(group~gene, scales = "free") +
      geom_violin() +
      ggtitle(make_plot_title(
        "Top 5 cluster enriched genes for each cluster")) +
      ggplot_set_theme_publication
  ggsave(paste0(out_graph, "cluster_enriched_violin.png")
    , height = 2+2*length(unique(cluster_enriched_de_tb$cluster))
    , width = 10)

}
################################################################################

### Percentage of cell types

make_plots_of_percent_of_cell_types <- function(seurat_obj){
  print("make_plots_of_percent_of_cell_types")

  seurat_obj[[]] %>%
    as_tibble() %>%
    select(cluster = cluster_ids, cell_type) %>%
    filter(cell_type != "NA") %>%
    group_by(cell_type) %>%
    summarise (n = n()) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    ggplot(aes(x = cell_type, y = percent)) +
      geom_bar(position = "dodge", stat = "identity") +
      coord_flip() +
      ggtitle(make_plot_title("Percent of cells by type"))
      ggsave(paste0(out_graph, "percent_cell_types_barplot.pdf"))

  seurat_obj[[]] %>%
    as_tibble() %>%
    filter(clinical_dx == "Control") %>%
    select(cluster = cluster_ids, cell_type) %>%
    filter(cell_type != "NA") %>%
    group_by(cell_type) %>%
    summarise (n = n()) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    ggplot(aes(x = cell_type, y = percent)) +
      geom_bar(position = "dodge", stat = "identity") +
      coord_flip() +
      ggtitle(make_plot_title("Percent of cells by type for control samples"))
      ggsave(paste0(out_graph, "percent_cell_types_control_barplot.pdf"))

  seurat_obj[[]] %>%
    as_tibble() %>%
    filter(clinical_dx == "Control") %>%
    select(cluster = cluster_ids, library_id, cell_type) %>%
    filter(cell_type != "NA") %>%
    group_by(library_id, cell_type) %>%
    summarise (n = n()) %>%
    group_by(library_id) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    ggplot(aes(x = cell_type, y = percent)) +
      geom_boxplot() +
      coord_flip() +
      ggtitle(make_plot_title("Percent of cells by type for each library for control samples"))
      ggsave(paste0(out_graph, "percent_cell_types_control_boxplot.pdf"))

  # number of cells in clinical_dx with least number of cells
  # this is for down-sampling
  min_number_cells <- seurat_obj[[c("clinical_dx")]] %>%
    table() %>%
    as_tibble() %>%
    rename(number_cells = "n", clinical_dx= ".") %>%
    filter(number_cells == min(number_cells)) %>%
    pull(number_cells)

  seurat_obj[[]] %>%
    as_tibble() %>%
    select(cluster = cluster_ids, clinical_dx, library_id, cell_type) %>%
    group_by(clinical_dx) %>%
    # down-sample to equal number of cells per clinical_dx
    sample_n(size = min_number_cells) %>%
    ungroup() %>%
    filter(cell_type != "NA") %>%
    group_by(library_id, cell_type, clinical_dx) %>%
    summarise(n = n()) %>%
    group_by(library_id) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    mutate(clinical_dx = factor(
      clinical_dx, levels = c("Control", "AD", "bvFTD", "PSP-S"))) %>%
    ggplot(aes(x = cell_type, y = percent, fill = clinical_dx)) +
      geom_boxplot() +
      coord_flip() +
      ggtitle(make_plot_title(
        paste0("Percent of cells by type for each library",
          "\nDown-sampled to equal number of cells per clinical dx")))
      ggsave(paste0(out_graph, "percent_cell_types_clinical_dx_boxplot.pdf"))

  dat <- seurat_obj[[]] %>%
    as_tibble() %>%
    select(cluster = cluster_ids, clinical_dx, library_id, region, cell_type) %>%
    group_by(clinical_dx) %>%
    # down-sample to equal number of cells per clinical_dx
    sample_n(size = min_number_cells) %>%
    ungroup() %>%
    filter(cell_type != "NA") %>%
    group_by(library_id, cell_type, clinical_dx, region) %>%
    summarise(n = n()) %>%
    group_by(library_id) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    mutate(clinical_dx = factor(
      clinical_dx, levels = c("Control", "AD", "bvFTD", "PSP-S")))

  my_comparisons <- list(c(1, 2), c(1, 3), c(1, 4))
  dat %>%
    ggplot(aes(x = clinical_dx, y = percent, fill = clinical_dx)) +
      geom_boxplot() +
      facet_grid(cell_type~region, scales = "free_y") +
      stat_compare_means(aes(group = clinical_dx)
        , comparisons = my_comparisons
        , method = "t.test"
        , p.adjust.method = "none"
        , label = "p.signif"
        # , label.y = max(dat$percent)
      ) +
      ggtitle(make_plot_title(
        paste0("Percent of cells by type for each library",
          "\nDown-sampled to equal number of cells per clinical dx")))
      ggsave(paste0(out_graph, "percent_cell_types_clinical_dx_region_boxplot.pdf"),
        width = 7, height = 12)

  dat <- seurat_obj[[]] %>%
    as_tibble() %>%
    select(cluster = cluster_ids, clinical_dx, library_id, region, cluster_cell_type_and_layer) %>%
    group_by(clinical_dx) %>%
    # down-sample to equal number of cells per clinical_dx
    sample_n(size = min_number_cells) %>%
    ungroup() %>%
    filter(cluster_cell_type_and_layer != "NA") %>%
    group_by(library_id, cluster_cell_type_and_layer, clinical_dx, region) %>%
    summarise(n = n()) %>%
    group_by(library_id) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    mutate(clinical_dx = factor(
      clinical_dx, levels = c("Control", "AD", "bvFTD", "PSP-S")))

  my_comparisons <- list(c(1, 2), c(1, 3), c(1, 4))
  dat %>%
    ggplot(aes(x = clinical_dx, y = percent, fill = clinical_dx)) +
      geom_boxplot() +
      facet_grid(cluster_cell_type_and_layer~region, scales = "free_y") +
      stat_compare_means(aes(group = clinical_dx)
        , comparisons = my_comparisons
        , method = "t.test"
        , p.adjust.method = "none"
        , label = "p.signif"
        , hide.ns = TRUE
        # , label.y = max(dat$percent)
      ) +
      ggtitle(make_plot_title(
        paste0("Percent of cells by type assigned by cluster annotation for each library",
          "\nDown-sampled to equal number of cells per clinical dx")))
      ggsave(paste0(out_graph, "percent_cell_types_and_layer_clinical_dx_region_boxplot.pdf"),
        width = 6, height = 22)

    dat <- seurat_obj[[]] %>%
      as_tibble() %>%
      select(cluster = cluster_ids, clinical_dx, library_id, region, cell_type_and_layer) %>%
      group_by(clinical_dx) %>%
      # down-sample to equal number of cells per clinical_dx
      sample_n(size = min_number_cells) %>%
      ungroup() %>%
      filter(cell_type_and_layer != "NA") %>%
      group_by(library_id, cell_type_and_layer, clinical_dx, region) %>%
      summarise(n = n()) %>%
      group_by(library_id) %>%
      mutate(percent = (n / sum(n)) * 100) %>%
      mutate(clinical_dx = factor(
        clinical_dx, levels = c("Control", "AD", "bvFTD", "PSP-S")))

    my_comparisons <- list(c(1, 2), c(1, 3), c(1, 4))
    dat %>%
      ggplot(aes(x = clinical_dx, y = percent, fill = clinical_dx)) +
        geom_boxplot() +
        facet_grid(cell_type_and_layer~region, scales = "free_y") +
        stat_compare_means(aes(group = clinical_dx)
          , comparisons = my_comparisons
          , method = "t.test"
          , p.adjust.method = "none"
          , label = "p.signif"
          , hide.ns = TRUE
          # , label.y = max(dat$percent)
        ) +
        ggtitle(make_plot_title(
          paste0("Percent of cells by type assigned by markers for each library",
            "\nDown-sampled to equal number of cells per clinical dx")))
        ggsave(paste0(out_graph, "percent_cell_types_and_layer_clinical_dx_region_boxplot_2.pdf"),
          width = 7, height = 16)

}
################################################################################

plot_number_of_cells_expressing_genes_box_jitter_plot <- function(seurat_obj, genes){

  print("plot_number_of_cells_expressing_genes_box_jitter_plot")

  FetchData(seurat_obj, vars = c(genes, "clinical_dx", "library_id")) %>%
   gather(key = "gene", value = "expression", -clinical_dx, -library_id) %>%
   mutate(expressed = (expression > 0)) %>%
   filter(expressed == TRUE) %>%
   count(clinical_dx, gene, library_id, expressed) %>%
   ggplot(aes(x = clinical_dx, y = n, fill = clinical_dx)) +
     geom_boxplot(outlier.shape = NA) +
     # geom_point(size = 0.5, position = position_jitterdodge(dodge.width = 0.2)) +
     facet_wrap(~gene, scale = "free_y") +
     geom_jitter(size = 0.5) +
     ylab("Number of cells expressing gene") +
     ggtitle(make_plot_title("Number of cells expressing gene by clinical_dx"))
     ggsave(paste0(out_graph, "number_of_cells_expressing_gene_box_jitter_plot.png"), width = 10)

  print("end of plot_number_of_cells_expressing_genes_box_jitter_plot")

}
################################################################################
