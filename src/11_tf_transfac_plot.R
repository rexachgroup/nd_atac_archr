liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "circlize", "scales", "RColorBrewer", "ComplexHeatmap")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

base_dir <- normalizePath("../data/archr/atac-2020-all")
proj_dir <- file.path(base_dir, "10_tf_transfac")
base_out_dir <- file.path(base_dir, "11_tf_transfac_plot")
cluster_args <- file.path(proj_dir, "args_tb.rds")

RESOURCES <- list(
    ncpus = 8,
    memory = 4 * 8,
    walltime = 86400
)

writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), str_glue(setBold, msg, setNorm), side = "both"))
}
alert <- function() { system("tput bel") }

main <- function() {
    dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)
    args_tb <- readRDS(cluster_args)
    args_tb <- args_tb %>% mutate(
        proj_dir = out_dir,
        out_dir = file.path(base_out_dir, project_names)
    )

    # load + reformat SummarisedExperiment into matrix
    args_tb <- mutate(args_tb, marker_signif_tb = pmap(args_tb, function(...) {
            cr <- list(...)
            marker_signif_list <- readRDS(file.path(cr$proj_dir, "motifMatrix_dx.rds"))
            marker_signif_tbs <- map(marker_signif_list, function(marker_signif) {
                signif_tb <- to_tb(marker_signif)
                signif_tb <- signif_tb %>%
                    mutate(DirPval = MeanDiff * -log(Pval))
            })
            bind_rows(marker_signif_tbs, .id = "dx")
        })
    )
    
    pmap(args_tb, function(...) {
        cr <- list(...)
        # reformat taipale names?
        signif_heatmap_tb <- cr$marker_signif_tb

        signif_heatmap_tb <- signif_heatmap_tb %>% 
            filter(Pval < 0.05) %>%
            mutate(name = fct_reorder(name, DirPval, .fun = function(x) mean(abs(x)), .desc = TRUE)) %>%
            arrange(name)

        signif_hclust <- signif_heatmap_tb %>%
            pivot_matrix("dx", "DirPval", "name") %>%
            replace(is.na(.), 0) %>%
            dist %>%
            hclust

        pdf(file.path(base_out_dir, paste0(cr$project_names, "_dirpval.pdf")),
            width = 10,
            height = 0.2 * length(unique(signif_heatmap_tb$name)))
        draw(heatmap_text_matrix(signif_heatmap_tb,
            "dx", "name", "DirPval", 
            label_name = "DirPval", row_hclust = signif_hclust,
            na_val_fill = 0
        ))
        graphics.off()
    })
}

to_tb <- function(marker_signif) {
    rd <- as.data.frame(rowData(marker_signif))
    marker_assays <- setNames(Reduce(cbind, assays(marker_signif)), nm = names(assays(marker_signif)))
    marker_tb <- as_tibble(bind_cols(rd, marker_assays))
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from) {
    tb_pivot <- tb %>%
        select(all_of(c(cols_from, values_from, rows_from))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = values_from) %>%
        ungroup

    tb_matrix <- select(tb_pivot, -rows_from) %>% as.matrix()

    rownames(tb_matrix) <- tb_pivot %>% 
        rowwise() %>%
        summarize(join_rowname = paste(c_across(rows_from), collapse = "|"), .groups = "drop") %>%
        pluck("join_rowname")

    glimpse(tb_matrix)
    return(tb_matrix)
}

heatmap_text_matrix <- function(tb, 
        row_str, col_str, val_str, txt_str = NULL, 
        row_split = NULL, col_split = NULL, 
        row_hclust = FALSE, col_hclust = FALSE,
        label_name = NULL, na_val_fill = NA, ...
){
    heatmap_val_matrix <- pivot_matrix(tb, row_str, val_str, col_str)
    if (!is.na(na_val_fill)) {
        heatmap_val_matrix <- replace(heatmap_val_matrix, is.na(heatmap_val_matrix), na_val_fill)
    }
    if (!is.null(txt_str)) {
        heatmap_text_matrix <- pivot_matrix(tb, row_str, txt_str, col_str)
    }
    row_annot <- tb %>%
        group_by(.data[[col_str]]) %>%
        slice_head %>% ungroup
    col_annot <- tb %>%
        group_by(.data[[row_str]]) %>%
        slice_head %>% ungroup
    
    if (is.logical(row_hclust) && !row_hclust) {
        plot_row_order <- rownames(heatmap_val_matrix)
    } else {
        plot_row_order <- rownames(heatmap_val_matrix)[row_hclust$order]
    }
    plot_row_label <- map_chr(plot_row_order, function(x) {
        if (str_starts(x, fixed("taipale"))) {
            ma <- str_match(x, "taipale-(\\w+)-(\\w+)-(\\w+)")
            return(str_glue("taipale-{ma[, 3]}-{ma[, 4]}"))
        } else { 
            return(x) 
        }
    })
    
    if (is.logical(col_hclust) && !col_hclust) {
        plot_col_order <- colnames(heatmap_val_matrix) #col_annot[[row_str]]
        plot_col_label <- as.character(plot_col_order)
    }
    
    plot_text_label <- function(j, i, x, y, w, h, col) {
        label <- heatmap_text_matrix[i, j]
        if (!is.na(label)) {
            grid.text(label, x, y)
        }
    }
    colormap <- colorRamp2(
        breaks = c(quantile(heatmap_val_matrix, 0.05, na.rm = TRUE), 0, quantile(heatmap_val_matrix, 0.95, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )

    Heatmap(
        heatmap_val_matrix,
        col = colormap,
        na_col = "grey75",
        cluster_rows = row_hclust,
        cluster_columns = col_hclust,
        row_order = if(is.logical(row_hclust) && !row_hclust) { plot_row_order } else { NULL },
        column_order = if(is.logical(col_hclust) && !col_hclust) { plot_col_order } else { NULL },
        row_label = plot_row_label,
        column_label = plot_col_label,
        cell_fun = if(!is.null(txt_str)) { plot_text_label } else { NULL },
        name = if(is.null(label_name)) { "" } else { label_name },
        row_title_rot = 0,
        column_title_rot = 90,
        row_names_max_width = unit(20, "cm"),
        row_dend_side = "right",
        row_dend_width = unit(7, "cm"),
        row_title_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 10),
        row_names_side = "left",
        ...
    )
}

if (!interactive()) {
    main()
}
