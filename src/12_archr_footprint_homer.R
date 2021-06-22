# Call peaks and test differentially per subcluster.
liblist <- c("tidyverse", "batchtools", "readxl", "Seurat", "ArchR", "patchwork", "ggrepel")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))

base_dir <- normalizePath("../data/archr/atac-2020-all")
base_out_dir <- file.path(base_dir, "12_archr_footprint_homer")
batchtools <- file.path(base_out_dir, "footprint_batchtools")

cluster_args <- file.path(base_dir, "10_tf_homer", "args_tb.rds")
writeMsg <- function(msg) {
    writeLines(str_pad(width = getOption("width"), msg, side = "both"))
}

RESOURCES <- list(
    ncpus = 8,
    memory = 64,
    walltime = 86400,
    partition = "bigmem"
)

motif_names <- c("NR3C1", "RUNX1", "RUNX2", "PU.1", "SPI1")

main <- function() {
    dir.create(base_out_dir)
    args_tb <- readRDS(cluster_args)

    args_tb <- args_tb %>%
        mutate(proj_dir = out_dir) %>%
        mutate(out_dir = file.path(base_out_dir, project_names))

    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(batchtools)
    } else {
        reg <- loadRegistry(batchtools, writeable = TRUE)
    }

    # footprint_worker ; getPositions / getFootprints
    clearRegistry()
    reg$packages <- liblist
    batchExport(mget(ls()))
    ids <- batchMap(footprint_worker, args = list(proj_dir = args_tb$proj_dir, out_dir = args_tb$out_dir), more.args = list(motif_names = motif_names))
    setJobNames(ids, args_tb$project_names)
    submitJobs(ids, RESOURCES)
    waitForJobs()
    system("tput bel")

    saveRDS(args_tb, file.path(base_out_dir, "args_tb.rds"))

    #     map(args_tb$out_dir, function(out_dir) {
    #         project <- loadArchRProject(out_dir)
    #         footprint_dx_list <- readRDS(file.path(out_dir, "Footprint", "footprint_dx.rds"))
    #         plotFootprints(
    #             seFoot = footprint_dx_list,
    #             ArchRProj = project,
    #             normMethod = "subtract",
    #             plotName = "archr_footprint",
    #             addDOC = TRUE,
    #             smoothWindow = 5,
    #             force = TRUE
    #         )
    #     })

    # plotFootprints
    args_tb <- mutate(args_tb, 
        footprint_dx_data = map(out_dir, function(out_dir) {
            footprint_path <- file.path(out_dir, "Footprint", "footprint_dx.rds")
            writeMsg(str_glue("load {footprint_path}")) 
            seFoot <- readRDS(footprint_path)
            get_footprint_data(seFoot, smooth_window = 15, flank_norm = 100)
        }),
        footprint_pairwise_data = map(out_dir, function(out_dir) {
            footprint_pairwise_path <- file.path(out_dir, "Footprint", "footprint_dx_pairwise.rds")
            writeMsg(str_glue("load {footprint_pairwise_path}")) 
            seFoot_list <- readRDS(footprint_pairwise_path)
            map(seFoot_list, get_footprint_data, smooth_window = 15, flank_norm = 100)
        })
    )

    pwalk(args_tb, function(...) {
        cr <- list(...)
        foot_tb <- map(cr$footprint_dx_data, ~.x$foot_df) %>% bind_rows(.id = "tf")
        bias_tb <- map(cr$footprint_dx_data, ~.x$bias_df) %>% bind_rows(.id = "tf")
        write_csv(foot_tb, file.path(base_out_dir, str_glue(cr$project_names, "-footprint.csv")))
        write_csv(bias_tb, file.path(base_out_dir, str_glue(cr$project_names, "-footprint-biased.csv")))
    })

    pwalk(args_tb, function(...) {
        cr <- list(...)
        out_path <- file.path(base_out_dir, str_glue(cr$project_names, "-footprint_dx.pdf"))
        writeMsg(out_path)
        plotlist <- plot_footprint_data(cr$footprint_dx_data)
        pdf(out_path, width = 4, height = 6)
        walk(plotlist, print)
        dev.off()
    })

    pwalk(args_tb, function(...) {
        cr <- list(...)
        iwalk(cr$footprint_pairwise_data, function(x, n) {
            out_path <- file.path(base_out_dir, str_glue(cr$project_names, "-footprint_{n}.pdf"))
            writeMsg(out_path)
            plotlist <- plot_footprint_data(x)
            pdf(out_path, width = 4, height = 6)
            walk(plotlist, print)
            dev.off()
        })
    })
}

footprint_worker <- function(proj_dir, out_dir, motif_names) {
    l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))
    startTime <- Sys.time()
    addArchRGenome("hg38")
    addArchRThreads(RESOURCES$ncpus)
    project <- loadArchRProject(proj_dir)  
    project <- saveArchRProject(project, out_dir, load = TRUE)
    annot_dir <- file.path(out_dir, "Footprint")
    plot_dir <- file.path(out_dir, "Plots")
    unlink(plot_dir)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)
    orig_dir <- getwd()
    setwd(out_dir)

    motif_pos <- getPositions(project)
    motif_subset <- unlist(lapply(motif_names, function(x) {
        grep(x, names(motif_pos), value = TRUE, fixed = TRUE)
    }))
    motif_granges <- motif_pos[motif_subset]

    writeMsg("footprint dx pairwise")
    footprint_dx_pairwise_list <- footprint_dx_pairwise(project, motif_granges)
    saveRDS(footprint_dx_pairwise_list, file.path(annot_dir, "footprint_dx_pairwise.rds"))
    writeMsg("footprint dx")
    footprint_dx_list <- getFootprints(project, positions = motif_granges, groupBy = "Clinical.Dx")
    saveRDS(footprint_dx_list, file.path(annot_dir, "footprint_dx.rds"))

}

footprint_dx_pairwise <- function(project, motif_granges) {
    dx <- c("AD", "bvFTD", "PSP-S") 
    dx_wk <- function(dx) {
        writeMsg(str_glue("footprint_dx_pairwise {dx}"))
        getFootprints(
            project,
            positions = motif_granges,
            groupBy = "Clinical.Dx",
            useGroups = c("Control", dx)
        )
    }
    setNames(map(dx, ~tryCatch(dx_wk(.), error = print)), dx)
}

# Get all data up to plotting -- basically plotFootprint / .ggFootprint without ggplot calls
get_footprint_data <- function(seFoot, smooth_window = 5, flank = 250, flank_norm = 50) {
    foot_names <- names(assays(seFoot))
    row_data <- rowData(seFoot)
    col_data <- colData(seFoot)
    
    foot_rows <- which(row_data$type == "footprint")
    bias_rows <- which(row_data$type == "bias")

    footprint_data_list <- map(foot_names, function(foot_name) {
        writeMsg(str_glue("data for {foot_name}"))
        foot_mat <- assay(seFoot, foot_name)[foot_rows, ]
        bias_mat <- assay(seFoot, foot_name)[bias_rows, ]
        foot_row_df <- row_data[foot_rows, ]
        bias_row_df <- row_data[bias_rows, ]

        writeMsg(str_glue("applying smoothing window={smooth_window} to footprint"))
        foot_mat <- apply(foot_mat, 2, function(x) .centerRollMean(x, smooth_window))
        bias_mat <- apply(bias_mat, 2, function(x) .centerRollMean(x, smooth_window))

        writeMsg("normalizing by flanking regions")
        idx <- which(abs(foot_row_df$x) >= flank - flank_norm)
        foot_mat <- t(t(foot_mat) / colMeans(foot_mat[idx, ,drop=FALSE]))
        bias_mat <- t(t(bias_mat) / colMeans(bias_mat[idx, ,drop=FALSE]))

        writeMsg("Tn5 bias division")
        foot_mat <- foot_mat / bias_mat

        writeMsg("mean / sd for footprint")
        foot_mat_mean <- .groupMeans(foot_mat, col_data$Group)
        foot_mat_sd <- .groupSds(foot_mat, col_data$Group)
        bias_mat_mean <- .groupMeans(bias_mat, col_data$Group)
        bias_mat_sd <- .groupSds(bias_mat, col_data$Group)
        smooth_foot <- rowMaxs(apply(foot_mat_mean, 2, function(x) .centerRollMean(x, 11)))

        plot_ids <- seq_len(nrow(foot_mat_mean))
        plot_foot_df <- lapply(seq_len(ncol(foot_mat_mean)), function(x) {
            data.frame(
                x = foot_row_df$x,
                mean = foot_mat_mean[, x],
                sd = foot_mat_sd[, x],
                group = colnames(foot_mat_mean)[x]
            )[plot_ids, ,drop = FALSE]
        }) %>% bind_rows %>% as_tibble

        plot_bias_df <- lapply(seq_len(ncol(bias_mat_mean)), function(x) {
            data.frame(
                x = bias_row_df$x,
                mean = bias_mat_mean[, x],
                sd = bias_mat_sd[, x],
                group = colnames(bias_mat_mean)[x]
            )[plot_ids, ,drop = FALSE]
        }) %>% bind_rows %>% as_tibble
        
        return(list(foot_df = plot_foot_df, bias_df = plot_bias_df, foot_smooth_df = smooth_foot))
    }) %>% setNames(., foot_names)
}

# HiddenUtils.R::.centerRollMean
.centerRollMean <- function(v = NULL, k = NULL) {
    o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
    if(k%%2==0){
        o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
    }else if(k%%2==1){
        o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
    }else{
        stop("Error!")
    }
    o2
}

# HiddenUtils.R::.groupMeans
.groupMeans <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups)==ncol(mat))
    gm <- lapply(unique(groups), function(x){
        if(sparse){
             Matrix::rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
        }else{
             rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
        }
    }) %>% Reduce("cbind",.)
    colnames(gm) <- unique(groups)
    return(gm)
}

# HiddenUtils.R::.groupSds
.groupSds <- function(mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups)==ncol(mat))
    gs <- lapply(unique(groups), function(x){
        if (sparse){
            matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
        }else{
            matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind",.)
    colnames(gs) <- unique(groups)
    return(gs)
}

plot_footprint_data <- function(footprint_list) {
    imap(footprint_list, function(data, n) {
        tryCatch({
            writeMsg(str_glue("plotting {n}"))

            pal <- paletteDiscrete(values = unique(data$foot_df$group))

            plot_max <- data$foot_df[order(data$foot_df$mean, decreasing = TRUE), ]
            plot_max <- plot_max[abs(plot_max$x) > 20 & abs(plot_max$x) < 50, ]
            plot_max <- plot_max[!duplicated(plot_max$group), ]

            plot_max <- data$foot_df %>%
                arrange(desc(mean)) %>%
                filter(abs(x) > 20 & abs(x) < 50) %>%
                group_by(group) %>% slice_head %>% ungroup %>%
                mutate(x = 25)

            gg_foot <- ggplot(data$foot_df, aes(x = x, y = mean, color = group)) +
                geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
                geom_line() +
                geom_label_repel(data = plot_max, aes(label = group), size = 3, xlim = c(75, NA), show.legend = FALSE) +
                scale_color_manual(values = pal) +
                scale_fill_manual(values = pal) +
                coord_cartesian(
                    expand = FALSE,
                    xlim = c(min(data$foot_df$x), max(data$foot_df$x)),
                    ylim = c(quantile(data$foot_df$mean, 0.0001), quantile(data$foot_smooth_df, 0.999))
                ) +
                labs(
                    x = "Distance to motif center (bp)",
                    y = str_glue("Tn5 Bias Division Normalized Insertions"),
                    title = n
                ) +
                theme_ArchR(baseSize = 6) +
                theme(legend.position = "bottom", legend.box.background = element_rect(color = NA))

            gg_bias <- ggplot(data$bias_df, aes(x = x, y = mean, color = group)) +
                geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
                geom_line() +
                scale_color_manual(values = pal) +
                scale_fill_manual(values = pal) +
                coord_cartesian(
                    expand = FALSE,
                    xlim = c(min(data$bias_df$x), max(data$bias_df$x)),
                    ylim = c(quantile(data$bias_df$mean, 0.0001), quantile(data$bias_df$mean, 0.999))
                ) +
                labs(
                    x = "Distance to motif center (bp)",
                    y = str_glue("Tn5 Bias Normalized Insertions"),
                    title = n
                ) +
                theme_ArchR(baseSize = 6) +
                theme(legend.position = "none", legend.box.background = element_rect(color = NA))

            return(wrap_plots(gg_foot, gg_bias, nrow = 2, ncol = 1, heights = c(2, 1)))
        },.error = function(e) {
            writeMsg("==== error in {n} ====")
            writeMsg(as.character(e))
        })
    })
    
}

if (!interactive()) {
    main()
}
