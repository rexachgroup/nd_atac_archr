# Seurat run.
library(Seurat)
library(Matrix)
library(data.table)
library(tidyverse)

mk_seurat <- function(path, project_name = "", downsample = FALSE, ds_n_genes = 6000, ds_n_cells = 3000) {
	expr <- Seurat::Read10X(data.dir = path)
	if (downsample) {
		expr <- expr[, sample(1:ncol(expr), ds_n_cells)]
	}
	seurat <- Seurat::CreateSeuratObject(counts = expr, project = project_name)
	return(seurat)
}

export_seurat <- function(path, seurat) {
	saveRDS(seurat, path, compress = FALSE)
}

add_metadata_seurat <- function(seurat, metadata) {
	
}

qc_seurat <- function(seurat, min_number_genes = 200, max_percent_mito = 5, min_cells_per_gene = 0) {
	cell_meta <- as.data.table(seurat@meta.data)
	gene_meta <- as.data.table(tibble::enframe(rowSums(GetAssayData(seurat, slot = "counts")), name = "gene", value = "number_cells"))
	gene_meta <- gene_meta[number_cells > min_cells_per_gene,]
	return(seurat)
}

process_seurat <- function(seurat, verbose = FALSE) {
	if(verbose) print("process_seurat")

	if(verbose) print("NormalizeData")
	seurat <- NormalizeData(object = seurat, verbose = verbose)

	if(verbose) print("FindVariableGenes")
	seurat <- FindVariableFeatures(seurat, selection.method = "vst", verbose = verbose)

	if(verbose) print("ScaleData")
	seurat <- ScaleData(seurat, do.scale = TRUE, do.center = TRUE, features = rownames(seurat), 
							verbose = verbose, num.cores = 8, do.par = TRUE)

	if(verbose) print("RunPCA")
	seurat <- RunPCA(seurat, online.pca = TRUE, npcs = 100, do.print = TRUE, pcs.print = 1:5, genes.print = 5, verbose = verbose)

	if(verbose) print("RunTSNE")
	seurat <- RunTSNE(seurat, dims.use = 1:100, do.fast = TRUE, nthreads = 8, 
						  tsne.method = "Rtsne", reduction = "pca", max_iter = 2000)

	if(verbose) print("RunUMAP")
	seurat <- RunUMAP(seurat, dims = 1:100)

	if(verbose) print("FindNeighbors")
	seurat <- FindNeighbors(object = seurat, reduction = "pca", dims = 1:100, nn.eps = 0, k.param = 30)

	if(verbose) print("FindClusters")
	seurat <- FindClusters(object = seurat, resolution = c(0.4,0.5,0.6,0.7,0.8), n.start = 100)

	return(seurat)
}

main <- function() {
	cellranger_base <- "../../cellranger"
	cellranger_runs <- c("aggr/PFC1_1")
	cellranger_paths <- file.path(cellranger_base, cellranger_runs, "outs/filtered_feature_bc_matrix")
	names(cellranger_paths) <- c("PFC1")
	dir.exists(cellranger_paths)

	out_base <- "../../data"
	raw_path <- file.path(out_base, "raw_seurat")

	raw_seurat <- sapply(names(cellranger_paths), function(name) {
		path <- cellranger_paths[[name]]
		mk_seurat(path = path, project_name = name, downsample = TRUE)
	}, simplify = FALSE)

	lapply(names(raw_seurat), function(name) {
		export_seurat(file.path(raw_path, paste0(name, ".rds")), raw_seurat[[name]])
	})

	#seurat_objs <- add_metadata
	seurat_objs <- lapply(raw_seurat, qc_seurat)
	seurat_objs <- lapply(seurat_objs, process_seurat, verbose = TRUE)
	
}

if (interactive()) main()
