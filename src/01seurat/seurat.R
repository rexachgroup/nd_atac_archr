# Seurat run.
library(Seurat)
library(Matrix)
library(data.table)

mk_seurat <- function(path, project_name = "", downsample = FALSE, ds_n_genes = 6000, ds_n_cells = 3000) {
	expr <- Seurat::Read10X(data.dir = path)
	seurat <- Seurat::CreateSeuratObject(counts = expr, project = project_name)
	return(seurat)
}

import_seurat <- function(path, project_name) {
	mk_seurat(path, project_name)
}

export_seurat <- function(path, seurat) {
	saveRDS(path, seruat, compress = FALSE)
}

main <- function() {
	cellranger_base <- "../../cellranger"
	cellranger_runs <- c("20200114/PFC1_1_ns/", "20200114/PFC1_3_ns/")
	cellranger_paths <- file.path(cellranger_base, cellranger_runs, "outs/filtered_feature_bc_matrix")
	names(cellranger_paths) <- c("PFC1_1", "PFC1_3")
	dir.exists(cellranger_paths)

	out_base <- "../../data"
	raw_path <- file.path(out_base, "raw_seurat")

	raw_seurat <- sapply(names(cellranger_paths), function(name) {
		path <- cellranger_paths[[name]]
		import_seurat(path = path, project_name = name)
	}, simplify = FALSE)

	lapply(names(raw_seurat), function(name) {
		print(file.path(raw_path, name))
	})
}

if (interactive()) main()
