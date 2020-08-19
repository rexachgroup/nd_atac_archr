CELLRANGER_BIN <- normalizePath(root$find_file("ext/cellranger-3.1.0/cellranger"))
CELLRANGER_REFERENCE <- normalizePath(root$find_file("ext/refdata-cellranger-GrCh38-3.0.0_premrna/"))
FASTQ_DIR <- normalizePath(root$find_file("data/00-cellranger-input"))

is_hoffman_cluster <- function() {
    grepl("hoffman", system("hostname -a", intern = TRUE))
}

is_nessie <- function() {
    grepl("nessie", system("hostname", intern = TRUE))
}

init_registry <- function(reg_path, clear = TRUE) {
    if (dir.exists(reg_path)) {
        reg <- loadRegistry(reg_path, writeable = TRUE)
    } else {
        reg <- makeRegistry(reg_path)
    }
    if (clear) clearRegistry()
    if (is_hoffman_cluster()) {
        reg$cluster.functions <- makeClusterFunctionsSGE(template =
            as.character(root$find_file("src/rscript_hoffman.tmpl")),
            fs.latency = 1)
        writeLines("Set cluster functions to SGE")
    } else {
        reg$cluster.functions <- makeClusterFunctionsMulticore()
        writeLines("Set cluster functions to multicore")
    }
    return(reg) 
}

bash_exec <- function(cmd) {
    system(paste("bash -c", shQuote(cmd)))
}
