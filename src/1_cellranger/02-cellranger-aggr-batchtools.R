lapply(c("rprojroot",
         "batchtools",
         "tidyverse"),
       function(x) suppressPackageStartupMessages(require(x, character.only = TRUE)))

root <- has_file(".rprojroot")

source(root$find_file("src/1_cellranger/00_CONSTANTS.R"))

IN_ROOT <- root$find_file("data/01-cellranger-individual")
OUT_ROOT <- root$find_file("data/02-cellranger-aggregate")
REGISTRY <- root$find_file("batchtools/02-cellranger-aggregate")
EXCLUDE_SAMPLE_IDS <- c("PFC3_6")
ID = "PFC_a"

MEM_LIMIT_TOTAL <- 128
NCORES <- 16
MEM_LIMIT_HDATA_GB <- paste0(ceiling(MEM_LIMIT_TOTAL / NCORES), "G")

make_aggregation_table <- function(path) {
    sample_dirs <- list.dirs(path, recursive = FALSE)
    library_id <- basename(sample_dirs)
    molecule_h5 <- file.path(sample_dirs, "outs", "molecule_info.h5")
    tibble(library_id = library_id, molecule_h5 = molecule_h5)
}

make_metrics_aggregation <- function(path) {
    sample_dirs <- list.dirs(path, recursive = FALSE)
    library_id <- basename(sample_dirs)
    metrics <- file.path(sample_dirs, "outs", "metrics_summary.csv")
    bind_rows(map2(library_id, metrics, function(key, value) {
        csv <- read_csv(value)
        csv$library_id <- key
        return(csv %>% select(library_id, everything()))
    }))
}

aggregation_path <- file.path(OUT_ROOT, ID, "aggregation.csv")

aggregation <- make_aggregation_table(IN_ROOT) %>% 
    filter(!(library_id %in% EXCLUDE_SAMPLE_IDS))

write_csv(aggregation, aggregation_path)

metrics <- make_metrics_aggregation(IN_ROOT) %>%
    filter(!(library_id %in% EXCLUDE_SAMPLE_IDS))

write_csv(metrics, file.path(OUT_ROOT, ID, "metrics.csv"))

aggr_cmd <- stringr::str_glue(
    "cd {OUT_ROOT} &&",
    "{CELLRANGER_BIN} aggr",
    "--id={ID}",
    "--csv={aggregation_path}",
    "--normalize=none",
    "--localmem={MEM_LIMIT_TOTAL}",
    "--localcores={NCORES}",
    "--nosecondary",
    .sep = " "
)

reg <- init_registry(REGISTRY)
batchMap(bash_exec, list(aggr_cmd), reg = reg)
submitJobs(reg = reg, resources = list(
    h_data_str = MEM_LIMIT_HDATA_GB,
    h_rt = "24:00:00",
    pe_shared = NCORES,
    highp = TRUE,
    measure.memory = TRUE
))
waitForJobs()
