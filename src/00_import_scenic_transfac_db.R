# Import and merge scenic results with transfact motif db.
liblist <- c("tidyverse", "TFBSTools")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))

scenic_motif_db <- normalizePath("../data/scenic.pfm.rda", mustWork = T)
motif_regulon_filter <- normalizePath("../data/excitatory_preCGregulons_genes_GRNBoost_weights.csv", mustWork = T)
transfac_motif_db <- normalizePath("../data/transfac_pwm.rds", mustWork = T)
base_dir <- normalizePath("../data/archr/atac-2020-all")
OUT_DIR <- file.path(base_dir, "00_import_scenic_transfac_db/")

main <- function() {
    dir.create(OUT_DIR)
    # load scenic
    load(scenic_motif_db)
    motif_tb <- read_csv(motif_regulon_filter)
    pwm <- TFBSTools::toPWM(pfm)
    scenic_tf <- grep_tf(pwm, motif_tb$TF)

    # load transfac motif db
    
    transfac_pwm <- readRDS(transfac_motif_db)
    transfac_tf <- grep_tf(transfac_pwm, motif_tb$TF)
    combined_motiflist <- c(scenic_tf, transfac_tf)
    stopifnot(names(combined_motiflist) == name(combined_motiflist))
    write.csv(name(combined_motiflist), file.path(OUT_DIR, "scenic_transfac_motifs.csv"))
    saveRDS(combined_motiflist, file.path(OUT_DIR, "scenic_transfac_db.rds"))
}

grep_tf <- function(pwmatrixlist, motif_names) {
    names(pwmatrixlist) <- name(pwmatrixlist)
    pwm_indices <- unlist(map(unique(motif_names), function(tf) {
        grep(tf, name(pwmatrixlist), fixed = TRUE)
    }))
    pwm_indices <- pwm_indices[!duplicated(pwm_indices)]
    return(pwmatrixlist[pwm_indices])
}

if (!interactive()) main()
