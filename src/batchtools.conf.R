library(rprojroot)
root <- rprojroot::has_file(".git/index")
ORION_SLURM_TMPL = root$find_file("src/slurm.tmpl")
#HOFFMAN_SGE_TMPL = "sge.tmpl"

is_orion <- function() dir.exists("/geschwindlabshares")

if (is_orion()) {
    compress = FALSE
    cluster.functions <- makeClusterFunctionsSlurm(template = ORION_SLURM_TMPL)
}
