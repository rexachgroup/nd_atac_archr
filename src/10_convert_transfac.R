library(universalmotif)
library(TFBSTools)

stopifnot("transfac2meme not found" = system("which transfac2meme > /dev/null") == 0)

system("transfac2meme ../data/matrix.dat > ../data/transfac_meme.dat")
transfac_meme <- read_meme("../data/transfac_meme.dat")
transfac_motifs <- convert_motifs(transfac_meme, "TFBSTools-PWMatrix")
transfac_pwm <- do.call(PWMatrixList, transfac_motifs)
saveRDS(transfac_pwm, "../data/transfac_pwm.rds")
