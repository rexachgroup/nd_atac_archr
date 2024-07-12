library(ArchR)
library(parallel)
addArchRThreads(threads = 1)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")
library(BSgenome.Hsapiens.UCSC.hg38)
pathToMacs2 <- "Your_path_to_macs2"

dx <- "AD" # Control PSP_S bvFTD
projss <- loadArchRProject( paste("projSS_by_DX",dx, sep = "/") )
selected_clusters <- c("C18","C19","C21", "C22", "C23", "C24", "C25")
projss <- projss[projss@cellColData$Clusters %in% selected_clusters,]

#Making Pseudo-bulk Replicates
projss <- addGroupCoverages(ArchRProj = projss, groupBy = "Clusters", force=T)

#MACS2 callpeaks
projss <- addReproduciblePeakSet(ArchRProj = projss, groupBy = "Clusters", pathToMacs2 = pathToMacs2)
projss <- addPeakMatrix(projss)

markersPeaks <- getMarkerFeatures(
  ArchRProj = projss, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

peakSet <- projss@peakSet
# extract markerPeaks
for (group in selected_clusters){
  gr.peak.group <- peakSet[names(peakSet) %in% group]
  
  MarkerPeakFile <- paste("markerPeak", dx, file_suffix,sep=".")
  markerPeakSet <- readRDS(file.path(path_peak, MarkerPeakFile))
  MarkerPeak <- getMarkers(markerPeakSet, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = F)
  markerpeak <- MarkerPeak[[group]]
  gr.markerpeak <- GRanges(seqnames = markerpeak$seqnames,
                           ranges = IRanges(markerpeak$start, end = markerpeak$end))
  dx.peaks <-c(gr.peak.group, gr.markerpeak) %>% unique(.)
  
  if(!isEmpty(dx.peaks)){
    fname <- paste(dx, group, "bed", sep = ".")
    df <- data.frame(chr = seqnames(dx.peaks), 
                     start = dx.peaks@ranges@start, 
                     end = dx.peaks@ranges@start + dx.peaks@ranges@width - 1,
                     id = paste(group, 1:length(dx.peaks), sep = "."))  
    write.table(df, file = fname, quote = F, col.names = F, row.names = F, sep = "\t")
  }
}
