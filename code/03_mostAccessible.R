library(ggplot2)
library(BuenColors)
library(data.table)
library(BuenColors)
library(GenomicRanges)
#library(mgatk)
library(chromVARmotifs)
library(motifmatchr)

counts <- data.matrix(fread("zcat < ../data/17aug_commonATAC_normalized.txt.gz"))

counts <- sweep(counts, 2, colSums(counts), "/")
maxx <- colnames(counts)[max.col(counts,ties.method="first")]
#gi <- giniRows(counts)

makePopPlot <- function(name, bin, color){
  plotdf.b <- data.frame(tSNE1 = tsne_peaks$Y[,1], tSNE2 = tsne_peaks$Y[,2], X = bin)
  p1 <- ggplot(shuf(plotdf.b), aes(tSNE1, tSNE2, color = X)) + geom_point(size = 0.1, alpha = 0.1) +
    scale_color_manual(values = c("black", color)) + pretty_plot() +
    ggtitle(name) + theme(legend.position = "none") +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  cowplot::ggsave(p1, file = paste0("../plots/",name,".png"), width = 10, height = 8, units = "in", dpi = 300)
}

makePopPlot("BCell", grepl("^B.", maxx), "goldenrod2")
makePopPlot("T", grepl("T.", maxx), "coral2")
makePopPlot("Progenitor", as.factor(as.numeric(grepl("LTHSC", maxx)) + as.numeric(grepl("STHSC", maxx)) + as.numeric(grepl("MMP", maxx))), "dodgerblue")
#makePopPlot("Xchromosome", seqnames(peaks) == "chrX")
makePopPlot("Mo", grepl("^Mo.", maxx), "darkorchid2")
makePopPlot("MF_microglia", grepl("^MF.microglia.", maxx), "deeppink")
#makePopPlot("NKT", grepl("^NKT", maxx))
#makePopPlot("T4", grepl("^T.4", maxx))
makePopPlot("DC", grepl("^DC", maxx), "deeppink")
#makePopPlot("ILC", grepl("^ILC", maxx))
makePopPlot("Ep", grepl("^Ep", maxx), "blue")
#makePopPlot("SLN", grepl("SLN", maxx))
#makePopPlot("Treg", grepl("Treg", maxx))


Gini = function (x, corr = FALSE, na.rm = TRUE) 
{
  if (!na.rm && any(is.na(x))) 
    return(NA_real_)
  x <- as.numeric(na.omit(x))
  n <- length(x)
  x <- sort(x)
  G <- sum(x * 1L:n)
  G <- 2 * G/sum(x) - (n + 1L)
  if (corr) 
    G/(n - 1L)
  else G/n
}

vec <- sapply((1:dim(counts)[1]), function(i){
  Gini(counts[i,])
})

# Make Gini index plot and CTCF plot
data("mouse_pwms_v2")
library(BSgenome.Mmusculus.UCSC.mm10)
peaksd <- fread("zcat < /Volumes/dat/Research/BuenrostroResearch/IMMGEN/immgen_final/data/ImmGenATAC1219.peak.bed.gz")
peaks <- makeGRangesFromDataFrame(peaksd, seqnames.field = "V1", start.field = "V2", end.field = "V3")
library(BSgenome.Mmusculus.UCSC.mm10)
ctcf <- matchMotifs(mouse_pwms_v2["ENSMUSG00000005698_LINE295_Ctcf_D"], peaks, BSgenome.Mmusculus.UCSC.mm10, p.cutoff = 1e-05)


pal <- "brewer_spectra"
plotdf.b <- data.frame(tSNE1 = tsne_peaks$Y[,1], tSNE2 = tsne_peaks$Y[,2], Gini = (vec),
                       CTCF = 1:dim(tsne_peaks$Y)[1] %in% (which(SummarizedExperiment::assays(ctcf)[["motifMatches"]][,1])))

p1 <- ggplot(plotdf.b[order(abs(plotdf.b$Gini - 0.5)),], aes(tSNE1, tSNE2, color = Gini)) +
  scale_color_gradientn(colors = jdb_palette(pal)) + geom_point(size = 0.1, alpha = 0.3) +
  pretty_plot() +
  theme(legend.position = "none") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
cowplot::ggsave(p1, file = paste0("../plots/","GINI",".png"), width = 10, height = 8, units = "in", dpi = 300)

p1 <- ggplot(plotdf.b %>% arrange(CTCF), aes(tSNE1, tSNE2, color = CTCF)) + geom_point(size = 0.1, alpha = 0.1) +
  scale_color_manual(values = c("black", "dodgerblue")) + pretty_plot() +
  theme(legend.position = "none") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + ggtitle("CTCF")
cowplot::ggsave(p1, file = paste0("../plots/","CTCF",".png"), width = 10, height = 8, units = "in", dpi = 300)



