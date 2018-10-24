library(BuenColors)
library(data.table)
library(dplyr)
load("../data/tSNE_peaks.rda")
promoter <- fread("../data/peakAnnoVector.txt")
mdf <- data.frame(tSNE1 = tsne_peaks$Y[,1], tSNE2 = tsne_peaks$Y[,2], TSS = promoter$peak.position == "TSS")
p1 <- ggplot(mdf %>% arrange(TSS), aes(x = tSNE1, y = tSNE2, color = TSS)) + geom_point(size = 0.1, alpha = .10) +
  scale_color_manual(values = c("black", "dodgerblue")) + pretty_plot() +
  theme(legend.position = "none") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
cowplot::ggsave(p1, file = paste0("../plots/", "TSS", ".png"), width = 10, height = 8, units = "in", dpi = 300)
