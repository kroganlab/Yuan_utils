evidence <- read.delim2(file = "/Users/yzhou/Downloads/PCA/test_data.txt", header = T, stringsAsFactors = F)

# PCA
library(dplyr)
library(package = "tidyr")
library(data.table)
rownames(evidence) <- evidence$PROTID
evidence <- evidence[, -1]
colnames(evidence) <- gsub("_SPC", "", colnames(evidence))
evidence.mat <- as.matrix(evidence)
pca.proteinRun <- prcomp(t(evidence.mat))
pcaDT <- as.data.table(pca.proteinRun$x, keep.rownames=TRUE)
pcaPercentDev <- round(100 * (pca.proteinRun$sdev)^2/sum((pca.proteinRun$sdev)^2), 1)

library(ggplot2)
library(ggrepel)
ggplot (pcaDT, aes(x=PC1, y=PC2, col = "red")) + 
  geom_point(alpha=0.8, size=7) + theme_bw() + 
  ggtitle ("PCA") + 
  geom_text(aes(label=rn), col = "black", size = 4)+ 
  xlab (sprintf ("PC1, %.1f%%", pcaPercentDev[1])) + 
  ylab (sprintf ("PC2, %.1f%%", pcaPercentDev[2]))
ggsave("PCA.pdf", width = 8, height = 6)


ggplot(pcaDT, aes(PC1, PC2)) + #, color=condition, shape=condition
  geom_point(alpha = .8, size = 5, na.rm = TRUE) +
  geom_label_repel(
    aes(label = rn),
    #label.size = 2,
    box.padding   = 0.8,
    point.padding = 0.1,
    segment.color = 'grey50', 
    max.overlaps = Inf
  ) +
  labs(title = "Principal Component Analysis",
       fill = "") +
  theme_linedraw()+ 
  xlab (sprintf ("PC1, %.1f%%", pcaPercentDev[1])) + 
  ylab (sprintf ("PC2, %.1f%%", pcaPercentDev[2]))
ggsave("PCA1.pdf", width = 8, height = 6)

