evidence <- read.table(file = "/Users/yzhou/Downloads/CoV01_HEK293T-SARS2_Nature-paper2/RAW_MaxQuant/combined/txt/evidence.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)
keys <- read.table(file = "/Users/yzhou/Downloads/CoV01_HEK293T-SARS2_Nature-paper2/keys_correct-order.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)

# CoV1
# heatmap for all viral protein
evidence_bait <- evidence[grep("SARS2_", evidence$Proteins), ]
evidence_bait <- merge(evidence_bait, keys[, c("RawFile", "BioReplicate")], by.x = "Raw file", by.y = "RawFile", all.x = T, all.y = F)
agg_intensity <- aggregate(evidence_bait$Intensity, by=list(RawFile = evidence_bait$`Raw file`, Experiment = evidence_bait$Experiment, Proteins=evidence_bait$Proteins, BioReplicate = evidence_bait$BioReplicate), FUN=sum, na.rm=TRUE, na.action=NULL)
agg_intensity <- agg_intensity[order(agg_intensity$RawFile), ]
agg_intensity$x <- log2(agg_intensity$x)
agg_intensity$BioReplicate <- paste(agg_intensity$RawFile, agg_intensity$BioReplicate, sep = "-")
agg_intensity <- agg_intensity[grep("ex", agg_intensity$RawFile, invert = T), ]
library(ggplot2)
ggplot(agg_intensity, aes(factor(BioReplicate, levels = unique(BioReplicate)), Proteins, fill= x)) + 
  geom_tile(color = "white")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  labs(title="Heatmap for all viral proteins", x ="", y = "", fill = "log2(Intensity)")+ 
  scale_fill_gradient(low="yellow",high="red") + 
  coord_equal()
setwd("/Users/yzhou/Downloads/CoV01_HEK293T-SARS2_Nature-paper2")
ggsave("Heatmap for all viral proteins.pdf", width = 16, height = 8)
