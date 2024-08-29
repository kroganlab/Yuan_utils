library (data.table)
library (ComplexHeatmap)
library (ggplot2)
library(readxl)
spectronautPeptideFile <- "/PATH_TO_PROJECT_DIR/MSstats.csv"

a <- fread (spectronautPeptideFile)

# histogram for log2 intensity
hist(log2(a$Intensity), breaks = 100)
abline(v = 5)

source (file.path("~/Downloads", "bp_utils", "spectronautFile2ArtMS.R"))

cf<- list()
# normalization method FALSE = no normalization; default is global medians which you can get my removing/commenting out all normalization lines
# cf$msstats$normalization_method = FALSE

#cf$msstats$normalization_method = "globalStandards"
#cf$msstats$normalization_reference <-  "P38398"

# should artms attempt to annotate proteins 1 = yes; 0 = no
cf$output_extras$annotate$enabled <- as.integer(1)
# should artms do extended QC 1 = yes; 0= no
cf$qc$extended <- as.integer(1)
cf$qc$basic <- as.integer(1)

# cf$output_extras$annotate$species <- "MOUSE"

# make files in artMS format
globalInput <- spectronautFile2ArtMS(spectronautPeptideFile, 
                                     outFilePrefix = "/PATH_TO_PROJECT_DIR", 
                                     artmsConfig = cf, contrastPatterns  = contrastPatterns)
evidence <- read.table(file = "/PATH_TO_PROJECT_DIR/evidence.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)
evidence_sub <- evidence[-which(is.na(evidence$Intensity)), ]
# check contaminate
contaminate <- c("O77727", "P00698", "P00761", "P00883", "P02662", "P02663", "P02666", "P02668", "P02769")
# check Leading proteins formate
if(any(grepl("sp\\|", evidence_sub$`Leading proteins`)))
{
  #evidence_sub$`Leading proteins` <- gsub("sp\\|", "", evidence_sub$`Leading proteins`)
  # evidence_sub$`Leading proteins`[grep("_HUMAN", evidence_sub$`Leading proteins`, invert = T)] <- paste("CON__", evidence_sub$`Leading proteins`[grep("_HUMAN", evidence_sub$`Leading proteins`, invert = T)], sep = "")
  #evidence_sub$`Leading proteins` <- gsub("\\|.*", "", evidence_sub$`Leading proteins`)
  
  evidence_sub$`Leading proteins` <- gsub('([a-z,0-9,A-Z,\\_]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                                   '\\2',
                                   evidence_sub$`Leading proteins`)
  
}
if(any(contaminate %in% evidence_sub$`Leading proteins`))
{
  evidence_sub$`Leading proteins`[which(evidence_sub$`Leading proteins` %in% contaminate)] <- 
    paste("CON__", evidence_sub$`Leading proteins`[which(evidence_sub$`Leading proteins` %in% contaminate)], sep = "")
}
write.table(evidence_sub, "/PATH_TO_PROJECT_DIR/evidence_sub.txt", sep = "\t", row.names = F, quote = F)

# QC
setwd("/PATH_TO_PROJECT_DIR")
artmsQualityControlEvidenceBasic(evidence_file = "/PATH_TO_PROJECT_DIR/evidence_sub.txt", 
                                 keys_file = "/PATH_TO_PROJECT_DIR/keys.txt", 
                                 prot_exp = "APMS")
artmsQualityControlEvidenceExtended(evidence_file = "/PATH_TO_PROJECT_DIR/evidence_sub.txt", 
                                    keys_file = "/PATH_TO_PROJECT_DIR/keys.txt", 
                                    plotPCA = FALSE)

# PCA
# MSstats
peptideDF <- fread(spectronautPeptideFile)
peptideDF$PeptideModifiedSequence <- peptideDF$PeptideSequence
peptideDF$IsotopeLabelType <- "L"

library(MSstats)
mssquant <- dataProcess(raw = as.data.frame(peptideDF),
                        MBimpute=F)
write.table(
  mssquant$FeatureLevelData,
  file = 'output/mss-FeatureLevelData.txt',
  eol = "\n",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
outputDirectory <- "/PATH_TO_PROJECT_DIR/output"
mss.normalized.txt <-file.path(outputDirectory, "mss-FeatureLevelData.txt") 
normPepInt <- fread (mss.normalized.txt)

normPepInt[, logCenteredIntensity := log2(INTENSITY/(median(INTENSITY, na.rm=TRUE))), by = PEPTIDE]
normInt.mat <- as.matrix(dcast(normPepInt, PEPTIDE~GROUP+SUBJECT, value.var = "logCenteredIntensity"), rownames = "PEPTIDE")

# subset to complete cases
normInt.mat <- normInt.mat[complete.cases(normInt.mat),]  # select rows with no missing values

colInfo <- data.table(colname = colnames(normInt.mat))
# something like, this depends on the structure of your condition names
# colInfo[,c("treat", "time", "rep") := tstrsplit(colname, "[_.]", keep = c(1,2,5)) ]
colInfo$treat <- gsub("_.*", "", colInfo$colname)
# colInfo$rep <- gsub(".*-", "", colInfo$colname)
colInfo$rep <- gsub(".*_", "", colInfo$colname)


title <- NULL
#PCA
pcaOut <- prcomp(t(normInt.mat))
pcaDT <- as.data.table(pcaOut$x, keep.rownames=TRUE)

pcaPercentVar <- round(100 * (pcaOut$sdev^2)/sum(pcaOut$sdev^2), 1)

pcaDT <- merge (pcaDT, colInfo, by.x = "rn", by.y = "colname", all.x = TRUE)

#plot first two components
p <- ggplot (pcaDT, aes(x=PC1, y=PC2, fill = treat, shape = rep)) + 
  geom_point(alpha=1.0, size=4) + 
  ggrepel::geom_text_repel(aes(label=rn), show.legend = FALSE, size = 3) +
  theme_bw() + 
  xlab (sprintf ("PC1, %.1f%%", pcaPercentVar[1])) + 
  ylab (sprintf ("PC2, %.1f%%", pcaPercentVar[2])) + 
  ggtitle (sprintf ("PCA %s using %d features (log intensity)", title, nrow(normInt.mat))) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  scale_shape_manual(values = 21:25) +
  #scale_fill_manual(values = c(`05` = "gray", `30` = "black")) +
  guides(fill = guide_legend(override.aes = list(shape =21) ) ,
         color = guide_legend(override.aes = list(shape =21) ) )
# BackupAsPDF(p, "PCA_Complete_Features")
cairo_pdf(file.path(outputDirectory, "PCA_Complete_Features.pdf"), width = 10, height = 7)
print (p)
dev.off()

# if you want to exclude runs, put them in this vector
# it is the "BioReplicate" column in the keys.txt file
runsToSkip <- c() #c("MDA_Control-1")


# timsTOF to SAINTexpress
spc <- data.table::fread(file = "/PATH_TO_PROJECT_DIR/reprint.spc.tsv")
# remove first row
spc <- spc[-1, ]
spc_long <- reshape2::melt(spc, id.vars = c("PROTID", "GENEID", "PROTLEN"))
spc_long$variable <- as.character(spc_long$variable)
spc_long$value <- as.integer(spc_long$value)
spc_long$variable <- gsub("_SPC", "", spc_long$variable)
spc_long$BioReplicate <- gsub(".*_", "", spc_long$variable)
spc_long$Condition <- sub('_[^_]*$', '', spc_long$variable)
spc_long$BioReplicate <- paste(spc_long$Condition, spc_long$BioReplicate, sep = "-")
# check contaminate and remove
if(length(contaminate) == 1){
  spc_long <- spc_long[-which(spc_long$PROTID == contaminate), ]
}else{
  spc_long <- spc_long[-which(spc_long$PROTID %in% contaminate), ]
}
# remove samples didn't pass QC
if(length(runsToSkip) > 0){
  if(length(runsToSkip) == 1){
    spc_long <- spc_long[-which(spc_long$BioReplicate == runsToSkip), ]
  }else{
    spc_long <- spc_long[-which(spc_long$BioReplicate %in% runsToSkip), ]
  }
  
}

# Bait
dir.create(paste("/PATH_TO_PROJECT_DIR", "/msspc", sep = ""))
setwd(paste("/PATH_TO_PROJECT_DIR", "/msspc", sep = ""))
Bait <- unique(spc_long$BioReplicate)
Bait <- gsub("-", "_", Bait)
Bait <- sub("_([^_]*)$", "-\\1", Bait)
Bait <- data.frame(BioReplicate = Bait, Condition = gsub("-.*", "", Bait))
Bait$SAINT <- "T"
Bait$SAINT[grep("Control", Bait$Condition)] <- "C"
write.table(Bait, "spectral_counts-saint-baits.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# interactions
interactions <- data.frame(BioReplicate = spc_long$BioReplicate, 
                           Condition = gsub("-.*", "", spc_long$BioReplicate), 
                           PROTID = spc_long$PROTID, 
                           spc = spc_long$value)
interactions$BioReplicate <- gsub("-", "_", interactions$BioReplicate)
interactions$BioReplicate <- sub("_([^_]*)$", "-\\1", interactions$BioReplicate)
interactions$Condition <- gsub("-.*", "", interactions$BioReplicate)

# remove 0 counts
interactions <- interactions[which(interactions$spc != 0), ]
write.table(interactions, "spectral_counts-saint-interactions.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# preys
preys <- data.frame(PROTID = spc_long$PROTID, 
                    PROTLEN = spc_long$PROTLEN, 
                    GENEID = spc_long$GENEID)
preys <- unique(preys)
write.table(preys, "spectral_counts-saint-preys.txt", sep = "\t", quote = F, row.names = F, col.names = F)

##########################
int <- data.table::fread(file = "/PATH_TO_PROJECT_DIR/reprint.int.tsv")
# remove first row
int <- int[-1, ]
int_long <- reshape2::melt(int, id.vars = c("PROTID", "GENEID"))
int_long$variable <- as.character(int_long$variable)
int_long$value <- as.numeric(int_long$value)
int_long$variable <- gsub("_INT", "", int_long$variable)
int_long$BioReplicate <- gsub(".*_", "", int_long$variable)
int_long$Condition <- sub('_[^_]*$', '', int_long$variable)
int_long$BioReplicate <- paste(int_long$Condition, int_long$BioReplicate, sep = "-")
# check contaminate and remove
if(length(contaminate) == 1){
  int_long <- int_long[-which(int_long$PROTID == contaminate), ]
}else{
  int_long <- int_long[-which(int_long$PROTID %in% contaminate), ]
}
# remove samples didn't pass QC
if(length(runsToSkip) > 0){
  if(length(runsToSkip) == 1){
    int_long <- int_long[-which(int_long$BioReplicate == runsToSkip), ]
  }else{
    int_long <- int_long[-which(int_long$BioReplicate %in% runsToSkip), ]
  }
  
}
# Bait
dir.create(paste("/PATH_TO_PROJECT_DIR", "/msint", sep = ""))
setwd(paste("/PATH_TO_PROJECT_DIR", "/msint", sep = ""))
Bait <- unique(int_long$BioReplicate)
Bait <- gsub("-", "_", Bait)
Bait <- sub("_([^_]*)$", "-\\1", Bait)
Bait <- data.frame(BioReplicate = Bait, Condition = gsub("-.*", "", Bait))
Bait$SAINT <- "T"
Bait$SAINT[grep("Control", Bait$Condition)] <- "C"
write.table(Bait, "Intensity-saint-baits.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# interactions
interactions <- data.frame(BioReplicate = int_long$BioReplicate, 
                           Condition = gsub("-.*", "", int_long$BioReplicate), 
                           PROTID = int_long$PROTID, 
                           int = int_long$value)
interactions$BioReplicate <- gsub("-", "_", interactions$BioReplicate)
interactions$BioReplicate <- sub("_([^_]*)$", "-\\1", interactions$BioReplicate)
interactions$Condition <- gsub("-.*", "", interactions$BioReplicate)

# remove 0 counts
interactions <- interactions[which(interactions$int != 0), ]
write.table(interactions, "Intensity-saint-interactions.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# preys
preys <- data.frame(PROTID = int_long$PROTID,  
                    GENEID = int_long$GENEID)
preys <- unique(preys)
file.copy(from="/PATH_TO_PROJECT_DIR/msspc/spectral_counts-saint-preys.txt", 
          to="/PATH_TO_PROJECT_DIR/msint", 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)















