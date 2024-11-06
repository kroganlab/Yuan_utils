# MSstats (normalizeByMedianPolish)
dir.create("/DIR_TO_PROJECT/Quantification_normalizeByMedianPolish")
setwd("/DIR_TO_PROJECT/Quantification_normalizeByMedianPolish")
library(MSstats)
evidence <- fread("/DIR_TO_PROJECT/evidence_sub.txt")
keys <- fread("/DIR_TO_PROJECT/keys_sub.txt")
keys$Run <- keys$RawFile
write.table(keys, "/DIR_TO_PROJECT/Quantification_normalizeByMedianPolish/keys.txt", sep = "\t", row.names = F, quote = F)
source("/Users/yzhou/Documents/GitHub/bp_utils/EvidenceFile2MSstatsInput.R")
MSstats_format <- prepareDataForMSStats(evidenceFile = "/DIR_TO_PROJECT/evidence_sub.txt", 
                                        keysFile = "/DIR_TO_PROJECT/Quantification_normalizeByMedianPolish/keys.txt", 
                                        outfile=NULL)
# Normalized by normalizeByMedianPolish
source("/Users/yzhou/Documents/GitHub/bp_utils/Normalization.R")
featureIntensitiesDT <- MSstats_format
featureIntensitiesDT$runID <- featureIntensitiesDT$Run
featureIntensitiesDT$proteinID <- featureIntensitiesDT$ProteinName
featureIntensitiesDT$featureID <- paste(featureIntensitiesDT$PeptideSequence, featureIntensitiesDT$PrecursorCharge, sep = "_")
featureIntensitiesDT$logIntensity <- log2(featureIntensitiesDT$Intensity)

standard_table <- data.frame(Condition = c("eL", "eNP", "eVP24", "eVP24_NES", "eVP30", "eVP35", "eVP40", 
                                           "mL", "mNP", "mVP24", "mVP30", "mVP35", "mVP40"), 
                             standards = c("AAD14589.1", "AAD14590.1", "AAD14588.1", "AAD14588.1", "AAM76036.1", "AAD14582.1", "AAD14583.1", 
                                           "ABA87130.2", "ABA87124.1", "CAA78119.1", "ABA87128.1", "ABA87125.1", "sp|P35260.3|VP40_MABVM"))
# remove GFP and Vector
C <- unique(unique(featureIntensitiesDT$Condition))
C <- C[-which(C %in% c("GFP", "Empty"))]
nomalized_featureIntensitiesDT <- data.frame()
for(c in C)
{
  featureIntensitiesDT_sub <- featureIntensitiesDT[which(featureIntensitiesDT$Condition == c), ]
  nomalized_featureIntensitiesDT_sub <- normalizeByMedianPolish(featureIntensitiesDT_sub, 
                                                                standards = standard_table$standards[which(standard_table$Condition == c)], 
                                                                doPlots = TRUE )
  
  nomalized_featureIntensitiesDT <- rbind(nomalized_featureIntensitiesDT, nomalized_featureIntensitiesDT_sub)
}

nomalized_featureIntensitiesDT$Intensity <- 2^nomalized_featureIntensitiesDT$normLogIntensity

# MSstats
nomalized_featureIntensitiesDT$PeptideModifiedSequence <- nomalized_featureIntensitiesDT$PeptideSequence
library(MSstats)
mssquant <- dataProcess(raw = as.data.frame(nomalized_featureIntensitiesDT),
                        MBimpute=F, normalization = F)
write.table(
  quantification(mssquant, type = "Sample"),
  file = 'mss-sampleQuant.txt',
  eol = "\n",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  mssquant$FeatureLevelData,
  file = 'mss-FeatureLevelData.txt',
  eol = "\n",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  mssquant$ProteinLevelData,
  file = 'mss-ProteinLevelData.txt',
  eol = "\n",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
# contrasts
contrastPatterns <- read.table(file = "/DIR_TO_PROJECT/Quantification/contrast.txt", header = F, sep = "\t", stringsAsFactors = F)

keys <- keys[-which(keys$Condition %in% c("Empty", "GFP")), ]

modelMat <- model.matrix(~0+Condition, data = keys) %>%
  as.data.frame() %>%
  setNames(nm = gsub(pattern = "Condition", replacement = "", names(.)))
rownames(modelMat) <- paste0(keys$RawFile, "_", keys$IsotopeLabelType)
# compares <- paste(contrast, collapse = "-", sep = "")
# compares <- paste(contrastDF[, 1], contrastDF[, 2], sep = "-")
contrastMat <- limma::makeContrasts(contrasts = contrastPatterns, levels = modelMat)

GroupComparison <- groupComparison(data = mssquant,
                                   contrast.matrix = t(contrastMat),
                                   use_log_file = FALSE)
write.table(
  GroupComparison$ComparisonResult,
  file = "result.txt",
  eol = "\n",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  GroupComparison$ModelQC,
  file = "ModelQC.txt",
  eol = "\n",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
