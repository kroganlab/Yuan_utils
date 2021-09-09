library(data.table)
library(dplyr)
library(tidyr)
MaxQtoMSstats <- function(evidencepath , keypath)
{
  evidenceDF <- fread(evidencepath, integer64 = "double")
  keyDF <- fread(keypath)
  
  # format evidence to be compatible for MSstats
  evidenceDF <- evidenceDF %>%
    # rename column to MSSstat format
    setNames(nm = make.names(names(.))) %>%
    dplyr::rename(RawFile         = Raw.file,
                  PeptideSequence = Sequence,
                  Sequence        = Modified.sequence) %>%
    # add missing column
    mutate(IsotopeLabelType = "L")
  
  keyDF$BioReplicate <- gsub("\\.", "-", keyDF$BioReplicate)
  
  # merge evidence file and keys file
  evidenceDF <- merge(x = evidenceDF, y = keyDF, by = c("RawFile", "IsotopeLabelType"))
  # remove missing Protein ids
  evidenceDF <- evidenceDF %>%
    filter(!is.na(Proteins))
  # remove protein groups
  evidenceDF <- evidenceDF %>%
    filter(!grepl(pattern = ";", Proteins))
  # remove contaminants 
  evidenceDF <- evidenceDF %>%
    filter(!grepl(pattern = "CON__|REV__", Proteins))
  
  # Prepare MSstats input file
  evidenceWideDF <- dcast.data.table( Proteins + Sequence + Charge ~ RawFile,
                                      data=data.table(evidenceDF), 
                                      value.var='Intensity', 
                                      fun.aggregate=sum, fill = NA)
  
  peptideDF <- evidenceWideDF %>%
    gather(RawFile, Intensity, -Proteins, -Sequence, -Charge) %>%
    merge(y = mutate(keyDF, RawFile = RawFile),
          by = "RawFile") %>%
    group_by(Sequence, Charge, IsotopeLabelType, Run, BioReplicate,
             Condition, Proteins) %>%
    summarize(Intensity = sum(Intensity), .groups = "keep") %>%
    dplyr::rename(ProteinName = Proteins,
                  PrecursorCharge = Charge,
                  PeptideSequence = Sequence) %>%
    mutate(ProductCharge = NA,
           FragmentIon = NA)
  
  peptideDF$PeptideSequence <- factor(peptideDF$PeptideSequence)
  peptideDF$ProteinName <- factor(peptideDF$ProteinName)
  #colnames(peptideDF)[which(colnames(peptideDF) == "IsotypeLabelType")] <- "IsotopeLabelType"
  
  write.table(peptideDF, "MSstats_intput.txt", row.names = F, quote = F, sep = "\t")
  
  return(peptideDF)
  
}

library(MSstats)
library(limma)
evidencepath <- "/Users/yzhou/Downloads/PCMI mutant/mutant/2021Apr_allMut_evidence.txt"
keypath <- "/Users/yzhou/Downloads/PCMI mutant/2021May_GRIN2B_mut_keys.txt"
evidenceDF <- fread(evidencepath, integer64 = "double")
keyDF <- fread(keypath)
outdir <- "/Users/yzhou/Downloads/PCMI mutant/Quantification"
contrastDF <- read.table(file = "/Users/yzhou/Downloads/PCMI mutant/GRIN2B_contrast.txt", header = F, sep = "\t")
ComparisonResult <- data.frame()

for(i in 1:dim(contrastDF)[1])
{
  contrast <- as.character(contrastDF[i])
  output_dir <- paste(outdir, paste(contrast, collapse = " vs "), sep = "/")
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  setwd(output_dir)
  keyDF_sub <- keyDF[which(keyDF$Condition %in% contrast), ]
  evidenceDF_sub <- evidenceDF[evidenceDF$`Raw file` %in% keyDF_sub$RawFile, ]
  write.table(keyDF_sub, paste("key", paste(contrast, collapse = " vs ", ".txt", sep = "")), row.names = F, quote = F, sep = "\t")
  write.table(evidenceDF_sub, paste("evidence", paste(contrast, collapse = " vs ", ".txt", sep = "")), row.names = F, quote = F, sep = "\t")
  peptideDF <- MaxQtoMSstats(evidencepath =  paste("evidence", paste(contrast, collapse = " vs ", ".txt", sep = "")), keypath = paste("key", paste(contrast, collapse = " vs ", ".txt", sep = "")))
  
  numThreads <- parallel::detectCores() - 4
  mssquant <- dataProcess(raw = as.data.frame(peptideDF),
                          featureSubset = "highQuality",
                          remove_uninformative_feature_outlier = TRUE,
                          #normalization = FALSE, # rely on complete case normalization above
                          clusters = numThreads, MBimpute=F)
  
  # contrasts
  modelMat <- model.matrix(~0+Condition, data = keyDF_sub) %>%
    as.data.frame() %>%
    setNames(nm = gsub(pattern = "Condition", replacement = "", names(.)))
  rownames(modelMat) <- paste0(keyDF_sub$RawFile, "_", keyDF_sub$IsotopeLabelType)
  compares <- paste(contrast, collapse = "-", sep = "")
  contrastMat <- makeContrasts(contrasts = compares, levels = modelMat)
  
  msstatsDE <- groupComparison(data = mssquant,
                               contrast.matrix = t(contrastMat))
  ########################################################
  #  Add replicate information for further post processing
  ########################################################
  # count and store replicates
  repCounts <- msstatsDE$ModelQC[!is.na(ABUNDANCE), 
                                 .(numBio = length(unique(SUBJECT_ORIGINAL)), 
                                   numTotal = length(unique(RUN))),
                                 by=.(PROTEIN, GROUP_ORIGINAL)]
  #unpack the contrast into group names
  setDT(msstatsDE$ComparisonResult)
  msstatsDE$ComparisonResult[,c("groupPos", "groupNeg") := tstrsplit(Label, split="-")]
  #merge counts for the positive condition
  setnames(repCounts, c("PROTEIN", "GROUP_ORIGINAL","bioRepPos", "totalRepPos"))
  msstatsDE$ComparisonResult <- merge (msstatsDE$ComparisonResult, 
                                       repCounts, by.x=c("Protein", "groupPos"), 
                                       by.y=c("PROTEIN", "GROUP_ORIGINAL"), all.x=TRUE)
  #merge counts for the negative condition
  setnames(repCounts, c("PROTEIN", "GROUP_ORIGINAL","bioRepNeg", "totalRepNeg"))
  msstatsDE$ComparisonResult <- merge (msstatsDE$ComparisonResult, 
                                       repCounts, by.x=c("Protein", "groupNeg"), 
                                       by.y=c("PROTEIN", "GROUP_ORIGINAL"), all.x=TRUE)
  
  
  write.csv(msstatsDE$ComparisonResult, paste("result", paste(contrast, collapse = " vs ", ".csv", sep = "")), row.names = F, quote = F)
  
  ComparisonResult <- rbind(ComparisonResult, msstatsDE$ComparisonResult)
  
  file.remove("msstats.log", "sessionInfo.txt")
}
setwd(outdir)
write.csv(ComparisonResult, "ComparisonResult.csv", row.names = F, quote = F)
