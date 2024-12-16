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
# change the last _ to -
spc_long$BioReplicate <- sub("_(?=[^_]*$)", "-", spc_long$variable, perl=TRUE)
spc_long$Condition <- gsub("\\-.*", "", spc_long$BioReplicate)
# spc_long$BioReplicate <- gsub(".*_", "", spc_long$variable)
# spc_long$Condition <- sub('_[^_]*$', '', spc_long$variable)
# spc_long$BioReplicate <- paste(spc_long$Condition, spc_long$BioReplicate, sep = "-")
# check contaminate and remove
if(length(contaminate) == 1){
  if(length(which(spc_long$PROTID == contaminate)) != 0){
    spc_long <- spc_long[-which(spc_long$PROTID == contaminate), ]
  }
}else{
  if(length(which(spc_long$PROTID %in% contaminate)) != 0){
    spc_long <- spc_long[-which(spc_long$PROTID %in% contaminate), ]
  }
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
# change the last _ to -
int_long$BioReplicate <- sub("_(?=[^_]*$)", "-", int_long$variable, perl=TRUE)
int_long$Condition <- gsub("\\-.*", "", int_long$BioReplicate)
# check contaminate and remove
if(length(contaminate) == 1){
  if(length(which(int_long$PROTID == contaminate)) != 0){
    int_long <- int_long[-which(int_long$PROTID == contaminate), ]
  }
}else{
  if(length(which(int_long$PROTID %in% contaminate)) != 0){
    int_long <- int_long[-which(int_long$PROTID %in% contaminate), ]
  }
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

#CompPASS
# preys=read.table(file.choose(),header=F,sep="\t")
# baits=read.table(file.choose(),header=F,sep="\t")
interactions=read.table(file="/PATH_TO_PROJECT_DIR/msspc/spectral_counts-saint-interactions.txt", header=F,sep="\t",stringsAsFactors = F)
interactions <- reshape2::dcast(interactions, V1 + V2 ~ V3, value.var = "V4", fill = 0)
interactions <- reshape2::melt(interactions, id.vars = c("V1", "V2"), measure.vars = colnames(interactions)[3:ncol(interactions)])
interactions_agg <- aggregate(interactions$value, by=list(V2 = interactions$V2, variable = interactions$variable), FUN=sum)
interactions_agg <- interactions_agg[which(interactions_agg$x != 0), ]
interactions_agg$V2_variable <- paste(interactions_agg$V2, interactions_agg$variable, sep = "-")
interactions$V2_variable <- paste(interactions$V2, interactions$variable, sep = "-")
interactions <- interactions[which(interactions$V2_variable %in% interactions_agg$V2_variable), ]
interactions <- interactions[, -which(colnames(interactions) == "V2_variable")]

#R shiny
input=interactions
colnames(input)=c("Replicate","Bait","Prey","Spectral.Count")
Experiment.ID=input$Bait
input=cbind(Experiment.ID,input)
input$Prey <- as.character(input$Prey)
write.table(input,paste("/PATH_TO_PROJECT_DIR", "/CompPASS_Rshiny_input.txt",sep = ""),sep="\t",row.names = F, quote = F)

library(cRomppass)
SPCwStats.results.noNaN <- comppass(input, stats = NULL, norm.factor = 0.98)

#it's just a few lines of code to do the percentile by bait analysis on the compPASS output file.
#Note, sometimes you get NaN values in the compPASS output file, these need to be replaced with Inf prior to running these few lines of code.
library(dplyr)
#SPCwStats.results.noNaN=read.table(file="/Users/yzhou/Downloads/Roche/JGZ04/RSV-P/comppass_scores.tsv",header=T,sep="\t")
SPCwStats.results.noNaN[is.na(SPCwStats.results.noNaN)]=Inf
# Create columns that give the Percentile for each WD and Z score
SPCwStats.results.noNaN.wide = SPCwStats.results.noNaN %>% mutate(wd_percentile = percent_rank(WD), z_percentile = percent_rank(Z) ) %>%
  
  # group by the Baits so we can get the Percentile PER Bait for WD and Z scores
  group_by(Bait) %>%
  mutate(wd_percentile_perBait = percent_rank(WD), z_percentile_perBait = percent_rank(Z) ) %>%
  data.frame(stringsAsFactors=F)

write.table(SPCwStats.results.noNaN.wide, paste("PATH_TO_PROJECT_DIR", "/spc_comppass_scores_percentile.txt",sep = ""), quote=F, row.names=F, sep='\t')












