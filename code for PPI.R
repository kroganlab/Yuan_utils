library(artMS)
library(seqinr)#read.fasta
library(readxl)
outdir <- "/Users/yzhou/Downloads/A3_PPI_Interactome"
evidence <- read.table(file = "/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/evidence.txt", header = T, sep = "\t", stringsAsFactors = F)

# sample_annotation <- read_excel("/Users/yzhou/Downloads/IBV/QX2/Vir Collaboration Log.xlsx")
# 
# a <- unique(evidence[, c("Raw.file", "Experiment")])
# sample_a <- merge(a, sample_annotation, by.x = "Raw.file", by.y = "File Name")
# sample_a <- sample_a[, c("Experiment", "Sample Name")]
# # sample_a$`Sample number` <- matrix(unlist(strsplit(sample_a$`Sample number`, split = " ")), ncol = 2, byrow = T)[, 2]
# # all(sample_a$Experiment == sample_a$`Sample number`)
# 
# 
# a <- unique(evidence[, c("Raw.file", "Experiment")])
# keys <- merge(a, sample_annotation, by.x = "Raw.file", by.y = "File Name")
# keys <- keys[, c("Raw.file", "Bait Name")]
# colnames(keys) <- c("RawFile", "Condition")
# keys$Condition <- gsub("-", "_", keys$Condition)
# keys$Condition <- gsub("/", "_", keys$Condition)
# keys$Condition <- gsub(" ", "_", keys$Condition)
# # keys$Condition[which(keys$Condition == "Empty_vector_1" | keys$Condition == "Empty_vector_2")] <- "Empty_vector"
# # keys$Condition[which(keys$Condition == "GFP_strep_1" | keys$Condition == "GFP_strep_2")] <- "GFP_strep"
# keys$IsotopeLabelType <- "L"
# keys$BioReplicate <- keys$Condition
# for(g in unique(keys$Condition)){
#   len <- sum(keys$Condition == g)
#   ind <- which(keys$Condition == g)
#   keys$BioReplicate[ind] <- paste(keys$BioReplicate[ind], seq(1,len), sep = "-")
# }
# keys$Run <- seq(1, dim(keys)[1])
# keys$SAINT <- "T"
# keys$SAINT[which(keys$Condition == "Empty_vector" | keys$Condition == "GFP_strep")] <- "C"
# write.table(keys, paste(outdir, "keys.txt", sep = "/"), sep = "\t", row.names = F, quote = F)


# remove VE20160521-75 (GFP_RNAse-3) : failed QC
# keys <- read.table(file = "/Users/yzhou/Downloads/A3_PPI_Interactome/Keys_File_noGFPRNAse_3.txt", header = T, sep = "\t", stringsAsFactors = F)

# a <- unique(evidence[, c("Raw.file", "Experiment")])
# b <- keys[, c("RawFile", "BioReplicate")]
# c <- merge(a, b, by.x = "Raw.file", by.y = "RawFile")


#QC
dir.create(paste(outdir, "/QC Extended", sep = ""))
setwd(paste(outdir, "/QC Extended", sep = ""))
artmsQualityControlEvidenceExtended(evidence_file="/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/evidence.txt",
                                    keys_file="/Users/yzhou/Downloads/A3_PPI_Interactome/Keys_File_noGFPRNAse_3.txt")

dir.create(paste(outdir, "/QC Basic", sep = ""))
setwd(paste(outdir, "/QC Basic", sep = ""))
artmsQualityControlEvidenceBasic(evidence_file="/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/evidence.txt",
                                 keys_file="/Users/yzhou/Downloads/A3_PPI_Interactome/Keys_File_noGFPRNAse_3.txt",
                                 prot_exp="APMS",output_name="A3", verbose = TRUE)


dir.create(paste(outdir, "/msspc", sep = ""))
setwd(paste(outdir,"/msspc", sep = ""))
artmsEvidenceToSaintExpress(evidence_file="/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/evidence.txt",
                            keys_file="/Users/yzhou/Downloads/A3_PPI_Interactome/Keys_File_noGFPRNAse_3.txt",
                            ref_proteome_file="/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/SwissProt.Human.canonical.2018.10.09.fasta",
                            quant_variable="msspc",output_file=paste("A3", "_spectral_counts.txt", sep = ""), verbose = TRUE)

dir.create(paste(outdir, "/msint", sep = ""))
setwd(paste(outdir, "/msint", sep = ""))
artmsEvidenceToSaintExpress(evidence_file="/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/evidence.txt",
                            keys_file="/Users/yzhou/Downloads/A3_PPI_Interactome/Keys_File_noGFPRNAse_3.txt",
                            ref_proteome_file="/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/SwissProt.Human.canonical.2018.10.09.fasta",
                            quant_variable="msint",output_file=paste("A3", "_MS_Intensity.txt", sep = ""), verbose = TRUE)


#MIST
my.load_evidencekey=function(evidence,key)
{
  data=evidence
  colnames(data) <- gsub('\\s','.',colnames(data))
  data <- subset(data, trimws(data$Proteins) != "") # remove white Proteins ids
  colnames(data)[grep(colnames(data), pattern="raw.file", ignore.case = TRUE)] <- "RawFile"
  data$Intensity[data$Intensity<1] <- NA
  if(!'IsotopeLabelType' %in% colnames(data)) data$IsotopeLabelType <- 'L'
  
  keys <- key
  
  # check
  unique_data <- unique(data$RawFile)
  unique_keys <- unique(keys$RawFile)
  keys_not_found = setdiff(unique_keys, unique_data)
  data_not_found = setdiff(unique_data, unique_keys)
  if (length(keys_not_found) > 0) {
    message(sprintf("Keys: %d/%d RawFiles not found: %s", length(keys_not_found), length(unique_keys), paste(keys_not_found, collapse=",")))
  }
  if (length(data_not_found) > 0) {
    message(sprintf("Data: %d/%d RawFiles not found: %s", length(data_not_found), length(unique_data), paste(data_not_found, collapse=",")))
  }
  
  # combine
  data = merge(data, keys, by=c('RawFile','IsotopeLabelType'))
  return(list(data=data,keys=keys))
}

filterMaxqData <- function(data){
  data_selected <- data[grep("CON__|REV__",data$Proteins, invert=T),]
  blank.idx <- which(data_selected$Proteins =="")
  if(length(blank.idx)>0)  data_selected = data_selected[-blank.idx,]
  return(data_selected)
}

filterData <- function(data){
  data_f = data
  #if(config$filters$protein_groups == 'remove') 
  data_f <- data_f[grep(";",data_f$Proteins, invert=T),]
  #if(config$filters$contaminants) 
  data_f <- filterMaxqData(data_f)
  msg <- sprintf("Data Filter: %d/%d (%s%%) records remained", 
                 nrow(data_f), nrow(data), round(100*nrow(data_f)/nrow(data),1))
  message(msg)
  return(data_f)
}

# ******************************************** #
# FUNCTIONS for AP-MS scoring
# ******************************************** #

# Filtering and aggregate intensity/spectral_count for each protein
## rm_itself: remove interactions of bait and itself
## fix_itself: for interactions of bait and itself, fix prey name
my.preprocessAPMS <- function(evidence,key, rm_itself = TRUE, fix_itself = TRUE) {
  # Load data/key
  #config <- my.load_config(dataset_id)
  data <- my.load_evidencekey(evidence,key)$data
  
  # Filter
  data_f <- data
  # if (!is.null(config$filters$resolve_baitgroup)) {
  #   data_f <- resolve_baitgroup(data=data_f, bait_pattern=trimws(config$filters$resolve_baitgroup))
  #   data_f <- remove_bait_contaminant(data_f, bait_pattern=trimws(config$filters$resolve_baitgroup))
  # }
  data_f <- filterData(data_f)
  colnames(data_f)[grep(pattern="ms.ms.count", x = colnames(data_f), ignore.case = TRUE)] <- 'spectral_counts'
  
  if (rm_itself) data_f <- subset(data_f, Condition != Proteins)
  if (!rm_itself & fix_itself) {
    idx <- which(data_f$Condition == data_f$Proteins)
    data_f[idx,"Proteins"] <- paste0(data_f[idx,"Proteins"], "prey")
  }
  
  # Aggregate Intensity
  data_f_agg <- aggregate(Intensity ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins+Sequence+Charge, data=data_f, FUN = max)
  data_f_agg <- aggregate(Intensity ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins, data=data_f_agg, FUN = sum)
  data_f_agg <- subset(data_f_agg, !is.na(Intensity))
  
  # Aggregate SPC
  data_f_spc <- aggregate(spectral_counts ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins+Sequence+Charge,data=data_f,FUN = max)
  data_f_spc <- aggregate(spectral_counts ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins,data=data_f_spc,FUN = sum)
  data_f_spc <- subset(data_f_spc, !is.na(spectral_counts))
  
  return( list(data_f = data_f, agg_intensity = data_f_agg, agg_spc = data_f_spc) )
}

my.MaxQToMIST <- function(evidence,key,ref_proteome_file, outdir = "/bf2/smb-data/tnguyen/projects/fluomics/tempdata", quant_variable="spc") {
  # Load and aggregate data/key
  datalist <- my.preprocessAPMS(evidence,key)
  
  quant_variable <- trimws(tolower(quant_variable))
  if (! quant_variable %in% c("spc","intensity")) stop("Please input quant_variable as 'spc' or 'intensity'")
  data_sel <- datalist$agg_spc
  data_sel$ms_spectral_counts <- data_sel$spectral_counts
  quant_col <- 'ms_spectral_counts'
  if (quant_variable == "intensity") {
    data_sel <- datalist$agg_intensity
    data_sel$ms_intensity <- data_sel$Intensity
    quant_col <- 'ms_intensity'
  }
  keysout <- unique(data_sel[,c("RawFile","Condition")])
  
  # Select columns
  data_sel$ms_unique_pep = ""
  data_sel <- data_sel[,c('RawFile','Proteins','ms_unique_pep', quant_col)]
  colnames(data_sel) = c('id','ms_uniprot_ac','ms_unique_pep', quant_col)
  
  # Uniprot annotate
  #conn <- my.bf2_connect()
  #uniprot <- dbGetQuery(conn, "select distinct Entry as ms_uniprot_ac, Length from view_uniprot")
  #dbDisconnect(conn)
  
  # library(dplyr)
  # uniprot <- distinct(uniprot,Entry,Length)
  # colnames(uniprot)[which(colnames(uniprot)=="Entry")]="ms_uniprot_ac"
  # d <- setdiff(data_sel$ms_uniprot_ac, uniprot$ms_uniprot_ac)
  # if (length(d)>0) {
  #   msg <- sprintf("These proteins are not found in uniprot db, please check: %s", paste(d, collapse=","))
  #   #stop(msg)
  #   message(msg)
  # }
  
  ref_proteome <- read.fasta(
    file = ref_proteome_file,
    seqtype = "AA",
    as.string = TRUE,
    set.attributes = TRUE,
    legacy.mode = TRUE,
    seqonly = FALSE,
    strip.desc = FALSE
  )
  p_lengths <- c()
  p_names <- c()
  for (e in ref_proteome) {
    p_lengths <- c(p_lengths, nchar(e[1]))
    p_names <- c(p_names, attr(e, 'name'))
  }
  ref_table <- data.table(names = p_names, lengths = p_lengths)
  ref_table[, uniprot_ac := gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                                 '\\2',
                                 names)]
  ref_table[, uniprot_id := gsub('([a-z,0-9,A-Z]+\\|{1})([a-z,0-9,A-Z]+\\|{1})([A-Z,a-z,0-9,_]+)',
                                 '\\3',
                                 names)]
  colnames(ref_table)[which(colnames(ref_table)=="uniprot_ac")]="ms_uniprot_ac"
  d <- setdiff(data_sel$ms_uniprot_ac, ref_table$ms_uniprot_ac)
  if (length(d)>0) {
    msg <- sprintf("These proteins are not found in fasta, please check: %s", paste(d, collapse=","))
    #stop(msg)
    message(msg)
  }
  
  # Get mass
  data_sel <- base::merge(data_sel, ref_table, by = "ms_uniprot_ac")
  data_sel$Mass <- 110*data_sel$lengths
  
  # Write out
  ## data
  outdir <- paste0(outdir, "/mist/", quant_variable)
  dir.create(outdir, show=FALSE, recursive = TRUE)
  data_outfile <- sprintf("%s/mist-data.%s.txt", outdir, quant_variable)
  write.table(data_sel, file=data_outfile, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  
  ## keys
  key_outfile <- sprintf("%s/mist-key.%s.txt", outdir, quant_variable)
  write.table(keysout, file=key_outfile, eol='\n', sep='\t', quote=F, row.names=F, col.names=F)
  
  return(list(data_file=data_outfile, keys_file=key_outfile))
}


setwd(outdir)
evidence <- read.table(file = "/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/evidence.txt", header = T, sep = "\t", stringsAsFactors = F)
keys <- read.table(file = "/Users/yzhou/Downloads/A3_PPI_Interactome/Keys_File_noGFPRNAse_3.txt", header = T, sep = "\t", stringsAsFactors = F)
#colnames(keys)[which(colnames(keys)=="Raw.file")] = "RawFile"
colnames(keys)[which(colnames(keys)=="SAINT")]="TestControl"
BaitName=keys$Condition
keys=cbind(keys,BaitName)
#uniprot <- read.delim(file="/Users/yzhou/Downloads/view_uniprot.txt",header=T,sep = "\t",stringsAsFactors=F)
library(seqinr)
library(data.table)
my.MaxQToMIST(evidence,keys,
              "/Users/yzhou/Downloads/A3_PPI_Interactome/MaxQuant_20200910/SwissProt.Human.canonical.2018.10.09.fasta", 
              outdir = outdir, quant_variable="spc")

library(yaml)
mist=read_yaml(file="/Users/yzhou/Downloads/Roche/mist.yaml")

mist$files$data="/Users/yzhou/Downloads/A3_PPI_Interactome/mist/spc/mist-data.spc.txt"
mist$files$keys="/Users/yzhou/Downloads/A3_PPI_Interactome/mist/spc/mist-key.spc.txt"
mist$files$remove="/Users/yzhou/Downloads/A3_PPI_Interactome/mist/spc/a.txt"
mist$files$collapse="/Users/yzhou/Downloads/A3_PPI_Interactome/mist/spc/a.txt"
mist$files$specificity_exclusions="/Users/yzhou/Downloads/A3_PPI_Interactome/mist/spc/a.txt"
mist$files$output_dir=paste(outdir, "/mist/spc/result", sep = "")

mist$preprocess$filter_contaminants=0
mist$preprocess$contaminants_file=NULL
mist$preprocess$id_colname="id"
mist$preprocess$prey_colname="ms_uniprot_ac"
mist$preprocess$pepcount_colname="ms_spectral_counts"
mist$preprocess$mw_colname="Mass"

mist$qc$matrix_file=NULL

mist$mist$matrix_file=NULL
mist$mist$weights="fixed"

mist$annotate$enabled=1
mist$annotate$species="HUMAN"
mist$annotate$uniprot_dir="/Users/yzhou/Downloads/uniprot-human-filtered-organism__Homo.tab"

mist$enrichment$enabled=0

write_yaml(mist,file=paste(outdir, "/mist/mist.yaml", sep = ""))

system(paste("Rscript /Users/yzhou/Downloads/private.mist/main.R -c ", paste(outdir, "/mist/mist.yaml", sep = ""), sep = ""))


#CompPASS
# preys=read.table(file.choose(),header=F,sep="\t")
# baits=read.table(file.choose(),header=F,sep="\t")
interactions=read.table(file="/Users/yzhou/Downloads/A3_PPI_Interactome/msspc/A3_spectral_counts-saint-interactions.txt",
                        header=F,sep="\t",stringsAsFactors = F)

#R shiny
input=interactions
colnames(input)=c("Replicate","Bait","Prey","Spectral.Count")
Experiment.ID=input$Bait
input=cbind(Experiment.ID,input)
write.table(input,paste(outdir, "/CompPASS_Rshiny_input.txt",sep = ""),sep="\t",row.names = F, quote = F)

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

write.table(SPCwStats.results.noNaN.wide, paste(outdir, "/spc_comppass_scores_percentile.txt",sep = ""), quote=F, row.names=F, sep='\t')
