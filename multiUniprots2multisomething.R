library(data.table)
source(file.path("/Users/yzhou/Downloads",  "bp_utils" ,"UniprotIDMapping.R"))
options(timeout=1000)

multiUniprots2multisomething <- function (uniprots, sep = ";", something = "ENSEMBL", species = "HUMAN", simplify = FALSE, useDatFile = FALSE, allowDups = FALSE){
  toGenes <- data.table(uniprots = uniprots)
  toGenes <- toGenes[,.(singleUniprot = unlist(strsplit(uniprots, sep))),by = uniprots]
  toGenes[,singleGene := translateUniprot2Something(singleUniprot, something = something, species = species, useDatFile = useDatFile)]
  toGenes <- toGenes[!is.na(singleGene), ]
  if (simplify == TRUE){
    simplify = function(x)unique(sort(x))
  }else if (simplify == FALSE){
    simplify = identity # do nothing
  }else if (! "function" %in% class (simplify)){
    stop("unexpected simplify format")
  }
  toGenes <- toGenes[, .(genes = paste(simplify (singleGene), collapse=sep)), by = uniprots]
  if(!allowDups){
    duplicatedGeneNames <- unique(toGenes[duplicated(genes)]$genes)
    toGenes[genes %in% duplicatedGeneNames, genes := paste0(genes, ".", uniprots)]
  }
  
  setkey(toGenes, uniprots)
  
  return(toGenes[uniprots, genes])
}