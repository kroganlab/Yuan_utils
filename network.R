# barplot
# Load score data
#NHBE
NHBE <- read.table(file = "/Users/yzhou/Downloads/Roche/Roche_SARS2/CoV20_NHBE-SARS2/result_filter.txt", header = T, sep = "\t", stringsAsFactors = F)
NHBE <- NHBE[grep("-qx", NHBE$Bait.x, invert = T), ]

#THP1
THP1 <- read.table(file = "/Users/yzhou/Downloads/Roche/Roche_SARS2/CoV25_THP1-SARS2/result_filter.txt", header = T, sep = "\t", stringsAsFactors = F)

#A549
A549 <- read.table(file = "/Users/yzhou/Downloads/Roche/Roche_SARS2/CoV27_A549-SARS2/result_filter.txt", header = T, sep = "\t", stringsAsFactors = F)

Homo <- read.delim2("/Users/yzhou/Downloads/uniprot-human-filtered-organism__Homo.tab")
Homo <- Homo[, c("Entry", "Gene.names...primary..")]

NHBE$Strain <- "NHBE"
THP1$Strain <- "THP1"
A549$Strain <- "A549"

result_filter <- rbind(NHBE, THP1, A549)
result_filter$BaitProtein <- result_filter$Bait.x


# Cytoscape
pchart<-'piechart: showlabels="false" range="0,1" arcstart="90" valuelist=".33,.33,.34" 
colorlist="up:#4292c6,zero:white,down:white;up:#807dba,zero:white,down:white;up:#41ab5d,zero:white,down:white" 
attributelist="NHBE,THP1,A549" labellist="NHBE,THP1,A549"'

result_filter$BaitProtein_Prey <- paste(result_filter$BaitProtein, result_filter$Prey.x, sep = "_")
sig_prey <- unique(result_filter$BaitProtein_Prey)

# read unthreshold score
result_NHBE <- read.table("/Users/yzhou/Downloads/Roche/Roche_SARS2/CoV20_NHBE-SARS2/result.txt", header = T, sep = "\t", stringsAsFactors = F)
result_NHBE$Strain <- "NHBE"

result_THP1 <- read.table("/Users/yzhou/Downloads/Roche/Roche_SARS2/CoV25_THP1-SARS2/result.txt", header = T, sep = "\t", stringsAsFactors = F)
result_THP1$Strain <- "THP1"
result_THP1 <- merge(result_THP1, Homo, by.x = "Prey.x", by.y = "Entry", all.x = T, all.y = F)

result_A549 <- read.table("/Users/yzhou/Downloads/Roche/Roche_SARS2/CoV27_A549-SARS2/result.txt", header = T, sep = "\t", stringsAsFactors = F)
result_A549$Strain <- "A549"
result_A549 <- merge(result_A549, Homo, by.x = "Prey.x", by.y = "Entry", all.x = T, all.y = F)

result <- rbind(result_NHBE, result_THP1, result_A549)

result$BaitProtein <- result$Bait.x

result$BaitProtein_Prey <- paste(result$BaitProtein, result$Prey.x, sep = "_")

result_filter_all <- result[which(result$BaitProtein_Prey %in% sig_prey), ]

result_filter_all$PreyGeneName <- result_filter_all$Gene.names...primary..
#result_filter_all$Strain <- unlist(sapply(strsplit(as.character(result_filter_all$Bait.x), "_"), `[`, 1))
#result_filter_all$BaitProtein <- unlist(sapply(strsplit(as.character(result_filter_all$Bait.x), "_"), `[`, 2))

mat <- reshape2::dcast(data=result_filter_all, BaitProtein + PreyGeneName ~ Strain, fun.aggregate = mean, value.var = "MIST", fill = 0)
mat$A549[which(is.na(mat$A549))] <- 0
mat$NHBE[which(is.na(mat$NHBE))] <- 0
mat$THP1[which(is.na(mat$THP1))] <- 0
mat$interaction <- "pp"
mat$type <- "prey"

# corum <- read.table(file = "/Users/yzhou/Downloads/corum_and_inverse.txt", header = T, sep = "\t", stringsAsFactors = F)
# 
# unisplit <- function(v) {
#   trimws(unlist(strsplit(v, split = "_")))
# }
# unilist <- lapply(corum$PPI.code, unisplit)
# 
# search_corum <- function(x) {
#   out <- list()
#   out2 <- c()
#   complex_df <- NULL
#   for(i in 1:length(unilist)) {
#     v <- sort(intersect(x, unilist[[i]])); if (length(v) == 0) next
#     
#     for (k1 in 1:length(v)) {
#       for (k2 in 1:length(v)) {
#         if (k2>k1) {
#           n <- paste(v[k1], v[k2], sep=".")
#           out[[n]] <- c(v[k1], v[k2], corum$complex_id[i], corum$complex_name[i], paste(v, collapse=" "))
#         }
#       }
#     }
#   }
#   return(out)
# }
# corum_list <- search_corum(unique(mat$Preys))
# corumdf <- t(as.data.frame(corum_list))
# colnames(corumdf) <- c("source","target", "corum_id", "corum_name", "protein_list")
# corumdf <- data.frame(corumdf, interaction = "corum")
# rownames(corumdf) <- NULL

# node
node_bait <- data.frame(id = unique(mat$BaitProtein), name = unique(mat$BaitProtein), type = "bait", A549 = NA, NHBE = NA,  THP1 = NA, pchart = NA)

node_prey <- unique(mat[,c("PreyGeneName","PreyGeneName", "type", "A549", "NHBE", "THP1")])
node_prey$pchart <- pchart
colnames(node_prey) <- colnames(node_bait)
node_df <- rbind(node_bait, node_prey)

# edge
edge_df <- unique(mat[,c("BaitProtein","PreyGeneName", "interaction")])
edge_df$corum_id <- edge_df$corum_name <- edge_df$protein_list <- NA
colnames(edge_df) <- c("source","target", "corum_id", "corum_name", "protein_list")
# colnames(edge_df) <- colnames(corumdf)
# edge_df <- rbind(edge_df, corumdf)

cytolist <- list()
cytolist[["NHBE_THP1_A549"]][["node"]] <- node_df
cytolist[["NHBE_THP1_A549"]][["edge"]] <- edge_df

## start to creat Cytoscape
library(biomaRt)
library(RCy3)
library(igraph)
for (cell in names(cytolist)) {
  cytolist[[cell]]$edge$target <- as.vector(cytolist[[cell]]$edge$target)
  cytolist[[cell]]$edge$source <- as.vector(cytolist[[cell]]$edge$source)
  # cytolist[[cell]]$node$id <- as.character(cytolist[[cell]]$node$id)
  # cytolist[[cell]]$node$name <- as.character(cytolist[[cell]]$node$name)
  # cytolist[[cell]]$node$type <- as.character(cytolist[[cell]]$node$type)
  nnet<-createNetworkFromDataFrames(edges = cytolist[[cell]]$edge, title=cell, collection="flu-PPIs")
  loadTableData(cytolist[[cell]]$node, "id", network = nnet)
  
  flu.style <- "flu-ppi-style"
  copyVisualStyle("BioPAX",flu.style)
  setNodeColorMapping("type",c("bait","prey"),c("#CCCCCC","#99CCFF"),"d","#CCCCCC",style.name = flu.style)
  setNodeShapeMapping("type",c("bait","prey"),c("ROUND_RECTANGLE","ELLIPSE"), "ELLIPSE", style.name = flu.style)
  setNodeSizeMapping("type",c("bait","prey"),c(30,40),"d", style.name = flu.style)
  #setNodeBorderWidthMapping("known",c("unknown","gold","silver","func"),c(1,8,8,8),"d",1,style.name = flu.style)
  #setNodeBorderColorMapping("known",c("unknown","gold","silver","func"),c("#888888","#feb24c","#6a51a3","#ef3b2c"),"d","#888888",style.name = flu.style)
  setEdgeColorMapping("interaction",c("pp","corum"),c("#969696","#41b6c4"),"d","#525252",style.name = flu.style)
  setEdgeLineStyleMapping("interaction",c("pp","corum"),line.styles=c("SOLID","DOT"),default.line.style="SOLID",style.name = flu.style)
  setEdgeLineWidthMapping("interaction",c("pp","corum"),widths=c(3,3),"d",default.width =3,style.name = flu.style)
  
  setEdgeLineWidthDefault(2, style.name = flu.style)
  pchart.map<-mapVisualProperty('node customgraphics 1','pchart','p')
  updateStyleMapping(flu.style, pchart.map)
  
  setVisualStyle(flu.style, network = nnet)
  layoutNetwork('force-directed defaultSpringCoefficient=.00006 defaultSpringLength=80', network = nnet)
}
#exportImage('/Users/yzhou/Downloads/horseshoe_bat/Bat_Nature_MIST>0.75','PDF')

corumdf_Genename <- merge(corumdf, node_df[, c("id", "name")], by.x = "source", by.y = "id")

corumdf_Genename <- corumdf_Genename[, c("source")]

pie(c(0.33, 0.33, 0.34),
    labels = c("NHBE", "THP1", "A549"), 
    col = c("#4292c6", "#807dba", "#41ab5d"), 
    init.angle = 90, cex = 5)
