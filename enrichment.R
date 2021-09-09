entichment_gene_list_dir <- function(data, label, direction)
{
  df <- read.csv(file = data, header = T, stringsAsFactor = F)
  df_sub <- df[which(df$Label %in% label), ]
  df_sub <- df_sub[-which(is.na(df_sub$log2FC)), ]
  
  #universe gene list
  gene_universe <- unique(df_sub$Protein)
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene_universe))
  if(length(ind) > 0)
  {
    gene_universe_ind <- as.character(gene_universe[ind])
    gene_universe_ind <- unlist(strsplit(gene_universe_ind, split = ";"))
    gene_universe <- gene_universe[-ind]
    gene_universe <- c(gene_universe, gene_universe_ind)
    # remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    # remove site
    a <- unlist(strsplit(gene_universe, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene_universe <- unique(b[, 1])
  }else{
    # remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    # remove site
    a <- unlist(strsplit(gene_universe, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene_universe <- unique(b[, 1])
  }
  
  #gene <- df_sub$Protein[which(abs(df_sub$log2FC)> 1 & df_sub$adj.pvalue < 0.05)]
  if(direction == "up")
  {
    gene <- df_sub$Protein[which((is.finite(df_sub$log2FC) & df_sub$log2FC > 1 & df_sub$adj.pvalue < 0.05) |
                                   (df_sub$log2FC == Inf & df_sub$totalRepPos >= 4))]
    message("There are ", sum(df_sub$log2FC == Inf & df_sub$totalRepPos >= 4), " genes have Inf")
  }else{
    gene <- df_sub$Protein[which((is.finite(df_sub$log2FC) & df_sub$log2FC < (-1) & df_sub$adj.pvalue < 0.05) |
                                   (df_sub$log2FC == -Inf & df_sub$totalRepNeg >= 8))]
    message("There are ", sum(df_sub$log2FC == -Inf & df_sub$totalRepNeg >= 8), " genes have -Inf")
  }
  
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene))
  if(length(ind) > 0)
  {
    gene_ind <- as.character(gene[ind])
    gene_ind <- unlist(strsplit(gene_ind, split = ";"))
    gene <- gene[-ind]
    gene <- c(gene, gene_ind)
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    # remove site
    a <- unlist(strsplit(gene, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene <- unique(b[, 1])
  }else{
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    # remove site
    a <- unlist(strsplit(gene, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene <- unique(b[, 1])
  }
  
  return(list(gene = gene, gene_universe = gene_universe))
}

library(clusterProfiler)
library(RColorBrewer)
library(tidyr)
# enrichment
enrichment <- function(gene, universe)
{
  out_tab <- NULL
  # GO enrichment
  ego <- enrichGO (gene = gene, 
                   #universe = universe, 
                   keyType = "UNIPROT", 
                   OrgDb = "org.Mm.eg.db", 
                   #ont = "ALL",
                   ont = "BP",  
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1, 
                   pool = F)
  ego_df <- as.data.frame(ego, stringsAsFactor = FALSE)
  
  out_tab <- rbind(out_tab, ego_df)
  
  # #KEGG enrichment
  # ke <- enrichKEGG (gene = gene, 
  #                   universe = universe, 
  #                   keyType = "uniprot", 
  #                   organism = "mmu",
  #                   pAdjustMethod = "BH", 
  #                   pvalueCutoff = 1,
  #                   qvalueCutoff = 1)
  # ke_df <- as.data.frame(ke, stringsAsFactor = FALSE)
  # ke_df <- data.frame(ONTOLOGY = "KEGG", ke_df, stringsAsFactors = F)
  # 
  # out_tab <- rbind(out_tab, ke_df)
  
  gene_ratio <- t(as.data.frame( lapply(out_tab$GeneRatio, function(x){as.numeric(trimws(unlist(strsplit(x, split = "/"))))}) ))
  bg_ratio <- t(as.data.frame( lapply(out_tab$BgRatio, function(x){as.numeric(trimws(unlist(strsplit(x, split = "/"))))}) ))
  odds_ratio <- (gene_ratio[,1]/gene_ratio[,2])/((bg_ratio[,1]-gene_ratio[,1])/(bg_ratio[,2]-gene_ratio[,2])); names(odds_ratio) <- NULL
  out_tab$odds_ratio <- odds_ratio
  
  # ## change ID in out_tab
  # out_tab$GOID <- out_tab$ID
  # out_tab$ID <- out_tab$Description
  # out_tab$ID <- gsub(" ","_" , paste("GO_", toupper(out_tab$ID), sep = ""))
  
  return(out_tab)
}

library(data.table)
simplifyEnrichBySimilarUniverseMembership <- function(enrichResultsTable, gmt, groupColumn="bait"){
  if (length(unique(enrichResultsTable$ID)) < 2){
    message ("Nothing to simplify")
    return (list (enrichResultsTable, data.frame()))
  }
  setDT(enrichResultsTable); setDT(gmt)
  ##Select Significant Terms
  target_overrep_sig <- enrichResultsTable[qvalue < 0.01]
  #target_overrep_sig <- enrichResultsTable[pvalue < 0.01]
  ##Prepare Significant GO Term Jaccard Similarity Matrix
  sig_go_terms <- unique(target_overrep_sig$ID)
  message ("Computing universal gene overlap between ", length(sig_go_terms), " significant GO terms from ", length(unique(enrichResultsTable[[groupColumn]])), " ", groupColumn, "(s)")
  gmt.subset <- gmt[ont %in% sig_go_terms, .(ont=factor(ont), gene=factor(gene))]
  termByGeneMat <-  Matrix::sparseMatrix(as.integer(gmt.subset$ont), as.integer(gmt.subset$gene), 
                                         dimnames=list(levels(gmt.subset$ont), 
                                                       levels(gmt.subset$gene)))
  go_dist_mat <- dist(termByGeneMat, method="binary")
  
  
  # ############################################
  # ##Prepare Significant GO Term Jaccard Similarity Matrix
  # #sig_go_terms <- unique(target_overrep_sig$ID)
  # sig_go_terms_int <- sig_go_terms[sig_go_terms %in% gmt$ont]
  # go_sim_mat <- matrix(NA,nrow=length(sig_go_terms_int),ncol=length(sig_go_terms_int),
  #                      dimnames=list(sig_go_terms_int,sig_go_terms_int))
  # for( i in 1:length(sig_go_terms_int))
  # {
  #   gs1 <- gmt$gene[gmt$ont == sig_go_terms_int[i]]
  #   for( j in 1:length(sig_go_terms_int))
  #   {
  #     gs2 <- gmt$gene[gmt$ont == sig_go_terms_int[j]]
  #     go_sim_mat[i,j] <- length(intersect(gs1,gs2))/length(union(gs1,gs2))
  # 
  #   }
  # }
  # ##Construct tree based on Jaccard Similarity and cut tree at specific level
  # ##We used h=0.99 as it seemed to give reasonable number of biologically meaningful clusters
  # go_sim_cluster <- hclust(1-as.dist(go_sim_mat))
  # go_sim_tree <- cutree(go_sim_cluster,h=0.99)
  # ##Select the broadest term in each cluster
  # sig_bait_go_pairs <-
  #   do.call('rbind', lapply(unique(as.character(
  #     target_overrep_sig$label
  #   )), function(bait) {
  #     gs <-
  #       as.character(target_overrep_sig[target_overrep_sig$label == bait, ]$ID)
  #     nr_gs <-
  #       as.character(tapply(gs, go_sim_tree[gs], function(i) {
  #         names(sort(c5_gene_sets_length[i], decreasing = TRUE))[1]
  #       }))
  #     df <-
  #       target_overrep_sig[target_overrep_sig$Bait == bait &
  #                            target_overrep_sig$ID %in% nr_gs, ]
  #   }))
  # ############################################
  
  hc <- hclust(go_dist_mat)
  clusters <- cutree(hc, h=0.99)
  clusters <- data.table (cluster = as.numeric(clusters), ID = attributes(clusters)$names )
  message ("GO terms clustered into ", max(clusters$cluster), " clusters")
  ## go gene set lengths
  gmt.setLengths <- gmt[,.(setSize = length(unique(gene))), by = ont]
  clusters <- merge (clusters, gmt.setLengths, by.x="ID", by.y="ont")
  ## local gene set lengths  -- this will only count those genes that make it into an enriched list
  id.gene.long <- enrichResultsTable[,.(ID, gene = unlist(strsplit(geneID, "/"))), by = seq_len(nrow(enrichResultsTable))]
  genesPerTerm <- id.gene.long[,.(localSetLength=length(unique(gene))), by = ID]
  clusters <- merge (clusters, genesPerTerm, by.x="ID", by.y="ID")
  #setorder(clusters, -localSetLength)  # the tie breaker
  clusterInfo <- merge (target_overrep_sig, clusters, by = "ID")
  clusterInfo[,maxSet := max (setSize), by = c("cluster", groupColumn)]
  #winners <- clusterInfo[,.SD[which(count == maxSet)], by = .(cluster, Bait)]  #keeps all ties...needs updating
  winners <- clusterInfo[,.SD[which.max(setSize),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
  message (length(unique(winners$ID)), " representative GO terms choosing the broadest significant term per GO-cluster per ", groupColumn)
  result <- enrichResultsTable[ID %in% winners$ID,]
  list(simplified = result, clusterInfo = clusterInfo)
}

library(RColorBrewer)
ms_enrichGOPlot <- function(c, tit = "H1N1 Enrichment", outfile = "H1N1_enrichment")
{
  terms <- unique(c$Description)
  c$odds_ratio <- log2(c$odds_ratio); c$odds_ratio[c$odds_ratio>6] <- 6
  comps <- unique(c$compare)
  comps_num <- 1:length(comps); names(comps_num) <- comps
  byterm <- ifelse(length(terms)>10, 1, 0.5)
  terms_num <- rev(seq(1,length(terms),byterm)[1:length(terms)]); names(terms_num) <- terms
  c$comp <- comps_num[c$compare]
  c$term <- terms_num[c$Description]
  c$p2 <- round(100*c$pvalue)
  col_vec <- rev(colorRampPalette(brewer.pal(9,"Blues"))(101))
  c$color <- col_vec[c$p2+1]
  
  c$cex <- c$odds_ratio/1.2
  #tit <- "PH"
  tcex <- ifelse(max(nchar(terms)) <= 80, 1, 0.7)
  ccex <- ifelse(max(nchar(comps)) <= 21, 1, 0.7)
  
  width <- 8 + round(length(unique(c$compare))/8)
  pdf(paste(outfile, ".pdf", sep = ""), width = width, height = 10)
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), widths = c(9,1.2))
  par(mar=c(10,30,3,1))
  plot(c$comp, c$term, pch = 20, cex = 0.001, col = "white", axes=F, xlab="", ylab="", main = "")
  for (n in comps_num) abline(v=n, col = "#f0f0f0", lwd=0.5)
  for (n in terms_num) abline(h=n, col = "#f0f0f0", lwd=0.5)
  points(c$comp, c$term, pch = 19, cex = c$cex, col = c$color, xpd = TRUE)
  text(1,max(terms_num), offset=2, labels=tit,pos=3,xpd = TRUE, cex=1.2)
  text(x=comps_num+0.1, par("usr")[3], labels =  names(comps_num), srt = 45, pos = 2, xpd = TRUE, cex = ccex)
  text(par("usr")[1], terms_num, labels = names(terms_num), pos = 2, xpd = TRUE, cex = tcex)
  
  par(mar=c(2,1,5,0.5))
  or_min <- ceiling(min(c$odds_ratio)); or_max <- ceiling(max(c$odds_ratio))
  or_num <- or_min:or_max
  or_cex <- or_num/1.2
  plot(x=rep(1,length(or_num)), y = 1:length(or_num), cex = or_cex, axes=F, pch=19, xlab="", ylab="", xpd = TRUE)
  axis(3, 1, "Odds Ratio\n(log2)", tick=F)
  text(1.1, 1:length(or_num), labels = or_num, pos =4, xpd = TRUE)
  
  par(mar=c(2,0.5,4,0.5))
  padj <- seq(min(c$p2), max(c$p2), 5)
  padj_col <- col_vec[padj+1]
  plot(x=rep(1,length(padj)), y = 1:length(padj), cex = 14, axes=F, pch='.', xlab="", ylab="", col = padj_col)
  axis(3, 1, "P value", tick=F)
  text(1.05, 1:length(padj), labels = padj/100, pos=4, xpd = TRUE)
  dev.off()
}

fixMsigdbGONames <- function(names){
  names <- gsub("^GO_","",names)
  names <- gsub("_"," ",names)
  names <- tolower(names)
  patterns <- c("\\bdna|dna\\b", 
                "\\brna|rna\\b", 
                "\\bgtp|gtp\\b",
                "\\batp|atp\\b"
  )
  replacements <- c("DNA", 
                    "RNA", 
                    "GTP",
                    "ATP"
  ) 
  for (i in seq_along(patterns)){
    names <- gsub (patterns[i], replacements[i], names)
  }
  #capitalize first letter
  substr(names,1,1) <- toupper(substr(names,1,1))
  return(names)
}

library (ComplexHeatmap)
enrichHeatmapBestPerGroup <- function(simplifiedEnrichTable, fullEnrichTable, groupColumn="bait", topN = 1, title="", cols = NULL, 
                                      negCols = NULL, reduceRedundantsAcrossGroups=TRUE,...){
  setorder(simplifiedEnrichTable, p.adjust)
  # bestTermPerBait <- simplifiedEnrichTable[p.adjust<0.01,.(ID=ID[1:topN]),by=groupColumn]
  
  bestTermPerBait <- simplifiedEnrichTable[,.(ID=ID[1:topN]),by=groupColumn]
  #bestTermPerBait <- simplifiedEnrichTable
  if(reduceRedundantsAcrossGroups){  
    #reduce redundant based on clusters in fullEnrichTable
    countsByID <- fullEnrichTable[ID %in% bestTermPerBait$ID, .(geneCount  = length(unique(unlist(strsplit(geneID, "/"))))), by = .(ID, cluster)]
    # get the term with most genes across whole dataset per term-cluster
    setorder(countsByID, -geneCount)
    bestTerms <- countsByID[,.SD[1],by=cluster]$ID
  } else bestTerms <- unique(bestTermPerBait$ID)
  
  main.wide <- dcast (fullEnrichTable[ID %in% bestTerms], as.formula(paste("Description", groupColumn, sep="~")), value.var="p.adjust")
  for(col in unique(c(cols,negCols))){
    if (is.null(main.wide[[col]])) main.wide[[col]] <- NA
  }
  
  main.mat <- -log10(as.matrix(main.wide, rownames = "Description"))
  main.mat[is.na(main.mat)] <- 0
  main.mat[main.mat > 5] <- 5
  rownames(main.mat) <- fixMsigdbGONames(rownames(main.mat))
  
  counts.wide <- dcast (fullEnrichTable[ID %in% bestTerms], as.formula(paste("Description", groupColumn, sep="~")), value.var="Count")
  for(col in unique(c(cols,negCols))){
    if (is.null(counts.wide[[col]])) counts.wide[[col]] <- NA
  }
  counts.mat <- as.matrix(counts.wide, rownames="Description")
  
  
  geneTable <- fullEnrichTable[ID %in% bestTerms, .(gene = unlist(strsplit(geneID, split="/"))),by = ID]
  geneTable[,cleanName := fixMsigdbGONames(ID)]
  
  if (!is.null(cols)){
    if (! all(cols %in% colnames(counts.mat) & cols %in% colnames(main.mat))){
      message ("Not all requested columns for heatmap found in data")
      message ("main.mat ", paste(colnames(main.mat), collapse=" "))
      message ("counts.mat ", paste(colnames(counts.mat), collapse=" "))
    }else{
      counts.mat<- counts.mat[,cols]
      main.mat<- main.mat[,cols]
    }
    
  }
  #print(str(counts.mat))
  
  # temporary for comparison with splitCircle
  #ddr <- as.dendrogram(hclust(dist(main.mat[,c(posCols, negCols)])))
  
  
  Blues = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  #Reds = colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
  colors <- Blues(100)
  
  heatmap_legend_param = list(legend_direction="horizontal", title = "-log10(adj.p)")
  
  if (!is.null(negCols)){
    colors <- circlize::colorRamp2 (breaks=seq(from=-max(main.mat), to = max(main.mat), length.out=101), colors =colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))
    main.mat[,negCols] = -main.mat[, negCols]
    heatmap_legend_param = c (heatmap_legend_param, list(at=c(-4,-2,0,2,4), labels = c(4,2,0,2,4)) )
  }
  
  
  ##Plot main figure heatmap
  hm <- ComplexHeatmap::Heatmap(main.mat, col = colors, border = TRUE, rect_gp = gpar(col = "grey", lwd = 1),
                                #cluster_rows = ddr,
                                column_title = title,
                                column_names_rot = 90, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                                show_row_dend = FALSE, show_column_dend = FALSE, heatmap_legend_param = heatmap_legend_param,
                                row_names_max_width = max_text_width(rownames(main.mat), gp = gpar(fontsize = 12)),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if (!is.na(counts.mat[i,j])){
                                    grid.text(sprintf("%.0f", counts.mat[i, j]), x, y, gp = gpar(fontsize=10, col="white"))
                                  }
                                }, ...)
  
  draw(hm,heatmap_legend_side="top")
  invisible(geneTable)
}


# # Simplify
# ## Select Significant Terms
# target_overrep_sig <- out_tab %>% filter(pvalue < 0.01)
# ## Prepare Significant GO Term Jaccard Similarity Matrix
# sig_go_terms <- unique(target_overrep_sig$ID)
# go_sim_mat <- matrix(NA,nrow=length(sig_go_terms),ncol=length(sig_go_terms),
#                      dimnames=list(sig_go_terms,sig_go_terms))
# 
# for( i in 1:length(sig_go_terms))
# {
#   gs1 <- c5_gene_sets_list[[sig_go_terms[i]]]
#   for( j in 1:length(sig_go_terms))
#   {
#     gs2 <- c5_gene_sets_list[[sig_go_terms[j]]]        
#     go_sim_mat[i,j] <- length(intersect(gs1,gs2))/length(union(gs1,gs2))
#     
#   }
# }



# Mouse_UB_H1N1
library(clusterProfiler)
Mouse_UB <- "/Users/yzhou/Downloads/flu2Yuan/output/mouse_UB_MSstats.csv"
label <- c("H1N1_D00-MOCK_D04",
           "H1N1_D01-MOCK_D04",
           "H1N1_D02-MOCK_D04",
           "H1N1_D03-MOCK_D04",
           "H1N1_D04-MOCK_D04")
#gmt <- read.gmt("/Users/yzhou/Downloads/c5.all.v7.1.symbols.gmt")
mouseGO <- clusterProfiler:::get_GO_data("org.Mm.eg.db", "BP", "UNIPROT")
gmt <- rbindlist(lapply (mouseGO$EXTID2PATHID, function(x) data.table(ont = x)), idcol="gene")
enrichResultsTable <- data.frame()
#library(openxlsx)
#wb <- createWorkbook()
for(i in 1:length(label))
{
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_UB, label[i], "up")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "up"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
  
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_UB, label[i], "down")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "down"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
}
write.table(enrichResultsTable, "Mouse_UB_H1N1.txt", sep = "\t", row.names = F, quote = F)
a <- simplifyEnrichBySimilarUniverseMembership(enrichResultsTable, gmt, groupColumn="groupColumn")
c <- a$simplified

enrichHeatmapBestPerGroup(a$simplified, a$clusterInfo, groupColumn="groupColumn", topN = 1, title="Mouse UB H1N1", 
                          cols = c("H1N1_D03-MOCK_D04_up", "H1N1_D04-MOCK_D04_up", "H1N1_D02-MOCK_D04_up", "H1N1_D00-MOCK_D04_up", "H1N1_D01-MOCK_D04_up", 
                                   "H1N1_D00-MOCK_D04_down", "H1N1_D04-MOCK_D04_down", "H1N1_D01-MOCK_D04_down", "H1N1_D02-MOCK_D04_down", "H1N1_D03-MOCK_D04_down"), 
                          negCols = c("H1N1_D00-MOCK_D04_down", "H1N1_D04-MOCK_D04_down", "H1N1_D01-MOCK_D04_down", "H1N1_D02-MOCK_D04_down", "H1N1_D03-MOCK_D04_down"), 
                          reduceRedundantsAcrossGroups=T)


# enrichHeatmapBestPerGroup(a$simplified, enrichResultsTable, groupColumn="groupColumn", topN = 1, title="Mouse UB H1N1", 
#                           cols = c("H1N1_D03-MOCK_D04_up", "H1N1_D04-MOCK_D04_up", "H1N1_D02-MOCK_D04_up", "H1N1_D00-MOCK_D04_up", "H1N1_D01-MOCK_D04_up", 
#                                    "H1N1_D00-MOCK_D04_down", "H1N1_D04-MOCK_D04_down", "H1N1_D01-MOCK_D04_down", "H1N1_D02-MOCK_D04_down", "H1N1_D03-MOCK_D04_down"), 
#                           negCols = c("H1N1_D00-MOCK_D04_down", "H1N1_D04-MOCK_D04_down", "H1N1_D01-MOCK_D04_down", "H1N1_D02-MOCK_D04_down", "H1N1_D03-MOCK_D04_down"), 
#                           reduceRedundantsAcrossGroups=F)


#, breaks = c(6, 200, 490)
c$omics <- c$label
c$tag <- c$dir
library(ggplot2)
library(tidyverse)
p1 <- c %>% 
  ggplot(aes(x = omics, y = Description, color = qvalue, size = odds_ratio, shape = tag)) +
  geom_point() +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  scale_size("Odds Ratio", range = c(2, 8)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'right', 
        legend.direction = 'vertical',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) + 
  ggtitle("H1N1 UB enrichment")

cairo_pdf("H1N1_UB.pdf", family="Arial Unicode MS", 4,8)
p1
dev.off()




c$compare <- c$groupColumn
ms_enrichGOPlot(c, tit = "H1N1 UB Enrichment", outfile = "H1N1_UB_enrichment")




library(ggplot2)
ggplot(c, aes(x = label, y = Description, color = pvalue, size = odds_ratio)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c("blue", 'white', "red")) +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  labs(color = expression('log'[2]*'(abundance'['MUT']*'/abundance'['WT']*')')) +
  scale_size("FDR", range = c(5.2, 2.4), breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1)) +
  theme_classic() +
  ylab('Gsp1 point mutant') +
  xlab('\nIntersect of Gsp1 interaction partners identified by AP-MS of Gsp1 strains with both tags') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))



  ggplot(c) +
  geom_point(aes(wt-0.027, mpg), shape="\u25D6", colour="red", size=3) +
  geom_point(aes(wt+0.027, mpg), shape="\u25D7", colour="blue", size=3) +
  theme_bw()



# Mouse_UB_H5N1
library(clusterProfiler)
Mouse_UB <- "/Users/yzhou/Downloads/flu2Yuan/output/mouse_UB_MSstats.csv"
label <- c("H5N1_D00-MOCK_D04",
           "H5N1_D01-MOCK_D04",
           "H5N1_D02-MOCK_D04",
           "H5N1_D03-MOCK_D04",
           "H5N1_D04-MOCK_D04")
gmt <- read.gmt("/Users/yzhou/Downloads/c5.all.v7.1.symbols.gmt")
enrichResultsTable <- data.frame()
#library(openxlsx)
#wb <- createWorkbook()
for(i in 1:length(label))
{
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_UB, label[i], "up")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "up"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
  
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_UB, label[i], "down")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "down"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
}
a <- simplifyEnrichBySimilarUniverseMembership(enrichResultsTable, gmt, groupColumn="groupColumn")
c <- a$simplified

pdf(file = paste("Heatmap.best1PerBait", "pdf", sep="."), width=6, height=18)
enrichHeatmapBestPerGroup(a$simplified, enrichResultsTable, groupColumn="groupColumn", topN = 1, title="Mouse UB H5N1", 
                          cols = unique(a$simplified$groupColumn), 
                          negCols = c("H5N1_D00-MOCK_D04_down", "H5N1_D04-MOCK_D04_down", "H5N1_D01-MOCK_D04_down", "H5N1_D02-MOCK_D04_down", "H5N1_D03-MOCK_D04_down"), 
                          reduceRedundantsAcrossGroups=F)
dev.off()



#, breaks = c(6, 200, 490)
c$omics <- c$label
c$tag <- c$dir
library(ggplot2)
library(tidyverse)
p1 <- c %>% 
  ggplot(aes(x = omics, y = Description, color = qvalue, size = odds_ratio, shape = tag)) +
  geom_point() +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  scale_size("Odds Ratio", range = c(2, 8)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'right', 
        legend.direction = 'vertical',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) + 
  ggtitle("H5N1 UB enrichment")

cairo_pdf("H5N1_UB.pdf", family="Arial Unicode MS", 4,9)
p1
dev.off()



c$compare <- c$groupColumn
ms_enrichGOPlot(c, tit = "H5N1 UB Enrichment", outfile = "H5N1_UB_enrichment")


#################
#PH
# Mouse_PH_H1N1
library(clusterProfiler)
Mouse_PH <- "/Users/yzhou/Downloads/flu2Yuan/output/mouse_PH_MSstats.csv"
label <- c("H1N1_D00-MOCK_D04",
           "H1N1_D01-MOCK_D04",
           "H1N1_D02-MOCK_D04",
           "H1N1_D03-MOCK_D04",
           "H1N1_D04-MOCK_D04")
gmt <- read.gmt("/Users/yzhou/Downloads/c5.all.v7.1.symbols.gmt")
enrichResultsTable <- data.frame()
#library(openxlsx)
#wb <- createWorkbook()
for(i in 1:length(label))
{
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_PH, label[i], "up")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "up"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
  
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_PH, label[i], "down")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "down"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
}
a <- simplifyEnrichBySimilarUniverseMembership(enrichResultsTable, gmt, groupColumn="groupColumn")
c <- a$simplified

pdf(file = paste("Heatmap.best1PerBait", "pdf", sep="."), width=6, height=18)
enrichHeatmapBestPerGroup(a$simplified, enrichResultsTable, groupColumn="groupColumn", topN = 1, title="Mouse PH H1N1", 
                          cols = unique(a$simplified$groupColumn), 
                          negCols = c("H1N1_D00-MOCK_D04_down", "H1N1_D04-MOCK_D04_down", "H1N1_D01-MOCK_D04_down", "H1N1_D02-MOCK_D04_down", "H1N1_D03-MOCK_D04_down"), 
                          reduceRedundantsAcrossGroups=F)
dev.off()




#, breaks = c(6, 200, 490)
c$omics <- c$label
c$tag <- c$dir
library(ggplot2)
library(tidyverse)
p1 <- c %>% 
  ggplot(aes(x = omics, y = Description, color = qvalue, size = odds_ratio, shape = tag)) +
  geom_point() +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  scale_size("Odds Ratio", range = c(2, 8)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'right', 
        legend.direction = 'vertical',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) + 
  ggtitle("H1N1 PH enrichment")

cairo_pdf("H1N1_PH.pdf", family="Arial Unicode MS", 4,5)
p1
dev.off()




c$compare <- c$groupColumn
ms_enrichGOPlot(c, tit = "H1N1 UB Enrichment", outfile = "H1N1_UB_enrichment")




library(ggplot2)
ggplot(c, aes(x = label, y = Description, color = pvalue, size = odds_ratio)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c("blue", 'white', "red")) +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  labs(color = expression('log'[2]*'(abundance'['MUT']*'/abundance'['WT']*')')) +
  scale_size("FDR", range = c(5.2, 2.4), breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1)) +
  theme_classic() +
  ylab('Gsp1 point mutant') +
  xlab('\nIntersect of Gsp1 interaction partners identified by AP-MS of Gsp1 strains with both tags') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))



ggplot(c) +
  geom_point(aes(wt-0.027, mpg), shape="\u25D6", colour="red", size=3) +
  geom_point(aes(wt+0.027, mpg), shape="\u25D7", colour="blue", size=3) +
  theme_bw()



# Mouse_PH_H5N1
library(clusterProfiler)
Mouse_PH <- "/Users/yzhou/Downloads/flu2Yuan/output/mouse_PH_MSstats.csv"
label <- c("H5N1_D00-MOCK_D04",
           "H5N1_D01-MOCK_D04",
           "H5N1_D02-MOCK_D04",
           "H5N1_D03-MOCK_D04",
           "H5N1_D04-MOCK_D04")
gmt <- read.gmt("/Users/yzhou/Downloads/c5.all.v7.1.symbols.gmt")
enrichResultsTable <- data.frame()
#library(openxlsx)
#wb <- createWorkbook()
for(i in 1:length(label))
{
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_PH, label[i], "up")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "up"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
  
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_list_dir(Mouse_PH, label[i], "down")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "down"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
}
a <- simplifyEnrichBySimilarUniverseMembership(enrichResultsTable, gmt, groupColumn="groupColumn")
c <- a$simplified


pdf(file = paste("Heatmap.best1PerBait", "pdf", sep="."), width=6, height=15)
enrichHeatmapBestPerGroup(a$simplified, enrichResultsTable, groupColumn="groupColumn", topN = 1, title="Mouse PH H5N1", 
                          cols = unique(a$simplified$groupColumn), 
                          negCols = c("H5N1_D00-MOCK_D04_down", "H5N1_D04-MOCK_D04_down", "H5N1_D01-MOCK_D04_down", "H5N1_D02-MOCK_D04_down", "H5N1_D03-MOCK_D04_down"), 
                          reduceRedundantsAcrossGroups=F)
dev.off()



#, breaks = c(6, 200, 490)
c$omics <- c$label
c$tag <- c$dir
library(ggplot2)
library(tidyverse)
p1 <- c %>% 
  ggplot(aes(x = omics, y = Description, color = qvalue, size = odds_ratio, shape = tag)) +
  geom_point() +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  scale_size("Odds Ratio", range = c(2, 8)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'right', 
        legend.direction = 'vertical',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) + 
  ggtitle("H5N1 PH enrichment")

cairo_pdf("H5N1_PH.pdf", family="Arial Unicode MS", 4,8)
p1
dev.off()



c$compare <- c$groupColumn
ms_enrichGOPlot(c, tit = "H5N1 UB Enrichment", outfile = "H5N1_UB_enrichment")






## Mouse_AB

entichment_gene_wosite_dir <- function(data, label, direction)
{
  df <- read.csv(file = data, header = T, stringsAsFactor = F)
  df_sub <- df[which(df$Label %in% label), ]
  if(length(which(is.na(df_sub$log2FC))) > 0)
  {
    df_sub <- df_sub[-which(is.na(df_sub$log2FC)), ]
  }
  
  #universe gene list
  gene_universe <- unique(df_sub$Protein)
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene_universe))
  if(length(ind) > 0)
  {
    gene_universe_ind <- as.character(gene_universe[ind])
    gene_universe_ind <- unlist(strsplit(gene_universe_ind, split = ";"))
    gene_universe <- gene_universe[-ind]
    gene_universe <- c(gene_universe, gene_universe_ind)
    #remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    gene_universe <- unique(gene_universe)
  }else{
    #remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    gene_universe <- unique(gene_universe)
  }
  
  if(direction == "up")
  {
    gene <- df_sub$Protein[which((is.finite(df_sub$log2FC) & df_sub$log2FC > 1 & df_sub$adj.pvalue < 0.05))]
    message("There are ", sum(df_sub$log2FC == Inf), " genes have Inf")
  }else{
    gene <- df_sub$Protein[which((is.finite(df_sub$log2FC) & df_sub$log2FC < (-1) & df_sub$adj.pvalue < 0.05))]
    message("There are ", sum(df_sub$log2FC == -Inf), " genes have -Inf")
  }
  
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene))
  if(length(ind) > 0)
  {
    gene_ind <- as.character(gene[ind])
    gene_ind <- unlist(strsplit(gene_ind, split = ";"))
    gene <- gene[-ind]
    gene <- c(gene, gene_ind)
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    gene <- unique(gene)
  }else{
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    gene <- unique(gene)
  }
  
  return(list(gene = gene, gene_universe = gene_universe))
}


# Mouse_AB_H1N1
library(clusterProfiler)
Mouse_AB <- "/Users/yzhou/Downloads/flu2Yuan/output/mouse_AB_H1N1_MSstatsTMT.csv"
label <- c("D00-MOCK_D04", 
           "D01-MOCK_D04", 
           "D02-MOCK_D04", 
           "D03-MOCK_D04", 
           "D04-MOCK_D04")
gmt <- read.gmt("/Users/yzhou/Downloads/c5.all.v7.1.symbols.gmt")
enrichResultsTable <- data.frame()
#library(openxlsx)
#wb <- createWorkbook()
for(i in 1:length(label))
{
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_wosite_dir(Mouse_AB, label[i], "up")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "up"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
  
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_wosite_dir(Mouse_AB, label[i], "down")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "down"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
}
a <- simplifyEnrichBySimilarUniverseMembership(enrichResultsTable, gmt, groupColumn="groupColumn")
c <- a$simplified


pdf(file = paste("Heatmap.best1PerBait", "pdf", sep="."), width=6, height=15)
enrichHeatmapBestPerGroup(a$simplified, enrichResultsTable, groupColumn="groupColumn", topN = 1, title="Mouse AB H1N1", 
                          cols = unique(a$simplified$groupColumn), 
                          negCols = c("D00-MOCK_D04_down", "D03-MOCK_D04_down"), 
                          reduceRedundantsAcrossGroups=F)
dev.off()




# Mouse_AB_H5N1
library(clusterProfiler)
Mouse_AB <- "/Users/yzhou/Downloads/flu2Yuan/output/mouse_AB_H5N1_MSstatsTMT.csv"
label <- c("D00-MOCK_D04", 
           "D01-MOCK_D04", 
           "D02-MOCK_D04", 
           "D03-MOCK_D04", 
           "D04-MOCK_D04")
gmt <- read.gmt("/Users/yzhou/Downloads/c5.all.v7.1.symbols.gmt")
enrichResultsTable <- data.frame()
#library(openxlsx)
#wb <- createWorkbook()
for(i in 1:length(label))
{
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_wosite_dir(Mouse_AB, label[i], "up")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "up"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
  
  #enri <- entichment_gene_list(Mouse_UB, label[i])
  enri <- entichment_gene_wosite_dir(Mouse_AB, label[i], "down")
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrich_out$dir <- "down"
  enrich_out$groupColumn <- paste(enrich_out$label, enrich_out$dir, sep = "_")
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
  
}
a <- simplifyEnrichBySimilarUniverseMembership(enrichResultsTable, gmt, groupColumn="groupColumn")
c <- a$simplified


pdf(file = paste("Heatmap.best1PerBait", "pdf", sep="."), width=6, height=15)
enrichHeatmapBestPerGroup(a$simplified, enrichResultsTable, groupColumn="groupColumn", topN = 1, title="Mouse AB H5N1", 
                          cols = c("D00-MOCK_D04_down", 
                                   "D01-MOCK_D04_down", 
                                   "D02-MOCK_D04_down", 
                                   "D03-MOCK_D04_down", 
                                   "D04-MOCK_D04_down", 
                                   "D00-MOCK_D04_up", 
                                   "D01-MOCK_D04_up", 
                                   "D02-MOCK_D04_up", 
                                   "D03-MOCK_D04_up", 
                                   "D04-MOCK_D04_up"), 
                          negCols = c("D00-MOCK_D04_down", 
                                      "D01-MOCK_D04_down", 
                                      "D02-MOCK_D04_down", 
                                      "D03-MOCK_D04_down", 
                                      "D04-MOCK_D04_down"), 
                          reduceRedundantsAcrossGroups=F, cluster_columns=FALSE)
dev.off()




















# Mouse_UB
library(clusterProfiler)
Mouse_UB <- "/Users/yzhou/Downloads/flu2Yuan/output/mouse_UB_MSstats.csv"
label <- c("H1N1_D00-MOCK_D04",
           "H1N1_D01-MOCK_D04",
           "H1N1_D02-MOCK_D04",
           "H1N1_D03-MOCK_D04",
           "H1N1_D04-MOCK_D04", 
           "H5N1_D00-MOCK_D04",
           "H5N1_D01-MOCK_D04",
           "H5N1_D02-MOCK_D04",
           "H5N1_D03-MOCK_D04",
           "H5N1_D04-MOCK_D04")

gmt <- read.gmt("/Users/yzhou/Downloads/c5.all.v7.0.symbols.gmt")
enrichResultsTable <- data.frame()
#library(openxlsx)
#wb <- createWorkbook()
for(i in 1:length(label))
{
  enri <- entichment_gene_list(Mouse_UB, label[i])
  enrich_out <- enrichment(enri$gene, enri$gene_universe)
  #addWorksheet(wb = wb, sheetName = label[i], gridLines = T)
  #writeDataTable(wb = wb, sheet = label[i], x = enrich_out)
  enrich_out$label <- label[i]
  enrichResultsTable <- rbind(enrichResultsTable, enrich_out)
}
a <- simplifyEnrichBySimilarUniverseMembership(enrichResultsTable, gmt, groupColumn="label")
c <- a$simplified
c$compare <- c$label
ms_enrichGOPlot(c, tit = "UB Enrichment", outfile = "UB_enrichment")








#saveWorkbook(wb,"/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx", overwrite = T)
library(readxl)
excel_sheets("/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx")
D00 <- read_excel("/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx", sheet = "H1N1_D00-MOCK_D04")
D01 <- read_excel("/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx", sheet = "H1N1_D01-MOCK_D04")
D02 <- read_excel("/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx", sheet = "H1N1_D02-MOCK_D04")
D03 <- read_excel("/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx", sheet = "H1N1_D03-MOCK_D04")
D04 <- read_excel("/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx", sheet = "H1N1_D04-MOCK_D04")

D00 <- D00[, c("ONTOLOGY", "ID", "Description", "pvalue", "odds_ratio")]
D01 <- D01[, c("ONTOLOGY", "ID", "Description", "pvalue", "odds_ratio")]
D02 <- D02[, c("ONTOLOGY", "ID", "Description", "pvalue", "odds_ratio")]
D03 <- D03[, c("ONTOLOGY", "ID", "Description", "pvalue", "odds_ratio")]
D04 <- D04[, c("ONTOLOGY", "ID", "Description", "pvalue", "odds_ratio")]

e <- merge(D00, D01, by= c("ONTOLOGY", "ID", "Description"), all = T)
e <- merge(e, D02, by= c("ONTOLOGY", "ID", "Description"), all = T)
e <- merge(e, D03, by= c("ONTOLOGY", "ID", "Description"), all = T)
e <- merge(e, D04, by= c("ONTOLOGY", "ID", "Description"), all = T)
colnames(e) <- c("ONTOLOGY", "ID", "Description", paste(rep(c("pvalue", "odds_ratio"), 5), rep("_D", 10), rep(0:4, each = 2), sep = ""))

# remove unsignificant terms across all time points
e1 <- e
#a <- as.matrix(e1[, c(4, 6, 8, 10, 12)])
e1[, c(4, 6, 8, 10, 12)] <- apply(e1[, c(4, 6, 8, 10, 12)], 2, function(x){x[which(is.na(x))] = 1
return(x)})
ind <- which(e1$pvalue_D0>0.05 & e1$pvalue_D1>0.05 & e1$pvalue_D2>0.05 & e1$pvalue_D3>0.05 & e1$pvalue_D4>0.05 )
e1 <- e[-ind, ]
addWorksheet(wb = wb, sheetName = "all", gridLines = T)
writeDataTable(wb = wb, sheet = "all", x = e1)
saveWorkbook(wb,"/Users/yzhou/Downloads/flu2Yuan/Enrichment dot plot/enrich_Mouse_UB_H1N1.xlsx", overwrite = T)


mgoSim 


s <- simplify(
  ego,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)





library (igraph)
library (GO.db)


groupGOTerms <- function(terms){
  #get all edges:
  edgesAsList <- lapply (terms, FUN=function(t) intersect(terms,GO.db::GOBPPARENTS[[t]]))
  names(edgesAsList) <- terms
  
  .child_parent_table <- function (c){
    if(length(edgesAsList[[c]])==0) parent <- "none"
    else parent <- edgesAsList[[c]]
    data.frame(child = c, parent = parent)
  }
  edgesAsTable <- rbindlist(lapply (names(edgesAsList), FUN=.child_parent_table ))
  dag <- graph_from_data_frame(edgesAsTable)
  dag <- delete.vertices(dag,"none")
}











entichment_gene_list_dir <- function(data, label, direction)
{
  df <- read.csv(file = data, header = T, stringsAsFactor = F)
  df_sub <- df[which(df$Label %in% label), ]
  df_sub <- df_sub[-which(is.na(df_sub$log2FC)), ]
  
  #universe gene list
  gene_universe <- unique(df_sub$Protein)
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene_universe))
  if(length(ind) > 0)
  {
    gene_universe_ind <- as.character(gene_universe[ind])
    gene_universe_ind <- unlist(strsplit(gene_universe_ind, split = ";"))
    gene_universe <- gene_universe[-ind]
    gene_universe <- c(gene_universe, gene_universe_ind)
    #remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    # remove site
    a <- unlist(strsplit(gene_universe, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene_universe <- unique(b[, 1])
  }else{
    #remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    # remove site
    a <- unlist(strsplit(gene_universe, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene_universe <- unique(b[, 1])
  }
  
  
  #gene list
  # remove Inf/-Inf
  df_sub <- df_sub[-which(df_sub$log2FC == Inf), ]
  df_sub <- df_sub[-which(df_sub$log2FC == -Inf), ]
  
  gene <- c()
  if(direction == "up")
  {
    for(g in unique(df_sub$Protein))
    {
      ind <- which(df_sub$Protein == g)
      a <- df_sub[ind, ]
      if(sum(a$log2FC> 1 & a$adj.pvalue < 0.05, na.rm = T) >=3)
        gene <- c(gene, g)
    }
  }else{
    for(g in unique(df_sub$Protein))
    {
      ind <- which(df_sub$Protein == g)
      a <- df_sub[ind, ]
      if(sum(a$log2FC< (-1) & a$adj.pvalue < 0.05, na.rm = T) >=3)
        gene <- c(gene, g)
    }
  }
  
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene))
  if(length(ind) > 0)
  {
    gene_ind <- as.character(gene[ind])
    gene_ind <- unlist(strsplit(gene_ind, split = ";"))
    gene <- gene[-ind]
    gene <- c(gene, gene_ind)
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    # remove site
    a <- unlist(strsplit(gene, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene <- unique(b[, 1])
  }else{
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    # remove site
    a <- unlist(strsplit(gene, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene <- unique(b[, 1])
  }
  
  return(list(gene = gene, gene_universe = gene_universe))
}

entichment_gene_wosite_dir <- function(data, label, direction)
{
  df <- read.csv(file = data, header = T, stringsAsFactor = F)
  df_sub <- df[which(df$Label %in% label), ]
  if(length(which(is.na(df_sub$log2FC))) > 0)
  {
    df_sub <- df_sub[-which(is.na(df_sub$log2FC)), ]
  }
  
  #universe gene list
  gene_universe <- unique(df_sub$Protein)
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene_universe))
  if(length(ind) > 0)
  {
    gene_universe_ind <- as.character(gene_universe[ind])
    gene_universe_ind <- unlist(strsplit(gene_universe_ind, split = ";"))
    gene_universe <- gene_universe[-ind]
    gene_universe <- c(gene_universe, gene_universe_ind)
    #remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    gene_universe <- unique(gene_universe)
  }else{
    #remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    gene_universe <- unique(gene_universe)
  }
  
  
  #gene list
  # remove Inf/-Inf
  if(length(which(df_sub$log2FC == Inf)) > 0)
  {
    df_sub <- df_sub[-which(df_sub$log2FC == Inf), ]
  }
  if(length(which(df_sub$log2FC == -Inf)) > 0)
  {
    df_sub <- df_sub[-which(df_sub$log2FC == -Inf), ]
  }
  
  gene <- c()
  if(direction == "up")
  {
    for(g in unique(df_sub$Protein))
    {
      ind <- which(df_sub$Protein == g)
      a <- df_sub[ind, ]
      if(sum(a$log2FC> 1 & a$adj.pvalue < 0.05, na.rm = T) >=3)
        gene <- c(gene, g)
    }
  }else if(direction == "down"){
    for(g in unique(df_sub$Protein))
    {
      ind <- which(df_sub$Protein == g)
      a <- df_sub[ind, ]
      if(sum(a$log2FC< (-1) & a$adj.pvalue < 0.05, na.rm = T) >=3)
        gene <- c(gene, g)
    }
  }else{
    for(g in unique(df_sub$Protein))
    {
      ind <- which(df_sub$Protein == g)
      a <- df_sub[ind, ]
      if(sum(abs(a$log2FC)> 1 & a$adj.pvalue < 0.05, na.rm = T) >=3)
        gene <- c(gene, g)
    }
  }
  
  
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene))
  if(length(ind) > 0)
  {
    gene_ind <- as.character(gene[ind])
    gene_ind <- unlist(strsplit(gene_ind, split = ";"))
    gene <- gene[-ind]
    gene <- c(gene, gene_ind)
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    gene <- unique(gene)
  }else{
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    gene <- unique(gene)
  }
  
  return(list(gene = gene, gene_universe = gene_universe))
}

library(clusterProfiler)
library(RColorBrewer)
# enrichment entrez
enrichment_entrez <- function(gene, universe)
{
  out_tab <- NULL
  # GO enrichment
  ego <- enrichGO (gene = gene, 
                   universe = universe, 
                   keyType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db", 
                   ont = "ALL",
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1, 
                   pool = F)
  ego_df <- as.data.frame(ego, stringsAsFactor = FALSE)
  
  out_tab <- rbind(out_tab, ego_df)
  
  #KEGG enrichment
  ke <- enrichKEGG (gene = gene, 
                    universe = universe, 
                    keyType = "ncbi-geneid", 
                    organism = "mmu",
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
  ke_df <- as.data.frame(ke, stringsAsFactor = FALSE)
  ke_df <- data.frame(ONTOLOGY = "KEGG", ke_df, stringsAsFactors = F)
  
  out_tab <- rbind(out_tab, ke_df)
  
  gene_ratio <- t(as.data.frame( lapply(out_tab$GeneRatio, function(x){as.numeric(trimws(unlist(strsplit(x, split = "/"))))}) ))
  bg_ratio <- t(as.data.frame( lapply(out_tab$BgRatio, function(x){as.numeric(trimws(unlist(strsplit(x, split = "/"))))}) ))
  odds_ratio <- (gene_ratio[,1]/gene_ratio[,2])/((bg_ratio[,1]-gene_ratio[,1])/(bg_ratio[,2]-gene_ratio[,2])); names(odds_ratio) <- NULL
  out_tab$odds_ratio <- odds_ratio
  
  return(out_tab)
}


entichment_gene_list <- function(data, label)
{
  df <- read.csv(file = data, header = T, stringsAsFactor = F)
  df_sub <- df[which(df$Label %in% label), ]
  df_sub <- df_sub[-which(is.na(df_sub$log2FC)), ]
  
  #universe gene list
  gene_universe <- unique(df_sub$Protein)
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene_universe))
  if(length(ind) > 0)
  {
    gene_universe_ind <- as.character(gene_universe[ind])
    gene_universe_ind <- unlist(strsplit(gene_universe_ind, split = ";"))
    gene_universe <- gene_universe[-ind]
    gene_universe <- c(gene_universe, gene_universe_ind)
    # remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    # remove site
    a <- unlist(strsplit(gene_universe, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene_universe <- unique(b[, 1])
  }else{
    # remove viral protein
    gene_universe <- gene_universe[grep(pattern = "H1N1|H3N2|H5N1", gene_universe, invert = T)]
    # remove site
    a <- unlist(strsplit(gene_universe, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene_universe <- unique(b[, 1])
  }
  
  
  #gene list
  # remove Inf/-Inf/NA
  df_sub <- df_sub[-which(df_sub$log2FC == Inf), ]
  df_sub <- df_sub[-which(df_sub$log2FC == -Inf), ]
  
  gene <- df_sub$Protein[which(abs(df_sub$log2FC)> 1 & df_sub$adj.pvalue < 0.05)]
  
  # resolve group prtein
  ind <- which(grepl(pattern = ";", gene))
  if(length(ind) > 0)
  {
    gene_ind <- as.character(gene[ind])
    gene_ind <- unlist(strsplit(gene_ind, split = ";"))
    gene <- gene[-ind]
    gene <- c(gene, gene_ind)
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    # remove site
    a <- unlist(strsplit(gene, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene <- unique(b[, 1])
  }else{
    # remove viral protein
    gene <- gene[grep(pattern = "H1N1|H3N2|H5N1", gene, invert = T)]
    # remove site
    a <- unlist(strsplit(gene, split = "_"))
    b <- matrix(a, ncol = 2, byrow = T)
    gene <- unique(b[, 1])
  }
  
  return(list(gene = gene, gene_universe = gene_universe))
}

















