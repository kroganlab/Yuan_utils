corum <- read.table(file = "/Users/yzhou/Downloads/CorumPairs.csv", header = T, sep = ",", stringsAsFactors = F)
simplifyCORUM <- function(ppi_result, bait_name = "Bait", prey_name = "Prey", corum, cutHeight = 0.99, minProteinPerComplex = 3, hclustMethod = "complete"){
  ppi_result$Bait <- ppi_result[, bait_name]
  ppi_result$Prey <- ppi_result[, prey_name]
  # get all baits' name
  Bait <- unique(ppi_result$Bait)
  # get corum complex per bait
  corumdf1 <- data.frame()
  for(b in Bait){
    a <- ppi_result[which(ppi_result$Bait == b), ]
    c <- unique(c(a$Bait, a$Prey))
    d <- corum[which(corum$gene1 %in% c), ]
    d <- d[which(d$gene2 %in% c), ]
    corumdf1 <- rbind(corumdf1, d)
  }
  corumdf1 <- unique(corumdf1)
  #Prepare Significant GO Term Jaccard Similarity Matrix
  sig_complex_terms <- unique(corumdf1$complex_name)
  
  corum_subset <- corum[corum$complex_name %in% sig_complex_terms, ]
  corum_subset1 <- unique(corum_subset[, c("gene1", "complex_name")])
  colnames(corum_subset1) <- c("gene", "complex_name")
  corum_subset2 <- unique(corum_subset[, c("gene2", "complex_name")])
  colnames(corum_subset2) <- c("gene", "complex_name")
  corum_subset <- unique(rbind(corum_subset1, corum_subset2))
  corum_subset$gene <- as.factor(corum_subset$gene)
  corum_subset$complex_name <- as.factor(corum_subset$complex_name)
  complexByGeneMat <-  Matrix::sparseMatrix(as.integer(corum_subset$complex_name), as.integer(corum_subset$gene), 
                                            dimnames=list(levels(corum_subset$complex_name), 
                                                          levels(corum_subset$gene)))
  
  complex_dist_mat <- dist(complexByGeneMat, method="binary")
  hc <- hclust(complex_dist_mat, method = hclustMethod)
  clusters <- cutree(hc, h=cutHeight)
  clusters <- data.table (cluster = as.numeric(clusters), ID = attributes(clusters)$names )
  
  # corum gene set lengths
  corum.setLengths <- as.data.frame(table(corum_subset$complex_name))
  colnames(corum.setLengths)[which(colnames(corum.setLengths) == "Freq")] <- "setSize"
  clusters <- merge (clusters, corum.setLengths, by.x="ID", by.y="Var1")
  
  # local gene set lengths -- this will only count those genes that make it into an corum list
  corumdf1_1 <- unique(corumdf1[, c("gene1", "complex_name")])
  colnames(corumdf1_1) <- c("gene", "complex_name")
  corumdf1_2 <- unique(corumdf1[, c("gene2", "complex_name")])
  colnames(corumdf1_2) <- c("gene", "complex_name")
  corumdf <- unique(rbind(corumdf1_1, corumdf1_2))
  genesPerTerm <- as.data.frame(table(corumdf$complex_name))
  colnames(genesPerTerm)[which(colnames(genesPerTerm) == "Freq")] <- "localSetLength"
  clusters <- merge (clusters, genesPerTerm, by.x="ID", by.y="Var1")
  
  # minProteinPerComplex
  clusters <- clusters[which(clusters$localSetLength >= minProteinPerComplex), ]
  
  # rank by localSetLength and setSize
  clusters1 <- clusters[, .SD[localSetLength == max(localSetLength)], by = cluster]
  clusters1 <- clusters1[, .SD[which.max(setSize)], by = cluster]
  winner <- unique(clusters1$ID)
  
  corumdf2 <- corumdf1[which(corumdf1$complex_name %in% winner), ]
  corumdf2 <- unique(corumdf2[, c("gene1", "gene2", "complex_name")])
  # still have duplications, collapse complex names
  corumdf2$gene1_gene2 <- paste(corumdf2$gene1, corumdf2$gene2, sep = "-")
  library(dplyr)
  corumdf2 <- corumdf2 %>%
    group_by(gene1_gene2) %>%
    summarise(
      gene1 = first(gene1),
      gene2 = first(gene2),
      complex_name = paste(unique(complex_name), collapse = ";"),
      .groups = "drop"
    )
  corumdf2 <- corumdf2[, -which(colnames(corumdf2) == "gene1_gene2")]
  return(corumdf2)
  
}
