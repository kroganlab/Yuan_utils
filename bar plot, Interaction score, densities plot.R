result <- read.table(file = "/Users/yzhou/Downloads/A3_PPI_Interactome/result_SAINTspc_CompPASS.txt", header = T, sep = "\t", stringsAsFactors = F)
# High-Confidence Interactions : SAINT BFDR <0.05 and a wd (per bait) > 0.9 & remove prey which appear in GFPs
result_filter <- result[which(result$BFDR < 0.05 & result$wd_percentile_perBait > 0.9), ]
remove_prey <- result_filter$Prey.x[which(result_filter$Bait.x == "GFP" | result_filter$Bait.x == "GFP_RNAse")]
result_filter <- result_filter[-which(result_filter$Bait.x == "GFP" | result_filter$Bait.x == "GFP_RNAse"), ]
ind_remove_prey <- which(result_filter$Prey.x %in% remove_prey)
result_filter <- result_filter[-ind_remove_prey, ]

# bar plot
result_filter_bar <- result_filter
result_filter_bar$RNAse <- grepl("_RNAse", result_filter_bar$Bait.x)
result_filter_bar$RNAse[which(result_filter_bar$RNAse == TRUE)] <- "RNAse +"
result_filter_bar$RNAse[which(result_filter_bar$RNAse == FALSE)] <- "RNAse -"
result_filter_bar$Bait <- result_filter_bar$Bait.x
result_filter_bar$Bait <- gsub("_RNAse", "", result_filter_bar$Bait)
count_df <- table(result_filter_bar$RNAse, result_filter_bar$Bait)
ord <- names(sort(apply(count_df, 2, sum), decreasing=T))
count_df <- count_df[,ord]
library(RColorBrewer)
strain_cols <- brewer.pal(n=3, name="Set1")
names(strain_cols) <- c("RNAse +","RNAse -")
barplot(count_df, beside = T, border = NA, col = strain_cols[rownames(count_df)], ylab = "Number of interactions",
        main = paste0("APOBEC3", ", total=",sum(count_df)), cex.names = 1,
        legend = rownames(count_df) )



# calculate Interaction Score (K) : Avg(WD, 1-BFDR)
result1 <- result
result1 <- result1[-which(is.na(result1$Bait.x)), ]
result1$RNAse <- grepl("_RNAse", result1$Bait.x)
result1$RNAse[which(result1$RNAse == TRUE)] <- "RNAse +"
result1$RNAse[which(result1$RNAse == FALSE)] <- "RNAse -"
result1$Bait <- gsub("_RNAse", "", result1$Bait.x)
result1$Bait_Prey1 <- paste(result1$Bait, result1$Prey.x, sep = "-")

result_wide_BFDR <- reshape2::dcast(result1, Bait_Prey1 + Bait + PreyGene ~ RNAse , value.var="BFDR")
colnames(result_wide_BFDR)[4:5] <- paste("BFDR", colnames(result_wide_BFDR)[4:5], sep = "_")

result_wide_wd_percentile_perBait <- reshape2::dcast(result1, Bait_Prey1 + Bait + PreyGene ~ RNAse , value.var="wd_percentile_perBait")
colnames(result_wide_wd_percentile_perBait)[4:5] <- paste("wd_percentile_perBait", colnames(result_wide_wd_percentile_perBait)[4:5], sep = "_")

result_wide_AvgSpec <- reshape2::dcast(result1, Bait_Prey1 + Bait + PreyGene ~ RNAse , value.var="AvgSpec")
colnames(result_wide_AvgSpec)[4:5] <- paste("AvgSpec", colnames(result_wide_AvgSpec)[4:5], sep = "_")

result_wide_AvePSM <- reshape2::dcast(result1, Bait_Prey1 + Bait + PreyGene ~ RNAse , value.var="AvePSM")
colnames(result_wide_AvePSM)[4:5] <- paste("AvePSM", colnames(result_wide_AvePSM)[4:5], sep = "_")

result_wide <- merge(result_wide_BFDR, result_wide_wd_percentile_perBait, by = c("Bait_Prey1", "Bait", "PreyGene"))
result_wide <- merge(result_wide, result_wide_AvgSpec, by = c("Bait_Prey1", "Bait", "PreyGene"))
result_wide <- merge(result_wide, result_wide_AvePSM, by = c("Bait_Prey1", "Bait", "PreyGene"))


# 1-BFDR
result_wide$`BFDR_RNAse +` <- 1 - result_wide$`BFDR_RNAse +`
result_wide$`BFDR_RNAse -` <- 1 - result_wide$`BFDR_RNAse -`

result1_filter <- result1[which(result1$BFDR < 0.05 & result1$wd_percentile_perBait > 0.9), ]
remove_prey <- result1_filter$Prey.x[which(result1_filter$Bait.x == "GFP" | result1_filter$Bait.x == "GFP_RNAse")]
# remove GFPs
result1_filter <- result1_filter[-which(result1_filter$Bait == "GFP"), ]
ind_remove_prey <- which(result1_filter$Prey.x %in% remove_prey)
result1_filter <- result1_filter[-ind_remove_prey, ]


result1_wide_filter <- result_wide[which(result_wide$Bait_Prey1 %in% result1_filter$Bait_Prey1), ]
result1_wide_filter$`BFDR_RNAse -`[is.na(result1_wide_filter$`BFDR_RNAse -`)] <- 0
result1_wide_filter$`BFDR_RNAse +`[is.na(result1_wide_filter$`BFDR_RNAse +`)] <- 0
result1_wide_filter$`wd_percentile_perBait_RNAse -`[is.na(result1_wide_filter$`wd_percentile_perBait_RNAse -`)] <- 0
result1_wide_filter$`wd_percentile_perBait_RNAse +`[is.na(result1_wide_filter$`wd_percentile_perBait_RNAse +`)] <- 0

result1_wide_filter$RNAse <- apply(result1_wide_filter[, c("BFDR_RNAse +", "wd_percentile_perBait_RNAse +")], 1, mean)
result1_wide_filter$NT <- apply(result1_wide_filter[, c("BFDR_RNAse -", "wd_percentile_perBait_RNAse -")], 1, mean)
result1_wide_filter$DIS <- result1_wide_filter$RNAse - result1_wide_filter$NT

# remove "_HUMAN" in PreyGene
result1_wide_filter$PreyGene <- gsub("_HUMAN", "", result1_wide_filter$PreyGene)
write.table(result1_wide_filter, "result_filter.txt", row.names = F, sep = "\t", quote = F)


# plot densities
library(sm)
result1_wide_filter$Bait_num <- as.factor(result1_wide_filter$Bait)
result1_wide_filter$Bait_num <- as.numeric(result1_wide_filter$Bait_num)
sm.density.compare(result1_wide_filter$DIS, result1_wide_filter$Bait_num, xlab="Differential Interaction Score (DIS)")

# add legend via mouse click
Bait.f <- factor(result1_wide_filter$Bait, levels= unique(result1_wide_filter$Bait_num), labels = unique(result1_wide_filter$Bait))
colfill<-c(2:(2+length(levels(Bait.f))))
legend(locator(1), levels(Bait.f), fill=colfill)

#box plot
boxplot(result1_wide_filter$DIS ~ result1_wide_filter$Bait, xlab = "Baits", ylab = "DIS")
