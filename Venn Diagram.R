genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)

library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

# Area Proportional
library(eulerr)
pdf("/Users/yzhou/Downloads/Jo/20240522_05mgADPrSYT/output/vsSVP0/Venn Diagram.pdf")
x <- list(
  # ADPr_DB = ADPr_DB, 
  ADPr_PARP = ADPr_PARP, 
  vsNeg = unique(result_vsNeg$Protein[which(result_vsNeg$Label == "Af02-Neg02")]), 
  vsSVP0 = unique(result_vsSVP0$Protein[which(result_vsSVP0$Label == "Af02-Af0")])
)
plot(euler(x, shape = "ellipse"), quantities = list(type = c('counts',"percent")), main = "SVP 0.2")
x <- list(
  # ADPr_DB = ADPr_DB, 
  ADPr_PARP = ADPr_PARP, 
  vsNeg = unique(result_vsNeg$Protein[which(result_vsNeg$Label == "Af05-Neg05")]), 
  vsSVP0 = unique(result_vsSVP0$Protein[which(result_vsSVP0$Label == "Af05-Af0")])
)
plot(euler(x, shape = "ellipse"), quantities = TRUE, main = "SVP 0.5")
x <- list(
  # ADPr_DB = ADPr_DB, 
  ADPr_PARP = ADPr_PARP, 
  vsNeg = unique(result_vsNeg$Protein[which(result_vsNeg$Label == "Af1-Neg1")]), 
  vsSVP0 = unique(result_vsSVP0$Protein[which(result_vsSVP0$Label == "Af1-Af0")])
)
plot(euler(x, shape = "ellipse"), quantities = TRUE, main = "SVP 1")
dev.off()