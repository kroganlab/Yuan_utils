#edgeR
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library(limma)
#mRNA
mRNA_nonnormalized=read.csv(file = "/Users/zhouyuan/Downloads/Snow/data/RNA_without_normalize.csv",header=T)
coldata=read.csv(file = "/Users/zhouyuan/Downloads/Snow/data/coldata.csv",header=T)
rownames(mRNA_nonnormalized)=mRNA_nonnormalized$Gene.Symbol
mRNA=mRNA_nonnormalized[,6:83]

#Only consider nHCC vs Adj.Norm
ind=which(coldata$labels=="nHCC"|coldata$labels=="Adj.Norm")
coldata_sub=coldata[ind,]
mRNA_sub=mRNA[,ind]

dge2<-DGEList(counts=mRNA_sub,group=coldata_sub$labels)
dge2<-calcNormFactors(dge2)
dge2<-estimateCommonDisp(dge2)
de.com2<-exactTest(dge2)
summary(decideTestsDGE(de.com2,adjust.method="BH"))
sum.de2<-summary(decideTestsDGE(de.com2,adjust.method="BH"))
sig.p2<-topTags(de.com2,n=sum.de2[1]+sum.de2[3],adjust.method="BH",sort.by="p.value")
sig.counts<-mRNA[row.names(sig.p2),]
sig.pre2<-dge2$pseudo.counts[row.names(sig.p2),]
#Take the Gene Annotation data of the selected significant genes
#sig.annot2<-rnaseq.count.annt[row.names(sig.p2),]
#Combine the count data and the analysis result of selected genes
sig.genes.complete.edgeR2<-cbind(sig.p2$table,sig.pre2)
#sig.genes.complete.edgeR2=sig.genes.complete.edgeR2[,-c(3:7)]
write.csv(sig.genes.complete.edgeR2,"sig.genes.complete.edgeR2.csv")
#write.csv(de.com2$table,"de.com2.csv")