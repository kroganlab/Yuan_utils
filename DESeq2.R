if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
filter_gene_counts=read.csv(file = "/Users/zhouyuan/Downloads/Snow/data/RNA_without_normalize.csv",header=T)
coldata=read.csv(file = "/Users/zhouyuan/Downloads/Snow/data/coldata.csv",header=T)
# Gene count should be integer
num_col<-ncol(filter_gene_counts)
for(i in 6:num_col){
  filter_gene_counts[,i]=as.integer(filter_gene_counts[,i])}
dim(filter_gene_counts)
filter_gene_counts=filter_gene_counts[,5:num_col]
rownames(filter_gene_counts)=filter_gene_counts[,1]
filter_gene_counts=filter_gene_counts[,-1]
colname=coldata[,1]
rownames(coldata)=colname
coldata=coldata[,2:11]

#nHCC AA v.s. Adj.Norm AA
subcoldata=coldata[coldata$Race=="AA",]
subcoldata1=subcoldata[subcoldata$labels=="nHCC"|subcoldata$labels=="Adj.Norm",]
row=rownames(subcoldata1)
subgene=filter_gene_counts[,colnames(filter_gene_counts)%in%row]
# DESeq2 need two datasets : count matrix (rownames : gene name; colnames : sample name) & table of sample information (rownames : sample name)
# The design indicates how to model the samples, here, that we want to measure the effect of the condition, controlling for batch differences. 
# The two factor variables batch and condition should be columns of coldata.
dds <- DESeqDataSetFromMatrix(countData = subgene,
                              colData = subcoldata1,
                              design= ~ PatientID+labels)
dds <- DESeq(dds)
res <- results(dds)
res
rld=rlog(dds)
plotPCA(rld,intgroup=c("labels"))
write.csv(as.data.frame(res),file="nHCC AA v.s. Adj.Norm AA.csv")
table(subcoldata1$labels)













