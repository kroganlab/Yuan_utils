individual.plot<-function(x,Group,Batch,cc=NULL,pp=NULL,
                          Group.label=NULL,main.text=NULL,y.lab="log Intensity",
                          type="mean",y.range=NULL,
                          axis.size=1,label.size=1,main.size=1,point.size=1)
{
  Group<-as.factor(Group)
  Batch<-as.factor(Batch)
  GL<-levels(Group)
  G<-length(GL)
  if(is.null(cc)==TRUE){cc<-1:G}
  if(is.null(pp)==TRUE){pp<-1:G}
  if(is.null(Group.label)==TRUE)
  {
    Group.label<-GL
  }
  if(is.null(y.range)==TRUE)
  {
    y.range<-range(x,na.rm=TRUE)
  }
  BL<-levels(Batch)
  B<-length(BL)
  x1=x2<-list()
  plot(c(0.5,B+0.5),range(x,na.rm=TRUE),type="n",xlab="",ylab=y.lab,ylim=y.range,
       main=main.text,
       cex.axis=axis.size,cex.lab=label.size,cex.main=main.size,xaxt="n")
  #title(main=title)
  for(i in 1:B)
  {
    
    z<-list()
    for(j in 1:G)
    {
      z[[j]]<-x[which(Group==GL[j]&Batch==BL[i])]
      if(type=="median")
      {
        zm<-median(z[[j]],na.rm=TRUE)
      }else
      {
        zm<-mean(z[[j]],na.rm=TRUE)
      }
      axis(1,i-0.5+1/(G+1)*j,Group.label[j],cex.axis=axis.size)
      u<-runif(length(z[[j]]),min=i-0.5+1/(G+1)*j-1/(G*6+1),max=i-0.5+1/(G+1)*j+1/(G*6+1))
      points(u,z[[j]],type="p",pch=pp[j],col=cc[j],cex=point.size)#data points
      s<-seq(i-0.5+1/(G+1)*j-1/(G*2.5+1),i-0.5+1/(G+1)*j+1/(G*2.5+1),0.001)
      #points(s,rep(zm,length(s)),pch=16,col=cc[j])
      lines(s,rep(zm,length(s)))#median or mean line
      #text(i-0.5+3/(2*(G+1)),max(x)-(max(x)-min(x))*0.05,BL[i])
    }}
  #lines(rep(1+0.5,2),y.range,type="l",lty=2)#vertical line
  
  # text(i-0.5+3/(2*(G+1))-0.1,max(x)-(max(x)-min(x))*0.1,"FC=")
  # text(i-0.5+3/(2*(G+1))+0.1,max(x)-(max(x)-min(x))*0.1,fc[i])
  # text(i-0.5+3/(2*(G+1))-0.1,max(x)-(max(x)-min(x))*0.15,"p=")
  # text(i-0.5+3/(2*(G+1))+0.1,max(x)-(max(x)-min(x))*0.15,p[i])
}

par(mfrow=c(1,2))
a=mRNA[,colnames(mRNA)=="MARCO"]
b=miRNA[,colnames(miRNA)=="hsa-miR-10b-5p"]
individual.plot(c(a,b),rep(y,2),c(rep("MARCO",98),rep("hsa-miR-10b-5p",98)),type="median",cc=c("red","blue"))
plot(b,a,col=rep(c("red","blue"),49),pch=rep(c(1,2),49),xlab="hsa-miR-10b-5p",ylab="MARCO")
abline(lm(a ~ b),col="black")
legend(12,12, legend=c("Tumor", "Normal"),col=c("red", "blue"),pch=c(1,2),cex=1,bty="n")
par(mfrow=c(1,1))
paste("MARCO","hsa-miR-10b-5p",sep = "_")


individual.plot(RNA_sub[,colnames(RNA_sub)=="LINC01296"],NULL,type="median")


pairs_selected_mRNAs=read.csv(file.choose(),header=T)

for(i in 1:52)
{
  a=mRNA[,colnames(mRNA)==pairs_selected_mRNAs$mRNA[i]]
  b=miRNA[,colnames(miRNA)==pairs_selected_mRNAs$miRNA[i]]
  mypath <- file.path("/Users","zhouyuan","Desktop","data","TCGA","mRNA_miRNA_fdr_fc","Model2a","log2",paste(as.character(pairs_selected_mRNAs$mRNA[i]),"_",as.character(pairs_selected_mRNAs$miRNA[i]) , ".jpg", sep = ""))
  jpeg(file=mypath,width = 8, height = 8, units = 'in', res = 300)
  individual.plot(c(a,b),rep(y,2),c(rep(as.character(pairs_selected_mRNAs$mRNA[i]),98),rep(as.character(pairs_selected_mRNAs$miRNA[i]),98)),type="median",cc=c("red","blue"),y.lab="log Intensity")
  dev.off()
  
  a=2^(mRNA[,colnames(mRNA)==pairs_selected_mRNAs$mRNA[i]])
  b=2^(miRNA[,colnames(miRNA)==pairs_selected_mRNAs$miRNA[i]])
  mypath <- file.path("/Users","zhouyuan","Desktop","data","TCGA","mRNA_miRNA_fdr_fc","Model2a","without log2",paste(as.character(pairs_selected_mRNAs$mRNA[i]),"_",as.character(pairs_selected_mRNAs$miRNA[i]) , ".jpg", sep = ""))
  jpeg(file=mypath,width = 8, height = 8, units = 'in', res = 300)
  individual.plot(c(a,b),rep(y,2),c(rep(as.character(pairs_selected_mRNAs$mRNA[i]),98),rep(as.character(pairs_selected_mRNAs$miRNA[i]),98)),type="median",cc=c("red","blue"),y.lab="Intensity")
  dev.off()
}