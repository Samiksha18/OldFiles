### RNAseq Expression Analysis 
#     Nick Waters
#       20150127
#        Version 0.1

### 

### http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
#http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
###  Check this out too for getting access to the finder-hiden directories:
#http://mrox.net/blog/2008/08/09/learning-the-terminal-on-the-mac-part-4-bringing-finder-and-terminal-together/
###Source:  http://bioinformatics.knowledgeblog.org/2011/06/21/volcano-plots-of-microarray-data/



library(ggplot2)
library(DESeq2)
library(Biobase)
library(plyr)
df<-read.csv('/Volumes/NICK3/R/data/RNAseq_sandbox.csv', header=T,sep=",", row.names=1)


#      Set Range of samples
countsTable <-data.matrix(df[0:2301,1:10])
samples<-names(df[1:10])
condition<-c(rep("ctrl",2),rep("A",2),rep("B",2),rep("C",2),rep("D",2))
pData = data.frame(cbind(samples, condition))
ddsfm <- DESeqDataSetFromMatrix(countData = countsTable, colData=pData, design=~condition)

dds<-DESeq(ddsfm)

res<-results(dds, contrast=c("condition", "B", "ctrl"))


######  Make a Proper Apply or Loop Function to analyze all

conds<-unique(condition)
aRes<-results(dds, contrast=c("condition", "A", "ctrl"))
bRes<-results(dds, contrast=c("condition", "B", "ctrl"))
cRes<-results(dds, contrast=c("condition", "C", "ctrl"))
dRes<-results(dds, contrast=c("condition", "D", "ctrl"))

###Source:  http://bioinformatics.knowledgeblog.org/2011/06/21/volcano-plots-of-microarray-data/
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
no_of_genes = length(rownames(dds))

cutFC<-2
cutPV<-0.05
res$threshold = as.factor(abs(res$log2FoldChange) > cutFC & res$pvalue < cutPV/no_of_genes)

###Plot limits
xlm<-c(-5, 5)
ylm<-c(0, 15)

volcanoTime<-function(inpdds){
  conds<-unique(condition)
  g<-ggplot(data=data.frame(res), aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold)) +
    geom_point(alpha=0.4, size=4) +
    theme(legend.position = "none") +
    labs(x=paste("log2 fold change"), y= paste("-log10 p-value"))+
    xlim(xlm) +
    ylim(ylm)
  g+geom_text(data=data.frame(res), aes(label=df$Gene[0:2301], x=log2FoldChange, y=-log10(pvalue)), size=3, colour="black", vjust=-1)+
    geom_hline(yintercept=cutPV, color="green")+
    geom_vline(xintercept=c(-cutFC, cutFC), color="red")
}

volcanoTime(aRes)

set<-volcanoTime(dds)

##  VOLCANO!!!!!
g = ggplot(data=data.frame(res), aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold)) +
  geom_point(alpha=0.4, size=4) +
  theme(legend.position = "none") +
  labs(x=paste("log2 fold change"), y= paste("-log10 p-value"))+
  xlim(xlm) +
  ylim(ylm) 
g

### above plot+ lines for cutoffs
g+geom_text(data=data.frame(res), aes(label=df$Gene[0:2301], x=log2FoldChange, y=-log10(pvalue)), size=3, colour="black", vjust=-1)+
  geom_hline(yintercept=cutPV, color="green")+
  geom_vline(xintercept=c(-cutFC, cutFC), color="red")


plotMA(dds,ylim=c(-2,2),main="DESeq2")

