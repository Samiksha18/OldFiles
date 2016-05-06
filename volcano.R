### Volcano Plot for Gene Expression
#     Nick Waters
#       20150121
#        Version 1

### http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
#http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
###  Check this out too for getting access to the finder-hiden directories:
#http://mrox.net/blog/2008/08/09/learning-the-terminal-on-the-mac-part-4-bringing-finder-and-terminal-together/



###Source:  http://bioinformatics.knowledgeblog.org/2011/06/21/volcano-plots-of-microarray-data/
library(ggplot2)
library(DESeq2)

library(Biobase)

df<-read.csv('/Volumes/NICK3/R/data/RNAseq_sandbox.csv', header=T,sep=",", row.names=1)


#####  Define Treatments
con1<-"WT_a"
con2<-'WT_a'
null1<-'null_a'
null2<-'null_b'
testA1<-'R61K_a'  
testA2<-'R61K_b'  
testB1<-'R61H_a' 
testB2<- 'R61H_b'  
testC1<-'F71Y_a' 
testC2<- 'F71Y_b'

treatments<-list(c(con1, con2, null1, null2, testA1, testA2, testB1, testB2, testC1,testC2))


####     Manual calculations of log2foldchanges
res$log2fc_null1<-log2(res[[null1]]/res[[con1]])
res$log2fc_null2<-log2(res[[null2]]/res[[con2]])
res$log2fc_testA1<-log2(res[[testA1]]/res[[con1]])
res$log2fc_testB1<-log2(res[[testB1]]/res[[con1]])
res$log2fc_testC1<-log2(res[[testC1]]/res[[con1]])
res$log2fc_testA2<-log2(res[[testA2]]/res[[con2]])
res$log2fc_testB2<-log2(res[[testB2]]/res[[con2]])
res$log2fc_testC2<-log2(res[[testC2]]/res[[con2]])


###      SET Range of samples
countsTable <-data.matrix(df[0:2301,1:10])
samples<-names(df[1:10])
condition<-c(rep("ctrl",2),rep("A",2),rep("B",2),rep("C",2),rep("D",2))
pData = data.frame(cbind(samples, condition))
ddsfm <- DESeqDataSetFromMatrix(countData = countsTable, colData=pData, design=~condition)

dds<-DESeq(ddsfm)

res<-results(dds)





###Source:  http://bioinformatics.knowledgeblog.org/2011/06/21/volcano-plots-of-microarray-data/
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
no_of_genes = length(rownames(dds))

cutFC<-2
cutPV<-0.05
res$threshold = as.factor(abs(res$log2FoldChange) > cutFC & res$pvalue < cutPV/no_of_genes)

###Plot limits
xlm<-c(-5, 5)
ylm<-c(0, 15)

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


