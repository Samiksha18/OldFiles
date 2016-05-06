#!/usr/bin/Rscript
####  Single File Shiny App   ####

# Workflow:
# Filter for count length
# RUV
# DESEQ2
# make df of results of edgeR analysis
# do qc plots
# make psudocounts from log2FoldChanges
# compute contrasts from pseudocounts
# plot



argv<- commandArgs()

library(ggplot2)
library(dplyr)
library(shiny)
library(reshape2)
library(gridExtra)
library(RUVSeq)
library(EDASeq)
library(DESeq2)
library(car)

print(argv)
print(paste("length of argv:", length(argv)))
conff<-(grepl("CONFIG", argv))
source(conff)
#source("appCONFIG.R")

print("usage: app.R appCONFIG.R")
angryMessage<-"where's your CONFIG file? What are you trying to pull?  
\n I'm not a mindreader..."
iff<-ifelse(conff==FALSE, stop(angryMessage), "input looks ok!")
iff

###############################################################################
#####3  fist, a bit of QC
pseudoCount<-log(countsTable + 1)
bpdf = melt(pseudoCount, value.name = "pseudocounts")
head(bpdf)
bpdf = data.frame(bpdf, Condition = gsub("^(.*)_(.*)", "\\1", bpdf$Var2))
ggplot(bpdf, aes(x = Var2, y = pseudocounts, fill = Condition)) + 
  geom_boxplot()+#geom_jitter(position = position_jitter(width = .25)) + 
  xlab("") +
  ylab(expression(log[2](count + 1)))# + scale_fill_manual(values = c("#619CFF", "#F564E3", "red", "blue", "green", "purple"))
prerle<-as.data.frame(countsTable)
#prele<-mutate(prerle, medi=rowMeans(prerle))
prele<-prerle%>%mutate(median=rowMeans(.)) # i know right?
pre2<-(prerle)
pre3<-log2(pre2+1)-log2(prele$median+1)
pre3$id<-row.names(pre3)
pre4<-melt(pre3, id.vars = "id")
pre5 = data.frame(pre4, Condition = gsub("^(.*)_(.*)", "\\1", pre4$variable))

ggplot(pre5, aes(x = variable, y = value, fill = Condition)) + 
  geom_point()+geom_jitter(position = position_jitter(width = .25)) + 
  geom_boxplot()+
  xlab("") +
  ylab("ratios") + scale_fill_manual(values = c("#619CFF", "#F564E3", "red", "blue", "green", "purple"))

###############################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")
#biocLite("EDASeq")
#library(RUVSeq)
#library(EDASeq)
filter <- apply(countsTable, 1, function(x) length(x[x>7])>=2)
filtered <- countsTable[filter,]
genes <- rownames(filtered)
print(colnames(filtered))
x<-as.factor(rep(unique(gsub("(.*)(_.*)", "\\1", colnames(countsTable))), each=techReps))
x<-relevel(x, "back.WT")
#x<-relevel(x, "WT")

#x<-as.factor(rep(c("G129D", "null", "R61E", "R61H", "R61K", "WT"), each=2))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered), conditions=x))
set
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
#plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=x)#, col=colors[x])
#plotPCA(set, col=colors[x], cex=1.2)
#set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)
############      RUVg Emperical#######
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:100]))]
########
set2 <- RUVg(set, empirical, k=1)
pData(set2)
#plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2, main="set2: RUVg")
#####################################  RUVs  ###################
differences <- matrix(data=c(1:2, 3:4, 5:6,7:8),byrow=T, nrow=4)#,  9:10, 11:12), byrow=TRUE, nrow=6)
differences <- matrix(data=c(1:4, 5:8),byrow=F, nrow=2)#,  9:10, 11:12), byrow=TRUE, nrow=6)
differences
set3 <- RUVs(set, genes, k=1, differences)
pData(set3)
plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[x], main="set3: RUVs")
plotPCA(set3, col=colors[x], cex=1.2, main="set3: RUVs")

######################################## RUVr ################################
#  Make this general
set4WT<-grep("back\\.WT", colnames(countsTable), value=T)
#set4WT<-grep("WT", colnames(countsTable), value=T)
set4Test<-grep("back\\.null", colnames(countsTable), value=T)
#set4Test<-grep("null", colnames(countsTable), value=T)
prex<-c(set4WT, set4Test)
seq <- set
# get control genes
DEsubset<-counts(seq)[,prex]
xcont<-as.factor(rep(unique(gsub("(.*)(_.*)", "\\1", prex)), each=techReps))
xcont<-relevel(xcont, "back.WT")
#xcont<-relevel(xcont, "WT")

designcont <- model.matrix(~xcont)
ycont <- DGEList(counts=DEsubset, group=xcont)
ycont <- calcNormFactors(ycont, method="upperquartile")
ycont <- estimateGLMCommonDisp(ycont, designcont)
ycont <- estimateGLMTagwiseDisp(ycont, designcont)
nonDEThresh<-.01
ex<-exactTest(ycont)$table
ex$FC<-10^ex$logFC
nonDE<-ex[round(ex$FC, 1)>=0.5,]
nonDE<-nonDE[round(nonDE$FC, 1)<=1.5,]
nonDE<-row.names(nonDE[nonDE$PValue<=.1,])
fitcont <- glmFit(ycont, designcont)
rescont<-residuals(fitcont, type="deviance")
################
x <- as.factor(rep(unique(gsub("(.*)(_.*)", "\\1", colnames(countsTable))), each=techReps))
x<-relevel(x, "back.WT")
#x<-relevel(x, "WT")

design <- model.matrix(~x)
y <- DGEList(counts=counts(seq), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
#head(res)
set4<-RUVr(seq, nonDE, k=1, residuals=res, center=T, round=TRUE, epsilon=1, tolerance=1e-8)
pData(set4<-RUVr(seq, nonDE, k=1, residuals=res, center=T, round=TRUE, epsilon=1, tolerance=1e-8))
plotRLE(set4, outline=FALSE, ylim=c(-4, 4), col=colors[x], main="set4: RUVr")
plotPCA(set4, col=colors[x], cex=1.2, main="set4: RUVr")


cds1<-as(set, "CountDataSet") # lane
cds2<-as(set2, "CountDataSet") # RUVg
cds3<-as(set3, "CountDataSet") # RUVs
cds4<-as(set4, "CountDataSet") # RUVr
#countsTable<-filtered
#counts(get(method))
###############################################################################
### Make DESeq Data Set
samples<-colnames(countsTable)
head(countsTable)
samples
co<-strsplit(x=samples, split="_", fixed = TRUE, perl = FALSE, useBytes = FALSE)
condition<-unlist(lapply(co,FUN=function(x){paste(x[1],sep="")}))
#condition<-c(rep("ctrl",2),rep("A",2),rep("B",2),rep("C",2),rep("D",2))
pData0 = data.frame(cbind(samples, condition), row.names = samples)
colnames(pData0)<-c("samples", "conditions")

#
pData<-pData(set3)
#colnames(pData)<-c("samples", "conditions", "RUV")
#pData$scaled<-pData$RUVg/max(pData$RUVg)
ddsfm <- DESeqDataSetFromMatrix(countData = filtered, colData=pData0, design=~conditions)
ddsfm1 <- DESeqDataSetFromMatrix(countData = filtered, colData=pData, design= ~ W_1+conditions)
#design(ddsfm1) <- ~ pData$RUVg + condition
dds<-DESeq(ddsfm1, betaPrior = betaPriorBool)
dds0<-DESeq(ddsfm, betaPrior = betaPriorBool)

head(counts(dds, normalized=T))
head(counts(dds0, normalized=T))

#plotMA(dds, ylim=c(2,-2), main="MA Plot of ")
#dds<-dds0


#########  Make a Proper Apply or Loop Function to analyze all, but for now this works
countsdf<-as.data.frame(get(choice))
countsdf<-as.data.frame(filtered) 
#test1Term<-"null"

c<-unique(unlist(gsub("(.*)_[abc123]", "\\1", colnames(countsdf[rangeOfCounts,]))))
cc<-c[grep(control, c, invert = T)]
cc<-relevel(as.factor(cc), test1Term)
cc<-cc[order(cc)]
cd<-cc[grep(test1Term, cc, invert = T)]
cd
cd<-as.character(cd)
cd
#cc<-c[grep("back.WT", c, invert = T)]
test1<-c[grep(test1Term, c, invert = F)]
#dd<-c[grep("back.WT", c, invert = F)]
print("non-control strain list:")
print(cc)
print("non-control, non-primary test strain list:")
print(cd)
#resList<-paste(cc, "res", sep="")

aRes<-results(dds, contrast=c("conditions", test1, control))# null
try(bRes<-results(dds, contrast=c("conditions", cd[1], control)), silent=T)
try(cRes<-results(dds, contrast=c("conditions", cd[2], control)), silent=T)
try(dRes<-results(dds, contrast=c("conditions", cd[3], control)), silent=T)
try(eRes<-results(dds, contrast=c("conditions", cd[4], control)), silent=T)
##############  Extract fold Change 
aRes$reg <- ifelse(aRes$log2FoldChange<0, "up", "down")
aRes$fc<-2^(aRes$log2FoldChange)
aRes$pseudo<-100*aRes$fc

try(bRes$reg <- ifelse(bRes$log2FoldChange<0, -1, 1), silent=T)
try(bRes$fc<-2^(bRes$log2FoldChange), silent=T)
try(bRes$pseudo<-100*bRes$fc, silent=T)

try(cRes$reg <- ifelse(cRes$log2FoldChange<0, -1, 1), silent=T)   
try(cRes$fc<-2^(cRes$log2FoldChange), silent=T)
try(cRes$pseudo<-100*cRes$fc, silent=T)

try(dRes$reg <- ifelse(dRes$log2FoldChange<0, -1, 1), silent=T)   
try(dRes$fc<-2^(dRes$log2FoldChange), silent=T)
try(dRes$pseudo<-100*dRes$fc, silent=T)

try(eRes$reg <- ifelse(eRes$log2FoldChange<0, -1, 1) , silent=T)  
try(eRes$fc<-2^(eRes$log2FoldChange), silent=T)
try(eRes$pseudo<-100*eRes$fc, silent=T) 
head(aRes)
#######

prebind<- grep(".{1}Res$", ls(), value=T) 
# prebind2<- grep(".{1}Res1$", ls(), value=T) 
#####   For plots, and other qc stuff
resdf<-data.frame()
#preStatTable<-data.frame(row.names(get(resList[1])))
for (x in prebind){
  t<-as.data.frame(get(x))
  t$comp<-paste(x)
  resdf<-rbind(resdf,t)
  resdf
}
#filelist <- mapply(cbind, filelist, "SampleID"=ID, SIMPLIFY=F)
key<-data.frame(comp=prebind, cond=c(test1Term, as.character(cd)))
key
resdf<-merge(resdf, key, all.x=T, by="comp")

# MA-plot

facet_MA<-ggplot(as.data.frame(resdf),aes(baseMean,log2FoldChange))+
  geom_point(aes(colour=padj),alpha=0.7)+
  scale_colour_gradient(low="red",high="green")+
  scale_x_log10()+
  facet_wrap(~cond)
facet_MA
####  From http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html
ntop = 500
Pvars <- genefilter::rowVars(assay(rlog(dds)))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,  length(Pvars)))]
PCA <- prcomp(t(assay(rlog(dds))[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    samples = colData(rlog(dds))$samples,
                    condition = colData(rlog(dds))$condition)

(qplot(PC1, PC2, data = dataGG, color =  condition, 
       main = "PC1 vs PC2, top variable genes", size = I(6))
 + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
        y = paste0("PC2, VarExp:", round(percentVar[2],4)))
 + scale_colour_brewer(type="qual", palette=2)
)
#install.packages("GGally")
library(GGally)
Pvars2 <- genefilter::rowVars(assay(rlog(dds2)))
select2 <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,  length(Pvars)))]
#PCA2 <- prcomp(t(assay(rlog(dds2))))#[select2, ]), scale = F)
ggpairs(with(PCA, data.frame(x[,1], x[,2], x[,3], x[,4])), title = "PCA Analysis")
#ggpairs(with(PCA2, data.frame(x[,1], x[,2], x[,3], x[,4])), title = "PCA Analysis")

#################



data <- plotPCA(rlog(dds), returnData=TRUE, intgroup=c("conditions"))
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2))+#, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

################  Histogram of pvals
facet_hist<-ggplot(resdf, aes(resdf$padj))+ geom_histogram(binwidth=.05)+facet_wrap(~cond)+scale_y_log10()
facet_hist
facet_hist_trim<-ggplot(resdf[resdf$baseMean > 1,], aes(resdf$pvalue[resdf$baseMean > 1]))+ geom_histogram(binwidth=.05)+facet_wrap(~cond)+scale_y_log10()
facet_hist_trim

c
d

# Test plots
#plotCounts(dds, gene=which.min(eRes$padj), intgroup="condition") # the smallest P-val
#plotMA(eRes, main="WT vs test1Term")  # test, if you wanna


#### What if we use the above to get list of DE genes, and then use simple contrasts?
# Note!!!  squaring removes notion of positive vs negative regulation. 
##  instead, threashold will be log2'd (ie if thresh = 8x, log2thresh=2 )
#Build in Contrasts separately
resList<-grep(".Res$", ls(), value=T)
length(resList)

preStatTable<-data.frame(row.names(get(resList[1])))
#preStatTable<-data.frame()
for (x in resList){
  t<-as.data.frame(get(x)[,c("log2FoldChange", "padj")])
  colnames(t)<-paste(x,colnames(t), sep="$")
  preStatTable<-cbind(preStatTable,t)
  preStatTable
}
preStatTable$row.names.get.resList.1...<-NULL

statTable<-preStatTable
#statTable<-data.frame(cbind(aRes$log2FoldChange, aRes$padj, bRes$log2FoldChange, bRes$padj, cRes$log2FoldChange, cRes$padj, dRes$log2FoldChange, dRes$padj, eRes$log2FoldChange, eRes$padj),  row.names=row.names(aRes))  
#colList<-c("aRes$log2FoldChange", "aRes$padj", "bRes$log2FoldChange", "bRes$padj", "cRes$log2FoldChange", "cRes$padj", "dRes$log2FoldChange", "dRes$padj", "eRes$log2FoldChange", "eRes$padj")
#colnames(statTable)<-colList
statTable$gene<-row.names(statTable)
print(head(statTable))

fcl2<-log2(fc_thresh)
fcl2tight<-log2(fc_thresh_tight)
####### Sig Genes
sigTable<- filter(statTable, aRes$padj<pval_thresh)

de_a<- filter(statTable, abs(aRes$log2FoldChange)>fcl2)
try(de_b<- filter(statTable, abs(bRes$log2FoldChange)>fcl2), silent=T)
try(de_c<- filter(statTable, abs(cRes$log2FoldChange)>fcl2), silent=T)
try(de_d<- filter(statTable, abs(dRes$log2FoldChange)>fcl2), silent=T)
try(de_e<- filter(statTable, abs(eRes$log2FoldChange)>fcl2), silent=T)
deList<-grep("^de_", ls(), value=T)
#
mergedfc<-suppressWarnings(rbind(get(deList[1]), 
                                 try(get(deList[2]), silent =T ), 
                                 try(get(deList[3]), silent =T ),
                                 try(get(deList[4]), silent =T ),
                                 try(get(deList[5]), silent =T ),
                                 try(get(deList[6]), silent =T )))

mergedfc<-mergedfc[ grep("Error", mergedfc[,2], invert = TRUE) , ]

#mergedfc2<-bind_rows(dec,ded,dee)
fclist<-unique(mergedfc$gene) ########## DE Genes
fcdf<-data.frame(as.character(fclist))
fcdf$thresh<-rep(fc_thresh, length(fclist))

###############  tight
tde_a<- filter(statTable, abs(aRes$log2FoldChange)>fcl2tight)
try(tde_b<- filter(statTable, abs(bRes$log2FoldChange)>fcl2tight), silent=T)
try(tde_c<- filter(statTable, abs(cRes$log2FoldChange)>fcl2tight), silent=T)
try(tde_d<- filter(statTable, abs(dRes$log2FoldChange)>fcl2tight), silent=T)
try(tde_e<- filter(statTable, abs(eRes$log2FoldChange)>fcl2tight), silent=T)
tightDeList<-grep("tde_", ls(), value=T)
#
mergedfcTight<-suppressWarnings(rbind(get(tightDeList[1]), 
                                      try(get(tightDeList[2]), silent =T ), 
                                      try(get(tightDeList[3]), silent =T ),
                                      try(get(tightDeList[4]), silent =T ),
                                      try(get(tightDeList[5]), silent =T ),
                                      try(get(tightDeList[6]), silent =T )))

mergedfcTight<-mergedfcTight[ grep("Error", mergedfcTight[,2], invert = TRUE) , ]


#mergedfcTight<-bind_rows(deat,debt,dect,dedt,deet)
fclistTight<-unique(mergedfcTight$gene) ########## DE Genes
fcdfTight<-data.frame(as.character(fclistTight))
fcdfTight$thresh<-rep(fc_thresh_tight, length(fclistTight))

###############
# Note: below dummy variables are different for a reason- want to avoid bad join
nextFree<-length(sigTable)+1  ##  next free column that is
sigTable[nextFree]<-"Value"    
colnames(sigTable)[nextFree]<-"sigdum"

#
colnames(fcdf)<-c("gene", "thresh")#, "dummy")
colnames(fcdfTight)<-c("gene", "thresh")#", "dummy")
#uniq<-fcdf[!fcdf$gene %in% fcdfTight$gene]
fcdf$bool<-!(unlist(fcdf$gene) %in% unlist(fcdfTight$gene))   ## start separating tight and loose for combining
fcAdd<-fcdf[fcdf$bool=="TRUE",] ####  if loose gene is not in tight, get it ready to add to the tight list but with the thresh being 2 
fcAdd$bool<-NULL #get rid of marker col
fcAll<-fcdfTight #  prep recipient df
fcAll<-rbind(fcAll, fcAdd)   #  rbind
#fcAll<-fcAll %>% group_by(gene)%>% mutate(mm=max(.$thresh))  # depreciated
fcAll$fcdum<-"VALUE"
#fcAll[1]
fcdf<-fcAll
#
fcdf$gene<-as.character(fcdf$gene)
bothmerged2<-full_join(fcdf, sigTable, by = "gene")
bothmerged3<-bothmerged2[!(bothmerged2$sigdum %in% NA),] 

difAndSig<-unlist(bothmerged3[!(bothmerged3$sigdum %in% NA),]$gene)   ###### Best of Both, optimistic
altDifAndSig<-intersect(de_a$gene, sigTable$gene)                   ###### Best of Both, pessimistic

##export CSV with Sig Genes
newdf=bothmerged2[,(c("sigdum", "fcdum"))]
row.names(newdf)<-bothmerged2$gene
colnames(newdf)<-c("sig", "de")
difsig_handle<-paste(difsig_handle, "_fc_",fc_thresh,"__", "pval_", pval_thresh,".csv", sep="")
print(paste("writing file: ", difsig_handle))
write.csv(newdf, difsig_handle)

#################

#print("dont go away, or freak out, or press control-c;  we are just calculating RPKM's")
#print("and it takes a while because SOMEONE is still using for-loops.  What a loser...")
###calc RPKM's for dataset            this takes FOREVER

merge$len<-unlist(lapply(as.character(merge$sequence), FUN = nchar))
preReadsPer<-as.data.frame(counts(set4))
#preReadsPer$row<-row.names(preReadsPer)
preReadsPer$qv_anno<-gsub("(.\\S*)(\\s.*)", "\\1",row.names(preReadsPer))
preReadsPer<-merge(preReadsPer, merge, by="qv_anno", all.x = T)
readsPer<-cbind(preReadsPer[,c("qv_anno", "usa300gene", "usa300_desc", paste(colnames(countsdf),".x", sep=""), "len")])
readsPer$row<-readsPer$qv
row.names(readsPer)<-paste(readsPer$qv_anno, readsPer$usa300gene)
#readsPer$len<-as.numeric(merge$len)
readsPer$usa300gene<-NULL
colnames(readsPer)<-gsub("(.*)(\\.x)","\\1", colnames(readsPer))
pre_rpkm<-grep("_\\d", colnames(readsPer), value = T)
rownamesReadsPer<-row.names(readsPer)

try(readsPer<-readsPer %>% mutate(R01=unlist((10^9*readsPer[pre_rpkm[1]])/(as.numeric(sum(readsPer[pre_rpkm[1]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R02=unlist((10^9*readsPer[pre_rpkm[2]])/(as.numeric(sum(readsPer[pre_rpkm[2]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R03=unlist((10^9*readsPer[pre_rpkm[3]])/(as.numeric(sum(readsPer[pre_rpkm[3]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R04=unlist((10^9*readsPer[pre_rpkm[4]])/(as.numeric(sum(readsPer[pre_rpkm[4]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R05=unlist((10^9*readsPer[pre_rpkm[5]])/(as.numeric(sum(readsPer[pre_rpkm[5]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R06=unlist((10^9*readsPer[pre_rpkm[6]])/(as.numeric(sum(readsPer[pre_rpkm[6]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R07=unlist((10^9*readsPer[pre_rpkm[7]])/(as.numeric(sum(readsPer[pre_rpkm[7]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R08=unlist((10^9*readsPer[pre_rpkm[8]])/(as.numeric(sum(readsPer[pre_rpkm[8]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R09=unlist((10^9*readsPer[pre_rpkm[9]])/(as.numeric(sum(readsPer[pre_rpkm[9]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R10=unlist((10^9*readsPer[pre_rpkm[10]])/(as.numeric(sum(readsPer[pre_rpkm[10]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R11=unlist((10^9*readsPer[pre_rpkm[11]])/(as.numeric(sum(readsPer[pre_rpkm[11]]))*readsPer["len"]))))
try(readsPer<-readsPer %>% mutate(R12=unlist((10^9*readsPer[pre_rpkm[12]])/(as.numeric(sum(readsPer[pre_rpkm[12]]))*readsPer["len"]))))
rang<-grep("^R\\d\\d$", colnames(readsPer))
colnames(readsPer)[rang]<-paste(pre_rpkm,"_RPKM", sep="")
row.names(readsPer)<-rownamesReadsPer

# calc_rpkm<-function(x){
#   col<-colnames(x)
#   col<-col[col!="len"]
#   for (co in col){
#     for (row in row.names(x))
# #      browser()
#       x[row, paste(co, "_RPKM", sep="")]<-(10^9*x[row, co])/(sum(unlist(x[co]))*(unlist(x[row,"len"])))
#     }
#   return(as.data.frame(x))  
# }
# rp<-calc_rpkm(readsPer)
# rpkmosdf_old<-rbind(rp[grep("RPKM", colnames(rp))])
rpkmosdf<-readsPer[,grep("RPKM", colnames(readsPer))]
#is this even needed? Double check here if error below
rpkmosdf<-rpkmosdf[,order(colnames(rpkmosdf))]
#row.names(rpkmosdf)<-row.names(readsPer)

print(c) #double check
i<-data.frame(sapply(c, function(y) grep(y, colnames(rpkmosdf)))) #Regex those colnames
print(i)
### Define strains  ##  IS THIS TRULY DEPRECIATED?????
# control<-c[6]
# test1Term<-c[2]
# st1<-c[5]
# st2<-c[4]
# st3<-c[3]
# st4<-c[1]
reallyAngryMessage<-"What did you do?!?!?!  Your cols look odd!"

stopifnot(identical(c,d))

print("if any of the below assertions are untrue, go back and reconfigure")
print(paste("Wild Type is:",  control))   ####   Time 
print(paste("null is:",  test1Term))        ####      to
print(paste("R61K is:",  st1))             ####       check
print(paste("R61H is:",  st2))             ####          your
print(paste("R61E is:",  st3))             ####           Sanity:
print(paste("G129D is:",  st4))            #### Did all these pairs match?

###Determine Up and Down Regulation

rpkmosdf$reg <- ifelse(rowMeans(rpkmosdf[unlist(i[test1Term])])>rowMeans(rpkmosdf[unlist(i[control])]), "Down", "Up")   
head(rpkmosdf)
upreg<-subset(rpkmosdf, rpkmosdf$reg == as.character("Up"))
downreg<-subset(rpkmosdf, rpkmosdf$reg == as.character("Down"))
rpkmosdf$reg<-NULL

#########################################################################################
####  Develop this into a function....
#sort up/down
#apply contrast
#merge back together
# countsdf[paste(c[1], "_contrast", sep="")]<-0
# countsdf[paste( test1Term, "_contrast", sep="")]<-100
# countsdf[paste(st3, "_contrast", sep="")]<-((rowMeans(countsdf[unlist(i[st3])])- rowMeans(countsdf[unlist(i[c[1]])]))/
#                                                 (rowMeans(countsdf[unlist(i[ test1Term])])- rowMeans(countsdf[unlist(i[c[1]])])))*100
# countsdf[paste(st2, "_contrast", sep="")]<-((rowMeans(countsdf[unlist(i[st2])])- rowMeans(countsdf[unlist(i[c[1]])]))/
#                                                 (rowMeans(countsdf[unlist(i[ test1Term])])- rowMeans(countsdf[unlist(i[c[1]])])))*100
# countsdf[paste(st2, "_contrast", sep="")]<-((rowMeans(countsdf[unlist(i[c[5]])])- rowMeans(countsdf[unlist(i[c[1]])]))/
#                                                 (rowMeans(countsdf[unlist(i[ test1Term])])- rowMeans(countsdf[unlist(i[c[1]])])))*100
# countsdf[paste(st2, "_contrast", sep="")]<-((rowMeans(countsdf[unlist(i[c[5]])])- rowMeans(countsdf[unlist(i[c[1]])]))/
#                                                (rowMeans(countsdf[unlist(i[ test1Term])])- rowMeans(countsdf[unlist(i[c[1]])])))*100
# #Build this
# contrasts<-function(x){}
#########################################################################################
### down
downreg$reg<-NULL
upreg$reg<-NULL
downreg[paste(control, "_contrast", sep="")]<-((rowMeans(downreg[unlist(i[control])])- rowMeans(downreg[unlist(i[control])]))/
                                                  (rowMeans(downreg[unlist(i[ test1Term])])- rowMeans(downreg[unlist(i[control])])))*100
downreg[paste( test1Term, "_contrast", sep="")]<-((rowMeans(downreg[unlist(i[test1Term])])- rowMeans(downreg[unlist(i[control])]))/
                                                   (rowMeans(downreg[unlist(i[ test1Term])])- rowMeans(downreg[unlist(i[control])])))*100
downreg[paste(st1, "_contrast", sep="")]<-((rowMeans(downreg[unlist(i[st1])])- rowMeans(downreg[unlist(i[control])]))/
                                             (rowMeans(downreg[unlist(i[ test1Term])])- rowMeans(downreg[unlist(i[control])])))*100
downreg[paste(st2, "_contrast", sep="")]<-((rowMeans(downreg[unlist(i[st2])])- rowMeans(downreg[unlist(i[control])]))/
                                             (rowMeans(downreg[unlist(i[ test1Term])])- rowMeans(downreg[unlist(i[control])])))*100
downreg[paste(st3, "_contrast", sep="")]<-((rowMeans(downreg[unlist(i[st3])])- rowMeans(downreg[unlist(i[control])]))/
                                             (rowMeans(downreg[unlist(i[ test1Term])])- rowMeans(downreg[unlist(i[control])])))*100
downreg[paste(st4, "_contrast", sep="")]<-((rowMeans(downreg[unlist(i[st4])])- rowMeans(downreg[unlist(i[control])]))/
                                             (rowMeans(downreg[unlist(i[ test1Term])])- rowMeans(downreg[unlist(i[control])])))*100

## up
upreg[paste(control, "_contrast", sep="")]<-((rowMeans(upreg[unlist(i[control])])- rowMeans(upreg[unlist(i[test1Term])]))/
                                                (rowMeans(upreg[unlist(i[control])])- rowMeans(upreg[unlist(i[test1Term])])))*100
upreg[paste( test1Term, "_contrast", sep="")]<-((rowMeans(upreg[unlist(i[test1Term])])- rowMeans(upreg[unlist(i[test1Term])]))/
                                                 (rowMeans(upreg[unlist(i[control])])- rowMeans(upreg[unlist(i[test1Term])])))*100
upreg[paste(st1, "_contrast", sep="")]<-((rowMeans(upreg[unlist(i[st1])])- rowMeans(upreg[unlist(i[test1Term])]))/
                                           (rowMeans(upreg[unlist(i[control])])- rowMeans(upreg[unlist(i[test1Term])])))*100
upreg[paste(st2, "_contrast", sep="")]<-((rowMeans(upreg[unlist(i[st2])])- rowMeans(upreg[unlist(i[test1Term])]))/
                                           (rowMeans(upreg[unlist(i[control])])- rowMeans(upreg[unlist(i[test1Term])])))*100
upreg[paste(st3, "_contrast", sep="")]<-((rowMeans(upreg[unlist(i[st3])])- rowMeans(upreg[unlist(i[test1Term])]))/
                                           (rowMeans(upreg[unlist(i[control])])- rowMeans(upreg[unlist(i[test1Term])])))*100
upreg[paste(st4, "_contrast", sep="")]<-((rowMeans(upreg[unlist(i[st4])])- rowMeans(upreg[unlist(i[test1Term])]))/
                                           (rowMeans(upreg[unlist(i[control])])- rowMeans(upreg[unlist(i[test1Term])])))*100

contrastdf<-rbind(downreg[grep("_contrast", colnames(downreg))], upreg[grep("_contrast", colnames(upreg))])
#contrastdf$gene<-NULL

#####  Make a nice sig and dif csv for clustering
# Use this if no anno is needed
sigcontrast<-contrastdf[intersect(row.names(contrastdf), altDifAndSig),]
#write.csv(sigcontrast, sigcontrast_handle)


#Else:
df$rows<-row.names(df)
adf<-data.frame(altDifAndSig)
adf <-merge(adf, bothmerged3, by.x = "altDifAndSig", by.y = "gene")
yy<-df[c("rows", "usa300_uni", "usa300_desc","vir")]
alt_desc<-merge(adf, yy, by.x = "altDifAndSig", by.y = "rows", all = F)

sigcontrast$rows<-row.names(sigcontrast)
sigcontr_anno<-merge(sigcontrast, alt_desc, by.x = "rows", by.y = "altDifAndSig", all = F)
row.names(sigcontr_anno)<-sigcontr_anno$rows
print(paste("writing file: ", sigcontrast_handle))
write.csv(sigcontr_anno, sigcontrast_handle)

# #Else:
# df$rows<-row.names(df)
# adf<-data.frame(altDifAndSig)
# adf <-merge(adf, bothmerged3, by.x = "altDifAndSig", by.y = "gene")
# adf$qv_anno<-gsub("(.*?)(\\s.*)","\\1",adf$altDifAndSig)
# yy<-df[c("rows", "usa300_uni", "vir")]
# alt_desc<-merge(adf, preReadsPer, by.x = "qv_anno", by.y = "qv_anno", all = F)
# 
# sigcontrast$rows<-row.names(sigcontrast)
# sigcontr_anno<-merge(sigcontrast, alt_desc, by.x = "rows", by.y = "altDifAndSig", all = F)
# row.names(sigcontr_anno)<-sigcontr_anno$rows
# print(paste("writing file: ", sigcontrast_handle))
# write.csv(sigcontr_anno, sigcontrast_handle)

# make a nice csv for goGetting

lisst<-gsub("(^.+\\d)(\\s.+)", "\\1", altDifAndSig)
#lisst<-gsub("(.*)(not)","\\1", lisst)
head(lisst)
row.names(merge)<-merge$qv_anno
merge$len<-NULL
withseq<-merge[intersect(merge[["qv_anno"]], lisst), ]
str(withseq)
goOut<-cbind(withseq[c("qv_anno","aaseq")])
colnames(goOut)<-c("qv_anno","aaseq")
print(paste("writing file: ", goOut_handle))
write.csv(goOut, goOut_handle)

#   Whoa hold on now...  That just aint right!
# ordering<-paste(c("back.WT_1","back.WT_2", "back.null_1","back.null_2", "no.WT_1","no.WT_2", "no.null_1", "no.null_2"), "_RPKM", sep='')
 contrastdf<-rpkmosdf
cols2<-colnames(contrastdf)
ordering<-c("back.null_1_RPKM", "back.null_2_RPKM", "back.WT_1_RPKM", "back.WT_2_RPKM",
            "no.null_1_RPKM" , "no.null_2_RPKM", "no.WT_1_RPKM", "no.WT_2_RPKM") 
cols2

#################################################################################### 
#                             Time to Shine
####################################################################################
# Define a server for the Shiny app
server<-function(input, output) {
  ig<-reactive(input$gene)
  output$slicePlot <- renderPlot({
    slice<-unlist(contrastdf[input$gene,])
    g<-ggplot(data.frame(slice))+
      geom_point(aes(x=cols2, y=slice), size=5)+
      scale_x_discrete(limits=ordering)+
      labs(title = paste(title_of_experiment, "Expression Overlay"), 
           y="% expression",
           x="CodY variant")
    
    g+stat_smooth(aes(x=cols2, y=slice, group=1))
  })
  output$volcanoPlot<-renderPlot({
    ###Source:  http://bioinformatics.knowledgeblog.org/2011/06/21/volcano-plots-of-microarray-data/
    ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    no_of_genes = length(rownames(dds))
    xlm<-c(-5, 5)
    ylm<-c(0, 10)
    aRes$threshold = as.factor((abs(aRes$log2FoldChange) > fcl2) & (aRes$padj < pval_thresh/no_of_genes))
    df<-data.frame(aRes)
    df<-na.omit(df)
    g<-ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      geom_point(alpha=0.4, size=4) +
      theme(legend.position = "none") +
      labs(x=paste("log2 fold change"), y= paste("-log10 BH adj. p-value"))
    #+xlim(xlm) +ylim(ylm)
    
    g+geom_point(data=df[input$gene,], aes(x=log2FoldChange, y=-log10(padj)), size=5.5, colour="black", vjust=0)+
      geom_hline(yintercept=pval_thresh, color="green")+
      geom_vline(xintercept=c(-fcl2, fcl2), color="red")+
      scale_y_log10()
  })
  output$geneInfo <- renderText({
    paste(df[input$gene, "usa300_desc"])
  })
  output$sigLev <- renderText({
    paste(sigcontr_anno[input$gene, "thresh"])
  })
}

dataset<-"test 1"

ui<-fluidPage(    
  titlePanel("Brinsmade Lab Gene Expression Browser"),
  sidebarLayout(          
    # Define the sidebar with one input
    sidebarPanel(
      selectInput("gene", "Gene of Interest:", 
                  row.names(contrastdf)),#choices=fclist)
      hr(),
      helpText("Gene:"),
      helpText(verbatimTextOutput("geneInfo")),
      helpText("Maximum Fold Change Threshold:" ),
      helpText(verbatimTextOutput("sigLev"))
    ),
    # Create a spot for the barplot
    mainPanel(
      plotOutput("slicePlot"),
      plotOutput("volcanoPlot")
    )
  )
)

sa<-shinyApp(ui=ui, server=server)
runApp(sa, launch.browser=T)

