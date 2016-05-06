#stupid fat microarray
# first time only ---------------------------------------------------------

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("saureuscdf")


# load libraries ----------------------------------------------------------


#biocLite("org.Dm.eg.db")
#biocLite("drosophila2.db")
library(affy)
library(affyPLM)
library(saureuscdf)
#library(org.Dm.eg.db)
#library(drosophila2.db)
setwd("~/GitHub/R/Majerczyk_microarray/")

# read the data
affydata <- ReadAffy()
affydata

# raw expression data
ed <- exprs(affydata)
sampleNames(affydata)<-(c("WT_1", "WT_2", "agr_1", "agr_2",
                           "codY_1", "codY_2", "codY_3", "codYagr_1", "codYagr_2"))
samp <- sampleNames(affydata)
probes <- featureNames(affydata)

ed[1:10,]
probes[1000:1100]
samp

# 2 Normalizing Data   
#
# The Affy package has implementations of a number of normalization methods
# for single-channel arrays. This includes (among others):
#   - mas5() - Affymetrix's own normalization program
#   - rma() - 'Robust Multi-Chip' average
#   - gcrma() - A bias-corrected RMA
# GCRMA is good but takes significantly longer than RMA, so RMA is the
# most commonly used
system.time(nvals <- rma(affydata))

# normalised expression data
ned <- exprs(nvals)

nsamp <- sampleNames(nvals)
nprobes <- featureNames(nvals)

ned[1:10,]
nprobes[1:10]
nsamp

# the normalized data is on the log-scale

# 3 visualize the data
# We can plot a case vs a control
#    1 GSM277408     Control Abd_1
#    7 GSM277414     Logjam Abd_1

# we can figure out which columns these are from colnames
colnames(ned)
plot(ned[,"WT_1"], ned[,"WT_2"])
plot(ned[,"codY_1"], ned[,"codY_2"])
plot(ned[,"WT_2"], ned[,"codY_2"])

# tidy it up a bit
plot(ned[,"WT_1"], ned[,"codY_1"], pch=".",
     xlab="WT_1", ylab="codY_1",
     main="CodY mutant vs WT")

# add a line of y=x
abline(0,1,col="blue")

# plot lines at two fold up and down regulation
# hint: two fold upregulation is 2, so we draw the intercept at log(2)
# hint: two fold downregulation is .5, so we draw the intercept at log(.5)
abline(log(2),1,col="red")
abline(log(.5),1,col="red")

# use the identify function to label points
#identify(ned[,"WT_1"], ned[,"codY_1"], nprobes)

# press escape to exit on windows, or right click to exit on Linux

# 4 Get some annotation
# get the package that contains the annotation for this array
csv<-
  read.csv("S_aureus.na35.annot.csv", skip = 15)

anno_table<-data.frame(row.names=nprobes, desc=paste(csv$Probe.Set.ID, csv$Target.Description, csv$SwissProt, sep="|"))
anno_table<-data.frame(row.names=nprobes, desc=paste(csv$Probe.Set.ID, csv$Transcript.Assignments, sep="|"))
#Transcript Assignments
# 5 Process and Export the data
# ned is our normalized expression matrix

# merge with our annotation
aned <- merge(ned, anno_table, by="row.names", all.x=TRUE, sort=FALSE)
row.names(aned)<-aned$desc
aned$desc<-NULL
aned$Row.names<-NULL
write.csv(aned, "micro_array_new_anno.csv", row.names=FALSE)


# save data
save.image()




aned<-ned

# do you wanna build a  Expression Set? #snowman?-----------------------------
condition<-gsub("(.*?)_(.*)","\\1", colnames(aned))
pData<-data.frame(row.names=colnames(aned), cond=condition)
metadata <- data.frame(labelDescription= c("gene of interest knocked out"),
                           row.names=c("condition"))
phenoData <- new("AnnotatedDataFrame",
                       data=pData, varMetadata=metadata)

experimentData <- new("MIAME", name="Nick Waters",lab="Brinsmade Lab",
                        contact="Dave Samuels;vits his fault",
                        title="reworking Charlotte's old data with decent annotations",
                        abstract="for the UAMS-1 RNA-seq paper",
                        url="www.lab.not.exist",
                        other=list(
                            notes="Created20151110"
                            ))
 exampleSet <- ExpressionSet(assayData=as.matrix(aned), phenoData=phenoData,
                                experimentData=experimentData,
                                annotation="saureus")
 featureNames(exampleSet)[100:115]
 # DE? ---------------------------------------------------------------------
 library(limma)
 
 
 design <- model.matrix(~ 0+factor(c(1,1,2,2,3,3,3,4,4)))
 colnames(design) <- c("WT", "agr", "codY", "codYagr")
 fit <- lmFit(exampleSet, design)
 contrast.matrix <- makeContrasts(agr-WT, codY-WT, codYagr-WT,codYagr-codY,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
 #A list of top genes differential expressed in group2 versus group1 can be obtained from
View(topTable(fit2, coef=2, adjust="BH", number = Inf, p.value = .05, lfc = 1))
 #The outcome of each hypothesis test can be assigned using
results <- decideTests(fit2, p.value = .01, lfc = 1)
# A Venn diagram showing numbers of genes significant in each comparison can be obtained from
vennDiagram(results)
 
plotMD(fit2, column=2)#, status=chrom, values=c("X","Y", "X|Y"))
#qqt(y=fit2$ , df.prior+fit$df.residual,pch=16,cex=0.2)

top<-topTable(fit2, coef=2, adjust="BH", number = Inf, p.value = .05, lfc = 1)
top_anno <- merge(top, anno_table, by="row.names", all.x=TRUE, sort=FALSE)


top_anno$P.Value<-round(top_anno$P.Value, 6)
top_anno$fc<-2^top_anno$logFC




write.csv(top_anno, "Mdata_lfc1_pval.05.csv")






