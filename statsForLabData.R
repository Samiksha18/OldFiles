#  Stats for lab data


library(XLConnect)   # load XLConnect package 
library(ggplot2)
library(dplyr)
library(reshape2)
library(pgirmess)    #install.packages("WriteXLS")
library(WriteXLS)    #install.packages("pgirmess")
library(car)
#install.packages("lubridate")
library(lubridate)
testPerl()  ##  needed for writing xls

###Expected Headers: Strain = "Strain", slope ="Val", and replicate/dilution = "Dil"
assayTitle1<-"Nuc assay"

## Format wide data
posControl<-"sarA"
df<-read.csv('../../Lab_Data/nuc_assay_agg.csv', header=T,sep=",")
tdf<-df[,1:5]
head(tdf)
mdf<-melt(tdf, value.name = "val", id.vars = c("strain", "sort"), variable.name = "rep")
head(mdf)
mdf<-mdf[complete.cases(mdf),]
mdf<-mdf%>%mutate(type=ifelse(!mdf$strain %in% posControl, "sample", "control"))

p<-ggplot(mdf, aes(reorder(strain, val), val))+
#  try(scale_x_discrete(limits=c("nuc","WT", "R61K", "R61H", "R61E", "G129D", "codY", "sarA")))+
  labs(x="Strain", y="Units of nuclease per OD600 per ml", title= assayTitle1)+
  geom_boxplot()+
  facet_grid(~type, scales = "free_x", space="free_x")
p
# no sarA, no nuc
samp<-subset(mdf, strain != "nuc")
samp<-subset(mdf, !strain %in% posControl)
anova<-aov(data=samp, val~strain) # ANOVA
anova
leveneTest(data=samp, val~strain)  ##check for homoscedasticity, should be greater than .05
shapiro.test(anova$residuals)    #### check for normality of residuals, should be more than .05

tukdf<-as.data.frame(TukeyHSD(anova)$strain)   #Tukey
names(tukdf)[names(tukdf)=="p adj"] <- "padj"
sig<-subset(tukdf, padj<= 0.05)

##  Add asterisks ala : http://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
mut<-samp%>% group_by(strain)%>% mutate(ave=mean(val))
mut<-mut[,!colnames(mut)=="Date"& !colnames(mut)=="val"]
colnames(mut)<-gsub("ave", "val", colnames(mut))
labeldf<-mut[!duplicated(mut$val),]
head(labeldf)
#  watch this:  this reduces to only significant ones, and allows for locating asterisks
sigloc<-merge(labeldf, data.frame(strain=gsub("-*f+", "", grep("fff", gsub("[cC]odY", "fff", row.names(sig)), value=T))), by = "strain", all.y=T)
p+ geom_text(data=sigloc, aes(x = strain, y= val+1), label = "***")
#
#
#
#
#

assayTitle2<-"Biofilm assay"
posControl<-NA
## format tall data
df<-read.csv('../../Lab_Data/4.2 Data Compilation 3 Replicates Revised 4.16_NWedits.csv', header=T,sep=",")
colnames(df)<-gsub("(.+)\\.", "\\1", colnames(df))
colnames(df)
mdf<-melt(df, value.name = "val", id.vars = c("Date", "Trial"), variable.name = "strain")
head(mdf)
mdf<-mdf[complete.cases(mdf),]
mut<-mdf%>% group_by(Date,strain)%>% mutate(ave=mean(val))
mut$Trial<-NULL
mut$val<-NULL
colnames(mut)<-gsub("ave", "val", colnames(mut))
mdf<-mdf%>%mutate(type=ifelse(!mdf$strain %in% posControl, "sample", "control"))

p<-ggplot(mdf, aes(reorder(strain, val), val))+
  #  try(scale_x_discrete(limits=c("nuc","WT", "R61K", "R61H", "R61E", "G129D", "codY", "sarA")))+
  labs(x="Strain", y="OD600", title= assayTitle2)+
  geom_boxplot()+
  facet_grid(~type, scales = "free_x", space="free_x")
p
# no sarA, no nuc
samp<-subset(mdf, !strain %in% posControl)
anova<-aov(data=samp, val~strain) # ANOVA
anova
leveneTest(data=samp, val~strain)  ##check for homoscedasticity, should be greater than .05
shapiro.test(anova$residuals)    #### check for normality of residuals, should be ess than .05? 

tuk<-TukeyHSD(anova)   #Tukey
tuk
tukdf<-as.data.frame(tuk$strain)
names(tukdf)[names(tukdf)=="p adj"] <- "padj"
sig<-subset(tukdf, padj<= 0.05)
 
##  Add asterisks ala : http://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
mut<-samp%>% group_by(strain)%>% mutate(ave=mean(val))
mut<-mut[,!colnames(mut)=="Date"& !colnames(mut)=="val"]
colnames(mut)<-gsub("ave", "val", colnames(mut))
labeldf<-mut[!duplicated(mut$val),]
head(labeldf)
#  watch this:  this reduces to only significant ones, and allows for locating asterisks
sigloc<-merge(labeldf, data.frame(strain=gsub("(.+-)(.+)", "\\2", grep("CodY", row.names(sig), value=T))), by = "strain", all.y=T)
p+ geom_text(data=sigloc, aes(x = strain, y= val+1), label = "***")

#
#
#
assayTitle3<-"protease assay"
df<-read.csv('../../Lab_Data/protease/protease_assay_agg.csv', header=T,sep=",")
#colnames(df)<-gsub("(.+)\\.", "\\1", colnames(df))
#colnames(df)<-gsub("X", "", colnames(df))
df<-subset(df, Date!="7/10/15")
colnames(df)<-gsub("(.+)_(.+)", "\\2", colnames(df))
colnames(df)
posControl<-c("sarA",  "Pro.K","SarA")
mdf<-melt(df, value.name = "val", id.vars = c("Date", "Trial"), variable.name = "strain")
head(mdf)
mdf<-mdf[complete.cases(mdf),]
mut<-mdf%>% group_by(Date,strain)%>% mutate(ave=mean(val))
mut$Trial<-NULL
mut$val<-NULL
colnames(mut)<-gsub("ave", "val", colnames(mut))
mdf<-mut[!duplicated(mut$val),]
mdf<-as.data.frame(mdf)%>%mutate(type=ifelse(!mdf$strain %in% posControl, "sample", "control"))


p<-ggplot(mdf, aes(reorder(strain, val), val))+
  #  try(scale_x_discrete(limits=c("nuc","WT", "R61K", "R61H", "R61E", "G129D", "codY", "sarA")))+
  labs(x="Strain", y="OD600", title= assayTitle3)+
  geom_boxplot()+
  geom_point()+
  facet_grid(~type, scales ="free", space="free_x")
p

samp<-subset(mdf, !strain %in% posControl)
anova<-aov(data=samp, val~strain)#ANOVA 
anova
leveneTest(data=samp, val~strain)  ##check for homoscedasticity
shapiro.test(anova$residuals)    #### check for normality of residuals 
anovaLog<-aov(data=samp, log(val)~strain)
leveneTest(data=samp, log(val)~strain)  ##check for homoscedasticity, log transformed
shapiro.test(anovaLog$residuals)    #### check for normality of residuals log transformed
tuk<-TukeyHSD(anova)
tuk
tukdf<-as.data.frame(tuk$strain)
names(tukdf)[names(tukdf)=="p adj"] <- "padj"
TukeyHSD<-tukdf
sig<-subset(tukdf, padj<= 0.05)
sig 

####   Add asterisks ala: http://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
mut<-samp%>% group_by(strain)%>% mutate(ave=mean(val))
mut<-mut[,!colnames(mut)=="Date"& !colnames(mut)=="val"]
colnames(mut)<-gsub("ave", "val", colnames(mut))
labeldf<-mut[!duplicated(mut$val),]
head(labeldf)
#  watch this:  this reduces to only significant ones, and allows for locating asterisks
sigloc<-merge(labeldf, data.frame(strain=gsub("(.+-)(.+)", "\\2", grep("odY", row.names(sig), value=T))), by = "strain", all.y=T)
p+ geom_text(data=sigloc, aes(x = strain, y= val+.1), label = "***")


#Kruskal (if assumtion of nornality is false) 
k<-kruskal.test(val ~ strain, data=samp)
k
kruskpost<- kruskalmc(val~strain, data=samp) 
kruskpost



posttest_p0.05_kruskal<-kruskpost$dif.com
ignore_nuc_and_sarA<-samp
WriteXLS(x =c("df", "mdf", "ignore_nuc_and_sarA", "TukeyHSD", "sig", "posttest_p0.05_kruskal"), 
         row.names = T,
         col.names = T,
         ExcelFileName = "stats_nuc_NW2010601_01.xlsx", 
         envir = .GlobalEnv)




