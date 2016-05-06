
library(XLConnect)               # load XLConnect package 
library(ggplot2)
library(dplyr)
library(reshape2)
library(pgirmess)    #install.packages("WriteXLS")
library(WriteXLS)    #install.packages("pgirmess")
testPerl()  ##  needed for writing xls
###Expected Headers: Strain = "Strain", slope ="Val", and replicate/dilution = "Dil"
df<-read.csv('../../Lab_Data/nuc_assay_agg.csv', header=T,sep=",")
tdf<-df[,1:5]
head(tdf)
#colnames(tdf)<-c("strain", "A", "B", "C" )
mdf<-melt(tdf, value.name = "val", id.vars = c("strain", "sort"), variable.name = "rep")
head(mdf)
mdf<-mdf[complete.cases(mdf),]
#mdf$strain <- factor(mdf$strain, levels = mdf$sort)

p<-ggplot(mdf, aes(strain, val))+geom_boxplot()+
  scale_x_discrete(limits=c("nuc", "WT", "R61K", "R61H", "R61E",
                            "G129D", "codY", "sarA"))
p
# no sarA, no nuc
samp<-subset(mdf, strain != "nuc")
samp<-subset(samp, strain != "sarA")
#ANOVA 
anova<-aov(data=samp, val~strain)
anova
#####  homoscedasticity
#install.packages("car")
library(car)
leveneTest(data=samp, val~strain)  ##check for homoscedasticity
shapiro.test(anova$residuals)    #### check for normality of residuals 
anovaLog<-aov(data=samp, log(val)~strain)
leveneTest(data=samp, log(val)~strain)  ##check for homoscedasticity, log transformed
shapiro.test(anovaLog$residuals)    #### check for normality of residuals log transformed

#Tukey
tuk<-TukeyHSD(anova)
tuk
tukdf<-as.data.frame(tuk$strain)
names(tukdf)[names(tukdf)=="p adj"] <- "padj"
TukeyHSD<-tukdf
sig<-subset(tukdf, padj<= 0.05)
sig 

sigp<-ggplot(samp, aes(strain, val))+geom_boxplot()+
  scale_x_discrete(limits=c("WT", "R61K", "R61H", "R61E",
                            "G129D", "codY"))
sigp


####   Add asterisks ala :
#     http://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
label.df <- data.frame(strain = c("R61E",
                                 "G129D",
                                 "codY"), #"WT", "R61K"
                       Value = c(12, 13, 17))

sigp + geom_text(data = label.df, aes(x = strain, y= Value), label = "***")


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