#  Pulldown qRT DAta
# Nicholas Waters, November, 2015
#  for processing qRT data
#  this assumes that you have manually check the negative control values and 
#  they are withing the expected range

#  The raw data should be the summary page from cfX manager.  
#  Add a "plate" column if you run multiple plates
#  Use a new "plate" for each batch of cDNA, and have a new normalizer run.  So if 
#  didnt do that, treat them as a single plate

library(reshape2)
library(dplyr)
library(ggplot2)
setwd("~/GitHub/R/qPCR/")
source("../uams1rnaseq/utilsV1.R")


# read in flat file and key
LACqRT<-read.csv("alecs_rot_pulldown.csv")
#____
LACqRT$Target<-gsub("rpoc", "rpoC",LACqRT$Target)
LACqRT$Rep   <-gsub("(.+)_(.)", "\\2",LACqRT$Sample)
LACqRT$Rep<-1
#LACqRT$Sample   <-gsub("(.+)_(.)", "\\1",LACqRT$Sample)
LACqRT$plate<-1


#____
curves<-read.csv("../../Lab_data/standards_post_nov15.csv", stringsAsFactors = F)
#  gene to normalize 
norm<-"rpoC"
#sample to calculate fold changes from

fold_norm<-"400_2"
# prepare control values -------------------------------------------------------------------
norm_wells<-LACqRT[LACqRT$Target==norm,]
# get rid of neg
norm_wells<-norm_wells[norm_wells$Content!="Neg Ctrl",]
#get control Cq means  (ignor the SEM calculations)

norm_sum<-summarySE(data = norm_wells, measurevar = "Cq", biorepvar = "Fluor", 
                    groupvars = c("Target", "Sample", "plate", "Rep"), 
                    na.rm = T, conf.interval = .95)
print("Enter 'print(normsum)' to view the standard deviations of the normalization sample")

norm_sum_reduced<-norm_sum[,c("Target", "Sample", "plate", "Rep", "mean")]
norm_sum_reduced$norm_mean_Cq<-norm_sum_reduced$mean
norm_sum_reduced$mean<-NULL
norm_merged<-merge(x=norm_sum_reduced, y=curves, by="Target", all.x=T)
norm_merged$norm_gene_transcripts<-10^((norm_merged$norm_mean_Cq-norm_merged$Y.Intercept)/norm_merged$Slope)
norm_reduced<-norm_merged[,c("Sample", "plate", "Rep", "norm_gene_transcripts")]


# Actual Calcs ------------------------------------------------------------

LACqRT$Sample<-as.character(LACqRT$Sample)
#add in STD values
with_std<-merge(x=LACqRT, y=curves, by="Target", all.x=T)
# add in norm values
with_norm<-merge(with_std, norm_reduced, by=c("Sample", "plate", "Rep"), all.x=T)
# get rid of neg
noNeg<-with_norm[with_norm$Content!="Neg Ctrl",]
# cal transcripts and transcripts for norm
noNeg$transcripts<-10^((noNeg$Cq-noNeg$Y.Intercept)/noNeg$Slope)
#Transcripts / normalizer transcripts
noNeg$normalized_transcripts<-noNeg$transcripts/noNeg$norm_gene_transcripts
# get rid of extraneous info
#noExtra<-noNeg[,!colnames(noNeg)%in% c("Cq", "Fluor", "Well", colnames(curves))]

# get Stats and averages
summary<-summarySE(data = noNeg, measurevar = "normalized_transcripts",
                   biorepvar = "Rep",groupvars = c("Target", "Sample"),
                   na.rm = T,conf.interval = .95)
fold_norm_key<-data.frame(summary[summary$Sample==fold_norm, c("Target", "mean")])
fold_norm_key$comp_mean<-fold_norm_key$mean
fold_norm_key$mean<-NULL
summary_merge<-merge(summary,fold_norm_key, by="Target", all.x=T )
summary_merge$foldChange<-round(summary_merge$mean/summary_merge$comp_mean, digits = 2)
#Clean up data

#write out data

#  plot like no one's watching








noNeg<-noNeg[!is.na(noNeg$Cq) ,]
tech_aves<-noNeg %>% group_by_("Rep", "plate", "Target", "Sample")%>% mutate(tech_ave=mean(Cq), sd=sd(Cq), length=length(Cq), sem=(sd(Cq))/sqrt(length(Cq))) 
tech_aves$Cq<-NULL
tech_aves$Well<-NULL
tech_aves$Fluor<-NULL
tech_aves$Content<-NULL
tech_aves<-tech_aves[!duplicated(tech_aves),]

tech_ave_std<-merge(x=tech_aves, y=curves, by="Target", all.x=T)
tech_ave_std$transcripts<-10^((tech_ave_std$tech_ave-tech_ave_std$Y.Intercept)/tech_ave_std$Slope)
rpoCdf<-tech_ave_std[tech_ave_std$Target=="rpoC", c("plate", "Rep", "Sample", "transcripts")]
rpoCdf$rpoC_transcripts<-rpoCdf$transcripts
rpoCdf$transcripts<-NULL

merged<-merge(tech_ave_std, rpoCdf, by=c("plate", "Rep", "Sample"), all.x=T)
merged$Sample<-as.character(merged$Sample)
merged$normed_transcripts<-merged$transcripts/merged$rpoC_transcripts)
#####  proceed with caution--  manual row removal
brnQ<-merged[merged$Target=="brnQ" & !merged$tech_ave>=30 ,]
brnQ$foldChange<-
  lukS<-merged[merged$Target=="lukS" & merged$plate=="B" ,]
brnQplot<-
  ggplot(brnQ, aes(x=Sample,y=normed_transcripts))+geom_bar(stat="identity")
ggplot(lukS, aes(x=Sample,y=normed_transcripts))+geom_bar(stat="identity")


