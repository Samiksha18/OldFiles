#  LAC qRT DAta
# Nicholas Waters, November, 2015
#  for processing qRT data
#  this assumes that you have manually check the negative control values and 
#  they are withing the expected range

#  The raw data should be the summary page from cfX manager.  
#  Add a "plate" column if you run multiple plates
#  Use a new "plate" for each batch of cDNA, and have a new normalizer run.  So if 
#  didnt do that, treat them as a single plate

#  This program expects the following headers:
# "plate"   "Well"    "Fluor"   "Target"  "Content" "Sample"  "Rep"     "Cq" 


library(reshape2)
library(dplyr)
library(ggplot2)
setwd("~/GitHub/R/qPCR/")
source("../uams1rnaseq/utilsV1.R")

# read in flat file and key
qRT<-read.csv("../../../Desktop/pulldown_20151120/20151120_050642_CT001651_RPOC ILVD  -  Quantification Summary_0.csv")
curves<-read.csv("../../Lab_data/standards_post_nov15.csv", stringsAsFactors = F)

qRT$Rep<-gsub("(.*)_(.*)", "\\2", qRT$Sample)
#qRT$Sample<-gsub("(.*)_(.*)", "\\1", qRT$Sample)

#  gene to normalize 
housekeep<-"rpoC"
#sample to calculate fold changes from
fold_norm<-"DNA"

# plate group?
qRT$plate<-qRT$Rep
plate="plate"   #  NULL #or null


# test the headers
expected<-c("plate", "Well", "Fluor", "Target", "Content", "Sample", "Rep", "Cq" )

if(!all(unlist(lapply(expected, function(x){ match(x, colnames(qRT), nomatch = 0) > 0})))){
  print("you've got column issues, buddy!")
}
# prepare control values -------------------------------------------------------------------
housekeep_wells<-qRT[qRT$Target==housekeep,]
# get rid of neg
housekeep_wells<-housekeep_wells[housekeep_wells$Content!="Neg Ctrl",]
#get control Cq means  (ignor the SEM calculations)
housekeep_sum<-summarySE(data = housekeep_wells, measurevar = "Cq", biorepvar = "Fluor", 
                      groupvars = c("Target", plate,
                                    "Sample", "Rep"), 
                      na.rm = T, conf.interval = .95)
print("Enter 'print(housekeepsum)' to view the standard deviations of the normalization sample")

housekeep_sum_reduced<-housekeep_sum[,c("Target", "Sample", plate,
                              "Rep", "mean")]
housekeep_sum_reduced$housekeep_mean_Cq<-housekeep_sum_reduced$mean
housekeep_sum_reduced$mean<-NULL
housekeep_merged<-merge(x=housekeep_sum_reduced, y=curves, by="Target", all.x=T)
housekeep_merged$housekeep_gene_transcripts<-10^((housekeep_merged$housekeep_mean_Cq-housekeep_merged$Y.Intercept)/housekeep_merged$Slope)
housekeep_reduced<-housekeep_merged[,c("Sample", plate, 
                             "Rep", "housekeep_gene_transcripts")]


 # Actual Calcs ------------------------------------------------------------

qRT$Sample<-as.character(qRT$Sample)
#add in STD values
with_std<-merge(x=qRT, y=curves, by="Target", all.x=T)
# add in norm values
with_housekeep<-merge(with_std, housekeep_reduced, by=c("Sample", plate, 
                                              "Rep"), all.x=T)
# get rid of neg
noNeg<-with_housekeep[with_housekeep$Content!="Neg Ctrl",]
# cal transcripts and transcripts for norm
noNeg$transcripts<-10^((noNeg$Cq-noNeg$Y.Intercept)/noNeg$Slope)
#Transcripts / normalizer transcripts
noNeg$normalized_transcripts<-noNeg$transcripts/noNeg$housekeep_gene_transcripts
# get rid of extraneous info
#noExtra<-noNeg[,!colnames(noNeg)%in% c("Cq", "Fluor", "Well", colnames(curves))]
# get Stats and averages
summary<-summarySE(data = noNeg, measurevar = "normalized_transcripts",
                   biorepvar = "Rep",groupvars = c("Target", "Sample"),
                   na.rm = T,conf.interval = .95)


summary$fold_norm_transcripts<-as.numeric(summary[summary$Sample==fold_norm,"transcripts"][1,]) 

# for (i in unique(summary$plate)){#print(i)}
#   for (i in length())
#   if(i==summary$fold_norm_transcripts<-summary[summary$Sample==fold_norm & summary$plate==i, 
#                                      "transcripts"]
# }
# #  fold enrichment
# unk/rpoc
# _____________
# Inp_dna /rpoC
summary$foldEnrichment<-(summary$normalized_transcipts/summary$housekeep_gene_transcripts)/
                        (summary$fold_norm_transcripts/summary$housekeep_gene_transcripts)

# get Stats and averages
summary<-summarySE(data = noNeg, measurevar = "normalized_transcripts",
                   biorepvar = "Rep",groupvars = c("Target", "Sample"),
                   na.rm = T,conf.interval = .95)

for (i in unique(noNeg$plate)){
  noNeg$fold_norm_transcripts<-noNeg[noNeg$Sample==fold_norm & noNeg$plate==i, 
                                     "mean"]
}
summary<-summary[,!colnames(summary)%in% c("Fluor", "R.2", "Well","Content", "Cq",
                                           "Efficiency..", "Slope","Y.Intercept", 
                                           "transcripts", "techrep_sd", "N",
                                           "housekeep_gene_transcripts","plate")]
fold_norm_key<-data.frame(summary[summary$Sample==fold_norm, c("Target", "mean")])
fold_norm_key$comp_mean<-fold_norm_key$mean
fold_norm_key$mean<-NULL
summary_merge<-merge(summary,fold_norm_key, by="Target", all.x=T )
summary_merge$foldChange<-round(summary_merge$mean/summary_merge$comp_mean, digits = 2)
summary_merge$ci<-NULL
summary_merge$se_perc<-(summary_merge$se/summary_merge$mean)*100
summary_merge<-summary_merge[order( summary_merge$Target, summary_merge$Sample), ]

summary_merge
#Clean up data

#cleaned<-summary_merge[summary_merge$Target!="ilvD",]
cleaned<-summary_merge
#write out data

##  plot like no one's watching
#fold change
ggplot(cleaned, aes(x=Sample, y=foldChange, fill=Target))+
  facet_grid(~Target, scales = "free")+
  geom_bar(stat="identity", position="dodge")+
  scale_x_discrete(limits=c("687","758", "746"))+
  labs(y="Fold Change compared to 687", title="Fold Change")


limits <- aes(ymax = mean + se, ymin=mean - se)
dodge <- position_dodge(width=0.9)

ggplot(cleaned, aes(x=Sample, y=mean, fill=Target))+
  facet_grid(~Target, scales = "free")+
  geom_bar(stat="identity", position="dodge")+
  geom_errorbar(limits, position=dodge, width=0.25)+
  scale_x_discrete(limits=c("687","758", "746"))+
  labs(y="transcripts/rpoC transcript", title="Normalized Transcripts")



