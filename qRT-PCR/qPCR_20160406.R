#!/usr/bin/Rscript
####       Processing qRT-PCR data      ####
####          Nick Waters 2016          ####
####              20160315              ####
#install.packages("dplyr")
# usage: cqsummary.csv, standard_curves.csv, "split"__or_key_file, normalizing_gene, rep_to_remove, condition_to_remove
neededPackages <- c("dplyr","ggplot")
for (i in neededPackages){
    if(! i %in% installed.packages()){
        install.package(i)
        library(i)
    } else{
        library(i)
    }
}
source("~/GitHub/R/uams1rnaseq/utilsV1.R")
DEBUG=T
SORT_CONDITIONS=F
MAKE_OUTPUT=F
test="hi"
print("USAGE:Rscript qRT.R qrt_data.csv standard_curves.csv sample_key.csv|split normalizer_gene rep_to_ignore condition_to_ignore")

if(DEBUG){
  args<-c("~/20160406/qPCR_aggregated.csv",
          "~/GitHub/Lab_data/dummy_standards.csv",
          "split", 
          "Adapters","2", "Stock")
} else{
  args<-commandArgs(T)
  if (args[1]=="-h"|args[1]=="--help"){
      print("use this script to analyze qRT data.  Best used interactively, honestly.  I'm not proud of that fact, but hey, what can you do?Try again")
      stop("one day I'll figure out how to get rid of this error")
  }
  #test if all exist
  if (!dir.exists(args[1])){}
  #test headers
}
#
rawdata   <-read.csv(args[1], stringsAsFactors = F, header = T)
rawdata<-rawdata[!rawdata$X =="#" &
                   rawdata$Cq!="NaN",colnames(rawdata)[colnames(rawdata)!="X"]]
std_curves<-read.csv(args[2][], stringsAsFactors = F)
merged<-merge(rawdata,y = std_curves[,c("Target", "m", "b")], by.x="Target", by.y="Target", all.x =T)
merged<-merged[!merged$Content %in% grep("Std.*",merged$Content , value=T ),]
print(table(merged$Target))
print(table(merged$Sample))
#average
merged<-merged %>%
  group_by(Content, Target, Sample) %>%
  mutate(mean=mean(Cq))
merged<-merged %>%
  group_by(Content, Target, Sample) %>%
  mutate(Transcripts=10^((mean-b)/m))


merged<-merged[!is.na(merged$Transcripts),]
kill_cols<-c("Well","mean", "Fluor", "Cq", "m","b", "Sq", "SQ")
for (i in kill_cols){
  merged[,i]<-NULL
}
merged<-as.data.frame(merged[!duplicated(merged),]) #cast from grouped dataframe to normal once

normlist<-c(args[4], "rpoC")
for (i in normlist){ #  this subsets out all the targets in normlist (line above) and aligns them by sample for comparison
  normalizer<-merged %>%
    filter(Target==i) %>%
    mutate(normtranscripts=Transcripts)
  normalizer<-
    normalizer[, colnames(normalizer)[!colnames(normalizer)%in%grep("_transcripts",colnames(normalizer),value=T)]]
  normalizer$Target<-NULL
  normalizer$Transcripts<-NULL
  colnames(normalizer)<-
    gsub("normtranscripts",paste(i,"_transcripts", sep=""), colnames(normalizer))
  merged<-merge(merged, normalizer, by=c( "Content", "Sample"), all.x = T)
  
}
# 
for (i in grep("_transcripts",colnames(merged),value=T)){ # this
 ratios<-merged[,"Transcripts"]/merged[,i]
 merged<-cbind(merged,ratios)
 colnames(merged)[grep("ratios",colnames(merged))]<-paste("transcripts_normalized_by_", gsub("(.*)_(.*)", "\\1",i), sep='')
}


############################  Normalize by input pool
inputnorm="Stock__A"
for (i in c(normlist)){
  normalizer<-merged %>%
    filter(Sample==inputnorm) %>%
    mutate(norm_stock_per=transcripts_normalized_by_rpoC) %>%
    select(Content, Sample, Target,norm_stock_per)
  normalizer$Sample<-NULL
  colnames(normalizer)[grep("norm_stock_per", colnames(normalizer))]<-
    paste("norm_stock_per_",i, sep="")
  merged<-merge(merged, normalizer, by=c( "Content", "Target"), all.x = T)
  
}

for (i in normlist){ # this makes the ratio for output/input
  normalized_column<-grep(paste("by_",i,sep=''),colnames(merged), value=T)
  ratios<-merged[,normalized_column]/merged[,grep(paste("norm_stock_per_",i,sep=""), colnames(merged))]
  merged<-cbind(merged,ratios)
  colnames(merged)[grep("ratios",colnames(merged))]<-
    paste("target_per_",i,"_output_input_ratio",sep="")
}
############################



#  get sample names assigned
if(args[3] != "split"){
  samples<-read.csv(args[3], stringsAsFactors = F)
  merged<-merge(merged, samples, by="Sample")
} else{
  merged$name<-merged$Sample
}
merged$rep<-gsub("(.*)_(.*)_(.*)", "\\3", merged$name)
merged$condition<-ifelse(all(gsub("(.*)_(.*)_(.*)", "\\2", merged$name)==""), "grouped",
                             gsub("(.*)_(.*)_(.*)", "\\2", merged$name))
if(SORT_CONDITIONS){
  print("sorting condition; ensure sorting levels are correct!")
  merged$condition<-factor(merged$condition, levels = c("early","mid","late","stat"))
} 
merged$strain<-gsub("(.*)_(.*)_(.*)", "\\1", merged$name)
merged$name<-NULL
merged$Sample<-NULL

#  WRITE AND SUBSET

#output<-"~/Desktop/ica_qRT_results/"
#dir.create(output)
#setwd(output)

#  transcript error percentage
expression_normalizer<-"rpoC"
enrichment_normalizer<-"Adapters"
#magick
# (merged$enrichment_ratio<-(merged$transcripts_normalized_by_rpoC)*
#   (merged$transcripts_normalized_by_Stock__A))
(merged$enrichment_ratio<-(merged$target_per_Adapters_output_input_ratio))

### quick and dirty
merged$strain<-factor(merged$strain,levels = c( "0", "0.4","2","10", "50", "Stock") )

quick<-merged%>%
  filter(
#    rep!="B",
#    rep!="A",
    !strain %in% c("Stock"),#, "50"),
    !Target %in% c("Adapters", "ilvD_CBS_u","icaB_CBS_u"))

ggplot(as.data.frame(quick),#[quick$Target!="rpoC",]),
       aes(x=strain, y=enrichment_ratio, fill=condition))+
  #,label=error_message))+
  geom_bar(position=dodge,stat="identity",fill="#5b5b5b", colour="#5b5b5b")+
  facet_wrap(~Target,scales= "free_y" ,nrow = 3)+
#  scale_x_discrete(labels=expression( "0", "0.4","2","10"))+
#  geom_errorbar(limits, position=dodge, width=.25)+
  labs(x="CodY (nM)")


write.csv(merged, "~/20160406/qPCR_prelim_output.csv")
merged$norm_transcripts<-merged[,grep(args[4],colnames(merged), value=T )[1]] #subset based on args
cols_to_keep<-colnames(merged)[!colnames(merged)%in%grep("_normalized_tr",colnames(merged),value=T)]
merged<-merged[, cols_to_keep ]
               
merged<-merged[merged$rep!=args[5] & merged$condition!=args[6],]
merged<-merged %>%
  group_by(Content, Target, condition, strain, rep)%>%
  mutate(normalized_transcripts=mean(Transcripts/norm_transcripts))

#####  Stats!
for_stats<-as.data.frame(merged[merged$Content!="NRT" & merged$rep!="3",c("strain", "normalized_transcripts", "Target")])
for_stats<-for_stats[for_stats$Target%in%c("icaB_CBS","ilvD_CBS"),]

h<-pairwise_sigs2(data = for_stats, strain = "strain", 
                  val = "normalized_transcripts", compareto = "50")


##### end stats

#get transcipt ratio
averaged<-merged %>%
  group_by(Content, Target, condition, strain)%>%
  mutate(normalized_transcripts_average=mean(normalized_transcripts))%>%
  mutate(normalized_transcripts_SEM=(sd(normalized_transcripts)/(sqrt(n()))))

averaged$Transcripts<-NULL
averaged$rep<-NULL
averaged$norm_transcripts<-NULL
averaged$normalized_transcripts<-NULL
RT<-averaged[averaged$Content!="NRT", colnames(averaged)[colnames(averaged)!="Content"]]
unique<-RT[!duplicated(RT),]

unique<-unique[complete.cases(unique),]

#unique<-unique[unique$Target!="16s",]
limits <- aes(ymax = normalized_transcripts_average + normalized_transcripts_SEM,
              ymin=normalized_transcripts_average - normalized_transcripts_SEM)

dodge <- position_dodge(width=0.9)

unique$percentage_error<-100*(unique$normalized_transcripts_SEM/unique$normalized_transcripts_average)
unique$error_message<-ifelse(unique$percentage_error>=30,"SEM > 30%", "")



unique$strain<-factor(unique$strain,levels = c( "0", "0.4","2","10", "50", "Stock") )
normalized_transcript_plot<-
  ggplot(unique, aes(x=strain, y=normalized_transcripts_average, fill=condition, label=error_message))+
  geom_bar(position=dodge,stat="identity",fill="#5b5b5b", colour="#5b5b5b")+
  facet_wrap(~Target,scales= "free_y" ,nrow = 3)+
  geom_text(vjust = 0, nudge_y = 0.5, angle=90, color="red")+
  geom_errorbar(limits, position=dodge, width=.25)
normalized_transcript_plot
##############

source("~/GitHub/R/uams1rnaseq/ggplot_themes.R")

# ica_min<-1640
# rpoC_wt<-201835.987
# min_detection<-ica_min/rpoC_wt

subset_normalized_transcript_plot<-
  ggplot(as.data.frame(unique[!unique$strain %in% c("Stock"),]),
         aes(x=strain, y=enrichment_ratio, fill=condition,label=error_message))+
  geom_bar(position=dodge,stat="identity",fill="#5b5b5b", colour="#5b5b5b")+
  facet_wrap(~Target,scales= "free_y" ,nrow = 3)+
#  scale_x_discrete(labels=expression( "0", "0.4","2","10"))+
#  geom_errorbar(limits, position=dodge, width=.25)+
  labs(x="CodY (nM)")
#  theme_assay
#  coord_cartesian(ylim=c(0,1.5))+
  # scale_fill_manual(values=c("#5b5b5b", "#5b5b5b","#5b5b5b","#5b5b5b","#5b5b5b")) +
  # scale_y_continuous(name = expression(paste(italic("icaA "),"Transcript Abundance")),
  #                    expand = c(0,0), limits = c(0, 1.1))
  # 
  # 
subset_normalized_transcript_plot

##############
baseline<-unique %>%
  filter(strain=="337", condition=="grouped") %>%
  mutate(baseline_transcripts=normalized_transcripts_average)
baseline<-baseline[,colnames(baseline)[!colnames(baseline) %in%
                                         c("normalized_transcripts_average","percentage_error", "error_message","strain","normalized_transcripts_SEM")]]
#  insert limit of detection for fold changes
baseline[baseline$Target=="icaA", "baseline_transcripts"]<-min_detection
unique[unique$Target=="icaA" &
         (unique$strain=="337" |
         unique$strain=="814"), "normalized_transcripts_average"]<-min_detection


fold_change_baseline<-merge(unique,baseline, by=c("Target" ,  "condition" ), all.x=T)
fold_change<-fold_change_baseline %>%
  group_by(Target, condition) %>%
  mutate(fold_change_from_WT=(normalized_transcripts_average/baseline_transcripts))
         #  for icaA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fold_change$fold_change_from_WT_error=fold_change$fold_change_from_WT*sqrt((fold_change$normalized_transcripts_SEM/fold_change$normalized_transcripts_average)^2+
                                          (fold_change[fold_change$strain=="337" &fold_change$Target=="icaA","normalized_transcripts_SEM"]/fold_change$baseline_transcripts)^2)
fold_change$normalized_transcripts_average<-NULL
fold_change$normalized_transcripts_SEM<-NULL
fold_change$baseline_transcripts<-NULL
fold_change<-fold_change[!duplicated(fold_change),]


fold_limits <- aes(ymax = fold_change_from_WT + fold_change_from_WT_error,
              ymin=fold_change_from_WT)

fold_dodge <- position_dodge(width=0.9)

fold_change_plot<-
  ggplot(fold_change[fold_change$Target=="icaA",], # +.1 is to overlap error bars
         aes(x=strain, y=fold_change_from_WT+.1, fill=condition))+ #remember to bunp up dot in actual figure
  geom_errorbar(fold_limits, position=fold_dodge, width=.25)+
  geom_bar(position=fold_dodge,stat="identity",fill="#5b5b5b", colour="#5b5b5b")+
  scale_x_discrete(labels=expression( "UAMS-1", "codY","icaR","codY icaR"))+
  # scale_y_continuous(#breaks=c(1,100,500,1000),
  #               name = expression(paste(italic("icaA "),"Transcript Fold Change")),
  #               expand = c(0,0), limits = c(0, 140))+
  scale_y_log10(#breaks=c(1,100,100),
                name = expression(paste(italic("icaA "),"Transcript Fold Change")),
                expand = c(0,0), limits = c(1, 500))+
  labs(x="Strain")+
  theme_assay+
#  coord_cartesian(ylim=c(0,20))+
  scale_fill_manual(values=c("#5b5b5b", "#5b5b5b","#5b5b5b","#5b5b5b","#5b5b5b"))

fold_change_plot

if(MAKE_OUTPUT){
  
  
  write.csv(fold_change, "fold_change_output.csv")
  pdf(file = "normalized_fold_changes.pdf",width = 5, height = 5)
  normalized_transcript_plot
  subset_normalized_transcript_plot
  fold_change_plot
  dev.off()
}
save.image("~/Desktop/ica_qRT_results/data.RData")

#  for actual, use with manuscript companion.R
######################################## figure 6
pdf(file = paste(new_path,"figs6a.pdf", sep=""),width = 5*(4.4/2)/2.8, height = 5)
fold_change_plot
dev.off()






