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
  args<-c("~/GitHub/Lab_data/ica_biofilm_qRT/aggregated_ica_bifilm_qRT.csv",
          "~/GitHub/Lab_data/std_curves_2016.csv",
          "split", 
          "rpoC","", "stat")
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
rawdata<-rawdata[rawdata$X !="#" &
                   rawdata$Cq!="NaN",colnames(rawdata)[colnames(rawdata)!="X"]]
std_curves<-read.csv(args[2][], stringsAsFactors = F)
merged<-merge(rawdata,y = std_curves[,c("Target", "m", "b")], by.x="Target", by.y="Target", all.x =T)
merged<-merged[!merged$Content %in% grep("Std.*",merged$Content , value=T ),]
#average
merged<-merged %>%
  group_by(Content, Target, Sample) %>%
  mutate(mean=mean(Cq))
merged<-merged %>%
  group_by(Content, Target, Sample) %>%
  mutate(Transcripts=10^((mean-b)/m))

kill_cols<-c("Well","mean", "Fluor", "Cq", "m","b")
for (i in kill_cols){
  merged[,i]<-NULL
}
merged<-as.data.frame(merged[!duplicated(merged),]) #cast from grouped dataframe to normal once
normlist<-c("rpoC")
for (i in normlist){
  normalizer<-merged %>%
    filter(Target==i) %>%
    mutate(norm_transcripts=Transcripts)
  normalizer<-
    normalizer[, colnames(normalizer)[!colnames(normalizer)%in%grep("_normalized_tr",colnames(normalizer),value=T)]]
  normalizer$Target<-NULL
  normalizer$Transcripts<-NULL
  colnames(normalizer)<-
    gsub("norm_transcripts",paste(i,"_normalized_transcripts", sep=""), colnames(normalizer))
  merged<-merge(merged, normalizer, by=c( "Content", "Sample"), all.x = T)
  
}
for (i in grep("_normalized_tr",colnames(merged),value=T)){
 ratios<-merged[,"Transcripts"]/merged[,i]
 merged<-cbind(merged,ratios)
 colnames(merged)<-
   gsub("ratios",paste(i,"_ratio", sep=""), colnames(merged))
 
}
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
#merged$error_percentage<-merged$
  

write.csv(merged, "~/GitHub/R/uams1rnaseq/biofilm/ica_biofilm_qrt_output.csv")
merged$norm_transcripts<-merged[,grep(args[4],colnames(merged), value=T )[1]]
cols_to_keep<-colnames(merged)[!colnames(merged)%in%grep("_normalized_tr",colnames(merged),value=T)]
merged<-merged[, cols_to_keep ]
               
merged<-merged[merged$rep!=args[5] & merged$condition!=args[6],]
merged<-merged %>%
  group_by(Content, Target, condition, strain, rep)%>%
  mutate(normalized_transcripts=mean(Transcripts/norm_transcripts))

#####  Stats!
for_stats<-as.data.frame(merged[merged$Content!="NRT" & merged$rep!="3",c("strain", "normalized_transcripts", "Target")])
for_stats<-for_stats[for_stats$Target%in%c("icaB_start","icaA"),]

h<-pairwise_sigs2(data = for_stats, strain = "strain", 
                  val = "normalized_transcripts", compareto = "813")


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

#unique<-unique[unique$Target!="16s",]
limits <- aes(ymax = normalized_transcripts_average + normalized_transcripts_SEM,
              ymin=normalized_transcripts_average - normalized_transcripts_SEM)

dodge <- position_dodge(width=0.9)

unique$percentage_error<-100*(unique$normalized_transcripts_SEM/unique$normalized_transcripts_average)
unique$error_message<-ifelse(unique$percentage_error>=30,"SEM > 30%", "")


normalized_transcript_plot<-
  ggplot(unique, aes(x=strain, y=normalized_transcripts_average, fill=condition, label=error_message))+
  geom_bar(position=dodge,stat="identity",fill="#5b5b5b", colour="#5b5b5b")+
  facet_wrap(~Target,scales= "free_y" ,nrow = 3)+
  geom_text(vjust = 0, nudge_y = 0.5, angle=90, color="red")+
  geom_errorbar(limits, position=dodge, width=.25)
normalized_transcript_plot
##############

source("~/GitHub/R/uams1rnaseq/ggplot_themes.R")

ica_min<-1640
rpoC_wt<-201835.987
min_detection<-ica_min/rpoC_wt

subset_normalized_transcript_plot<-
  ggplot(as.data.frame(unique[!unique$Target %in% c("rpoC", "icaB_start", "icaB_end"),]),
         aes(x=strain, y=normalized_transcripts_average, fill=condition,label=error_message))+
  geom_bar(position=dodge,stat="identity",fill="#5b5b5b", colour="#5b5b5b")+
  facet_wrap(~Target,scales= "free_y" ,nrow = 3)+
  scale_x_discrete(labels=expression( "UAMS-1", "codY","icaR","codY icaR"))+
  labs(x="Strain")+
  geom_errorbar(limits, position=dodge, width=.25)+
  theme_assay+
#  coord_cartesian(ylim=c(0,1.5))+
  scale_fill_manual(values=c("#5b5b5b", "#5b5b5b","#5b5b5b","#5b5b5b","#5b5b5b")) +
  scale_y_continuous(name = expression(paste(italic("icaA "),"Transcript Abundance")),
                     expand = c(0,0), limits = c(0, 1.1))
  
  
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






