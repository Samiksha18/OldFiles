#!/usr/bin/Rscript
####       Processing qRT-PCR data      ####
####          Nick Waters 2015          ####
####              20150416              ####
# install.packages("XLConnect")
argv<-commandArgs()

library(XLConnect)               # load XLConnect package 
library(ggplot2)
library(dplyr)
library(reshape2)

#source("demoCONFIG.R")
# print(argv)
# print(paste("length of argv:", length(argv)))
# conff<-grep("CONFIG", argv, value=T)
# source(conff)
# 


df<-primerBind
df$val<-as.numeric(df$val)
df$Cycle<-as.numeric(df$Cycle)
#df$well<-as.character(df$well)
d<-ggplot(df, aes(x=Cycle, y=val, group=well))+geom_line()+facet_wrap(~ps)+ ggtitle("Standard Curve")
d
#####  "Estimate" Threshold (ie, flat-out poach..)
#Caution
standard_cqdf<-standard_cqdf[complete.cases(standard_cqdf$Cq),]
#End Caution
sample_well=standard_cqdf[standard_cqdf$Cq<28,][1,]
print(paste("calibrating threshold from well", sample_well$well))
sample_Cq=sample_well$Cq
sample_subset<-subset(df, well==sample_well$well)

##Best Guess (if not overridden)
testInt<-approx(y = sample_subset$val, x = sample_subset$Cycle, xout = sample_Cq)
thresh<- testInt$y
thresh<- ifelse(exists("threshOverride"), threshOverride, thresh)
if(exists("threshOverride")){print("Brace yourselves: someone has overridden the threshold...")}
print(paste("hmmmm... looks like the threshold is", thresh))


d+geom_hline(y=testInt$y)

##  Calulate Cq and Cq_mean values
dfs<-df %>%group_by(well) %>% mutate(Cq= approx(y=Cycle, x=val, xout=thresh)$y)
df<-dfs[complete.cases(dfs),]
print(paste(length(unique(dfs$well))-length(unique(df$well)), " well(s) did not cross the threashold"))
df<-df %>%group_by(rep_group,ps)%>% mutate(Cq_mean= mean(Cq))
df<-df %>%group_by(rep_group,ps)%>% mutate(Cq_StDev= sd(Cq))
df<-df %>%group_by(rep_group,ps)%>% mutate(Cq_StErrMean= sd(Cq)/sqrt(length(Cq)))

# create summarized dataframe
altsum<-df
altsum$Cycle<-NULL
altsum$val<-NULL
sumdf<-altsum[!duplicated(altsum), ]
# sumdf<-as.data.frame(cbind(df$ps, df$rep_group, df$Cq_mean, df$Cq_StDev, df$Cq_StErrMean))
# colnames(sumdf)<-c("ps", "rep_group", "Cq_mean", "Cq_StDev","Cq_StErrMean" )
# sumdf<-sumdf[!duplicated(sumdf), ]



# Make a column for true concentration, write out results
sumdf$Cq_mean<-as.numeric(sumdf$Cq_mean)
sumdf$Cq_StDev<-as.numeric(sumdf$Cq_StDev)
sumdf$Cq_StDev<-as.numeric(sumdf$Cq_StErrMean)
sumdf<-subset(sumdf, rep_group!="H")#        Remove 0's
sumdf<-sumdf%>%group_by(rep_group,ps) %>% mutate(conc=get(as.character(rep_group)))
write.csv(sumdf, paste(new_path,title, "_calculated_Cq.csv", sep=""))
sumdf$logConc<-log(sumdf$conc)


coefps1<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps1)))
try(coefps2<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps2))), silent = T)
try(coefps3<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps3))), silent = T)
try(coefps4<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps4))), silent = T)
try(coefps5<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps5))), silent = T)
try(coefps6<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps6))), silent = T)
# r^2
lrps1<-summary(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps1)))$r.squared
try(lrps2<-summary(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps2)))$r.squared, silent = T)
try(lrps3<-summary(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps3)))$r.squared, silent = T)
try(lrps4<-summary(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps4)))$r.squared, silent = T)
try(lrps5<-summary(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps5)))$r.squared, silent = T)
try(lrps6<-summary(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps6)))$r.squared, silent = T)


coefs<- grep("coefps.", ls(), value=T) 
coefdf<-suppressWarnings(cbind(get(coefs[1]), 
                                   try(get(coefs[2]), silent =T ), 
                                   try(get(coefs[3]), silent =T ),
                                   try(get(coefs[4]), silent =T ),
                                   try(get(coefs[5]), silent =T ),
                                   try(get(coefs[6]), silent =T )))
coefdf<-suppressWarnings(as.data.frame(matrix(as.numeric(unlist(coefdf)),nrow=nrow(coefdf))))
coefdf<-coefdf[, colSums(is.na(coefdf)) != nrow(coefdf)]

lrpss<- grep("lrps\\d", ls(), value=T) 
lrpsdf<-suppressWarnings(cbind(get(lrpss[1]), 
                               try(get(lrpss[2]), silent =T ), 
                               try(get(lrpss[3]), silent =T ),
                               try(get(lrpss[4]), silent =T ),
                               try(get(lrpss[5]), silent =T ),
                               try(get(lrpss[6]), silent =T )))
lrpsdf<-suppressWarnings(as.list(matrix(as.numeric(unlist(lrpsdf)),nrow=nrow(lrpsdf))))
#lrpsdf<-lrpsdf["NA"<-NULL]



coefdfr<-pexist
rsqr<-lrpsdf
colnames(coefdf)<-coefdfr
coefdf<-rbind(coefdf, rsqr)
row.names(coefdf)<-c("intercept", 'log(conc, 10)', "r^2")
write.csv(coefdf, paste(new_path, title,"_coefdf.csv", sep=""))

###  Double check-  does green match black????????
plot_lr<-function(df){
   pdf(file =  paste(new_path, title, "_stdCurve_plots.pdf",sep=""))
   par(mfrow = c(2,2))
   for(p in pexist){
     q<-ggplot(subset(df, ps==p), aes(y=Cq_mean, x=log(conc, 10)))+
       geom_point()+
       ggtitle(p)+
       geom_text(aes(x=6, y=21), label =paste("r^2:",coefdf[[p]][3]))+
       geom_smooth(aes(group=ps), method="lm", se=FALSE, color="green")+  #True
       geom_abline(intercept =coefdf[[p]][1], slope = coefdf[[p]][2])
      
       
     print(q) 
    }
   dev.off()
}
plot_lr(sumdf)


#########################################  from test files
#all<-bind_rows(wk1, wk2)

###  Make data long, remove nonsense
all<-melt(all, id.vars = c("well", "strain", "ps", "condition"), value.name = "val")
all<-all[complete.cases(all),]
all$name<-all$strain# all<-all %>% group_by(row.names(all)) %>% mutate(name=get(as.character(strain)))

all$val<-as.numeric(all$val)
try(all$cycle<-as.numeric(gsub("(\\D)(.+)", "\\2", all$variable)), silent=T)
#try(all$temp<-as.numeric(gsub("(.*)\\((.*)\\)(.*)", "\\2", all$name)), silent=T)
all<-all[all$ps != "",]
#try(all$str<-gsub("(.*)\\((.+)", "\\1", all$name), silent=T)

####  Nice quick visual check
pdf(file = paste(new_path, title,"_sample_visuals.pdf", sep=""),width = 8, height = 4)
p<-ggplot(all, aes(x=cycle, y=val, group=well, color=name))+
  geom_line()+
  facet_grid(condition~ps)
p
dev.off()



#  Technical replicate 
all<-all %>%group_by(well, ps) %>% mutate(Cq= approx(y=cycle, x=val, xout=thresh)$y)
all<-all %>%group_by(strain,ps,condition)%>% mutate(tech_Cq_mean= mean(Cq, na.rm=T))
all<-all %>%group_by(strain,ps,condition)%>% mutate(tech_Cq_StDev= sd(Cq))
all<-all %>%group_by(strain,ps,condition)%>% mutate(tech_Cq_StErrMean= sd(Cq)/sqrt(length(Cq)))


##  Summary Trimming
summ<-all
summ$variable<-NULL
summ$val<-NULL
summ$row.names<-NULL
summ$well<-NULL
summ$Cq<-NULL
summ["row.names(all)"]<-NULL
summ$cycle<-NULL
summ<-summ[!duplicated(summ), ]
##  Did you check the NRT condition? Its ok? Alrght, move along...
summ<-subset(summ, name !="Water")

# 
# ## Biological replicates
summ<-summ %>%group_by(strain, ps,condition) %>% mutate(bio_Cq_mean= mean(tech_Cq_mean))
summ<-summ %>%group_by(strain, ps,condition) %>% mutate(bio_Cq_StDev= sd(tech_Cq_mean))
summ<-summ %>%group_by(strain, ps,condition) %>% mutate(bio_Cq_StErrMean= sd(tech_Cq_mean)/sqrt(length(tech_Cq_mean)))
#summ$name<-paste(summ$strain, "(", summ$temp,")", sep="")
summ[is.na(summ)] <- 0
summ$tech_Cq_mean<-NULL
summ$tech_Cq_StDev<-NULL
summ$tech_Cq_StErrMean<-NULL
summ$strain<-NULL
summ<-summ[!duplicated(summ), ]


####  
print("Sanity Check: if true, NA's exist (or something much more sinister)")
duplicated(x = summ$bio_Cq_mean)

#if(all(summ$bio_))

#coefdf<-read.csv("coefdf.csv")
ps1sumdf<-subset(summ, ps==ps1)
try(ps2sumdf<-subset(summ, ps==ps2), silent=T)
try(ps3sumdf<-subset(summ, ps==ps3), silent=T)
try(ps4sumdf<-subset(summ, ps==ps4), silent=T)
try(ps5sumdf<-subset(summ, ps==ps5), silent=T)
try(ps6sumdf<-subset(summ, ps==ps6), silent=T)


ps1sumdf<-ps1sumdf%>%mutate(transcripts=10^((bio_Cq_mean-coefdf[,ps1][1])/(coefdf[,ps1][2])), 
                            transcript_error=10^((bio_Cq_mean-coefdf[,ps1][1])/(coefdf[,ps1][2]))*bio_Cq_StErrMean)
try(ps2sumdf<-ps2sumdf%>%mutate(transcripts=10^((bio_Cq_mean-coefdf[,ps2][1])/(coefdf[,ps2][2])), 
                                transcript_error=10^((bio_Cq_mean-coefdf[,ps2][1])/(coefdf[,ps2][2]))*bio_Cq_StErrMean), silent = T)
try(ps3sumdf<-ps3sumdf%>%mutate(transcripts=10^((bio_Cq_mean-coefdf[,ps3][1])/(coefdf[,ps3][2])), 
                                transcript_error=10^((bio_Cq_mean-coefdf[,ps3][1])/(coefdf[,ps3][2]))*bio_Cq_StErrMean), silent = T)
try(ps4sumdf<-ps4sumdf%>%mutate(transcripts=10^((bio_Cq_mean-coefdf[,ps4][1])/(coefdf[,ps4][2])), 
                                transcript_error=10^((bio_Cq_mean-coefdf[,ps4][1])/(coefdf[,ps4][2]))*bio_Cq_StErrMean), silent = T)
try(ps5sumdf<-ps5sumdf%>%mutate(transcripts=10^((bio_Cq_mean-coefdf[,ps5][1])/(coefdf[,ps5][2])), 
                                transcript_error=10^((bio_Cq_mean-coefdf[,ps5][1])/(coefdf[,ps5][2]))*bio_Cq_StErrMean), silent = T)
try(ps6sumdf<-ps6sumdf%>%mutate(transcripts=10^((bio_Cq_mean-coefdf[,ps6][1])/(coefdf[,ps6][2])), 
                                transcript_error=10^((bio_Cq_mean-coefdf[,ps6][1])/(coefdf[,ps6][2]))*bio_Cq_StErrMean), silent = T)

backsum<-bind_rows(ps1sumdf,ps2sumdf, ps3sumdf, ps4sumdf, ps5sumdf, ps6sumdf)
#backsum$str<-gsub("(.*)\\((.+)", "\\1", backsum$name)
## test contrast with subsets xx and yy
#xx<-subset(backsum, name=="337(25)")
# 
# contrast_rpoC<-function(inp){
#   contTrans<-unlist(inp[inp["ps"]=="rpoC","transcripts"])
#   return(contTrans)
# }
# 
# contrast_polC<-function(inp){
#   contTran<-unlist(inp[inp["ps"]=="polC","transcripts"])
#   return(contTran)
# }
if(exists("norm_primer1")){
  backsum<-backsum %>% 
    group_by(name) %>%  
    mutate(Norm1_transcripts= transcripts*(transcripts/(transcripts[ps==norm_primer1])))
} else{
  backsum<-backsum %>% 
    group_by(name) %>%  
    mutate(Norm1_transcripts= transcripts*1)
}


# backsum<-try(backsum %>%  ###  This is on its way out!!!!
#   group_by(name) %>%  
#   mutate(Norm2_transcripts= transcripts*(transcripts/(transcripts[ps==norm_primer2]))), silent=T)

#yy<-subset(backsum, ps=="rsaD")
if(exists("foldChange0")){
  backsum<-backsum %>% 
    group_by(ps) %>% 
    mutate(Norm1_foldchange= Norm1_transcripts/(Norm1_transcripts[name==foldChange0]))
} else {
  backsum<-backsum %>% 
    group_by(ps) %>% 
    mutate(Norm1_foldchange= Norm1_transcripts/1)
}


# backsum<-try(backsum %>% 
#   group_by(ps) %>% 
#   mutate(Norm2_foldchange= Norm2_transcripts/(Norm2_transcripts[name==foldChange0])), silent=T)

pdf(file = paste(new_path, title, "_", "_transcript_barplots.pdf", sep=""))
norm1_bar<-ggplot(backsum, aes(x=name, 
                               #color=temp, 
                               #fill=str,
                               y=as.numeric(transcripts)))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=as.numeric(transcripts)-transcript_error, 
                    ymax=as.numeric(transcripts)+transcript_error), color="black")+
  scale_y_log10()+
  facet_wrap(~ps, scales = "fixed")+
#  xlab("Strain(Temp)")+
  ylab(paste("Transcripts (log scale)"))+
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ggtitle(title)
print(norm1_bar)
dev.off()




pdf(file = paste(new_path, title, "_", norm_primer1, "_Normalized_transcript_barplots.pdf", sep=""))
norm1_bar<-ggplot(backsum, aes(x=name, y=as.numeric(Norm1_transcripts), color=as.character(temp), fill=str))+
  geom_bar(stat="identity")+
#  scale_y_log10()+
  facet_wrap(~ps, scales = "free")+
  xlab("Strain(Temp)")+
  ylab(paste("Transcripts (normalized to ", norm_primer1, ")"))+
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ggtitle(title)
print(norm1_bar)
dev.off()


pdf(file = paste(new_path, title, "_", norm_primer1,"_foldchange_barplot.pdf", sep=""))
ratios<-ggplot(backsum, aes(x=name, y=as.numeric(Norm1_foldchange), color=as.character(temp), fill=str))+
  geom_bar(stat="identity")+
  #  scale_y_log10()+
  facet_wrap(~ps, scales = "free")+
  xlab("Strain(Temp)")+
  ylab("Fold Change")+
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ggtitle(title)
print(ratios)
dev.off()



write.csv(backsum, paste(new_path,title, "_summary.csv", sep=""))
write.csv(summ, paste(new_path,title, "_rel_summary.csv", sep=""))



