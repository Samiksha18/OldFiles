####   Analyzing qrt-PCR data ####
library(XLConnect)               # load XLConnect package 
library(ggplot2)
library(dplyr)
library(reshape2)
#define primer sets
ps1<-"rpoC"
ps2<-"polC"
ps3<-"srn1020"
ps4<-"rsaD"
### Samples
s1<-"337(37)a"  #row A
s2<-"372(37)a"  #row B
s3<-"337(37)b"  #row C
s4<-"372(37)b"  #row D
s5<-"337(25)"   #row E
s6<-"372(25)"   #row F

###  Load Summary
cqwk<- loadWorkbook("../../Lab_data/RT_PCR/rpoC_palC/rpoCpalC_samples -  Quantification Cq Results.xlsx")
cqdf<- readWorksheet(cqwk, sheet="0")
#    Load in one or more curve files
ampwk1<- loadWorkbook("../../Lab_data/RT_PCR/rpoC_palC/rpoCpalC_samples -  Quantification Amplification Results_monkey.xlsx")
wk1<-(readWorksheet(ampwk1, sheet="Wide")) #sheet="SYBR"))
ampwk2<- loadWorkbook("../../Lab_data/RT_PCR/srn1020_rsaD/rsaD_1020_samples -  Quantification Amplification Results.xlsx")
wk2<-(readWorksheet(ampwk2, sheet="Wide")) #sheet="SYBR"))

all<-bind_rows(wk1, wk2)



###  Make data long, remove nonsense
all<-melt(all, id.vars = c("well", "strain", "ps", "condition"), value.name = "val")
all<-all[complete.cases(all),]
all<-all %>% group_by(row.names(all)) %>% mutate(name=get(as.character(strain)))
all$cycle<-as.numeric(gsub("(\\D)(.+)", "\\2", all$variable))
all$temp<-as.numeric(gsub("(.*)\\((.*)\\)(.*)", "\\2", all$name))
all<-all[all$condition != "N",]
all$str<-gsub("(.*)\\((.+)", "\\1", all$name)

####  Nice quick visual check
p<-ggplot(all, aes(x=cycle, y=val, group=well, color=name))+
  geom_line()+
  facet_grid(condition~ps)
p
##Best Guess
thresh<- 408.7006
#
all<-all %>%group_by(well, ps) %>% mutate(Cq= approx(y=cycle, x=val, xout=thresh)$y)
all<-all %>%group_by(strain,ps,condition)%>% mutate(Cq_mean= mean(Cq))
all<-all %>%group_by(strain,ps,condition)%>% mutate(Cq_StDev= sd(Cq))

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
summ<-subset(summ, condition !="NRT")

####  Sanity Check
duplicated(x = summ$Cq_mean)


coefdf<-read.csv("coefdf.csv")
ps1sumdf<-subset(summ, ps==ps1)
ps2sumdf<-subset(summ, ps==ps2)
ps3sumdf<-subset(summ, ps=="sra1020")
ps4sumdf<-subset(summ, ps==ps4)


ps1sumdf$transcripts<-10^((ps1sumdf$Cq_mean-coefdf[,ps1][1])/(coefdf[,ps1][2]))
ps2sumdf$transcripts<-10^((ps2sumdf$Cq_mean-coefdf[,ps2][1])/(coefdf[,ps2][2]))
ps3sumdf$transcripts<-10^((ps3sumdf$Cq_mean-coefdf[,ps3][1])/(coefdf[,ps3][2]))
ps4sumdf$transcripts<-10^((ps4sumdf$Cq_mean-coefdf[,ps4][1])/(coefdf[,ps4][2]))

plot(y=ps4sumdf$transcripts, x=ps4sumdf$name)


backsum<-bind_rows(ps1sumdf,ps2sumdf, ps3sumdf, ps4sumdf)
backsum$str<-gsub("(.*)\\((.+)", "\\1", backsum$name)

ggplot(backsum, aes(x=name, y=as.numeric(transcripts), color=temp, fill=str))+
  geom_bar(stat="identity")+
  scale_y_log10()+
  facet_wrap(~ps, scales = "fixed")+
  xlab("Strain(Temp)")+
  ylab("transcripts")+
  ggtitle("sRNA qRT-PCR")
  

  



