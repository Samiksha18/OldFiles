####       Processing qRT-PCR data      ####
####          Nick Waters 2015          ####
####              20150416              ####
# library(gdata)
# install.packages("XLConnect")
library(XLConnect)               # load XLConnect package 
library(ggplot2)
library(dplyr)
library(reshape2)
#define primer sets
ps1<-"rpoC"
ps2<-"polC"
ps3<-"srn1020"
ps4<-"rsaD"

#Import cq and amp xlsx files
cqwk<- loadWorkbook("../../Lab_data/RT_PCR/stdCurvesSYBR20150408/standardCurvesRpocPalc1020Rsad -  Quantification Cq Results.xlsx") 
cqdf<- readWorksheet(cqwk, sheet="0")
ampwk<- loadWorkbook("../../Lab_data/RT_PCR/stdCurvesSYBR20150408//standardCurvesRpocPalc1020Rsad -  Quantification Amplification Results.xlsx")

ps1df<-readWorksheet(ampwk, sheet=ps1)
ps1df$ps<-ps1
ps1df<-melt(ps1df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well")
ps1df$rep_group<-gsub("(\\D)(.+)", "\\1", ps1df$well)

ps2df<-readWorksheet(ampwk, sheet="palC")
ps2df$ps<-ps2
ps2df<-melt(ps2df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well")
ps2df$rep_group<-gsub("(\\D)(.+)", "\\1", ps2df$well)

ps3df<-readWorksheet(ampwk, sheet=ps3)
ps3df$ps<-ps3
ps3df<-melt(ps3df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well")
ps3df$rep_group<-gsub("(\\D)(.+)", "\\1", ps3df$well)

ps4df<-readWorksheet(ampwk, sheet=ps4)
ps4df$ps<-ps4
ps4df<-melt(ps4df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well")
ps4df$rep_group<-gsub("(\\D)(.+)", "\\1", ps4df$well)

df<-rbind(ps1df, ps2df, ps3df, ps4df)

d<-ggplot(df, aes(x=Cycle, y=val, group=well))+geom_line()+facet_wrap(~ps)
d
##### Approximate Cq value based on results for H1=11.10
sample_well="H1"
sample_Cq=11.10
subset<-subset(df, well==sample_well)
testInt<-approx(y = subset$val, x = subset$Cycle, xout = sample_Cq)
print(paste("the sample well", sample_well, "has a Cq of", sample_Cq, ", therefore the threshhold is", testInt$y))
d+geom_hline(y=testInt$y)

##  Calulate Cq and Cq_mean values
df<-df %>%group_by(well) %>% mutate(Cq= approx(y=Cycle, x=val, xout=testInt$y)$y)
df<-df %>%group_by(rep_group,ps)%>% mutate(Cq_mean= mean(Cq))
df<-df %>%group_by(rep_group,ps)%>% mutate(Cq_StDev= sd(Cq))

# create summarized dataframe
sumdf<-as.data.frame(cbind(df$ps, df$rep_group, df$Cq_mean, df$Cq_StDev))
colnames(sumdf)<-c("ps", "rep_group", "Cq_mean", "Cq_StDev")
sumdf<-sumdf[!duplicated(sumdf), ]

###########revamp dilution notation 
#  define input genomic dna concentration (copies per ul), and dilution series used
inp_dna_concentration<-1.82*10^8
H<-inp_dna_concentration/1
G<-inp_dna_concentration/10
F<-inp_dna_concentration/100
E<-inp_dna_concentration/1000
D<-inp_dna_concentration/10000
C<-inp_dna_concentration/100000
B<-inp_dna_concentration/1000000
A<-inp_dna_concentration/10000000
#wells<-(c(A, B, C, D, E, F, G, H))

#  Make a column for true concentration
sumdf$Cq_mean<-as.numeric(sumdf$Cq_mean)
sumdf$Cq_StDev<-as.numeric(sumdf$Cq_StDev)
sumdf<-sumdf%>%group_by(rep_group,ps) %>% mutate(conc=get(as.character(rep_group)))
write.csv(sumdf, "synth_Cq.csv")
sumdf$logConc<-log(sumdf$conc)

coefps1<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps1)))
coefps2<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps2)))
coefps3<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps3)))
coefps4<-coef(lm(Cq_mean~log(conc, 10), data = subset(sumdf, ps==ps4)))

coefdf<-bind_cols(as.data.frame(coefps1),as.data.frame(coefps2),as.data.frame(coefps3),as.data.frame(coefps4))
coefdfr<-c(ps1, ps2, ps3, ps4)
colnames(coefdf)<-coefdfr
coefdf
write.csv(coefdf, "coefdf.csv")


coefps4

###  Double check-  does green match black????????
plot_lr<-ggplot(sumdf, aes(y=Cq_mean, x=log(conc, 10)))+
#  scale_x_log10()+
  geom_point()+
  geom_smooth(aes(group=ps), method="lm", se=FALSE, color="green")+  #True
  geom_abline(intercept =coefps1[1], slope = coefps1[2]) # Does black match rpoC Green?
  
plots<-sumdf %>% group_by(ps)%>% do(plots=plot_lr %+% . + facet_wrap(~ps)+geom_abline(intercept =coefdf$ps[1], slope = coefps1[2]))

  

