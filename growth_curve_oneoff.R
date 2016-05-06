### Growth Curve Analysis
library(dplyr)
library(reshape2)
library(ggplot2)
data<-read.csv("../../Lab_Data/lucON20150416.csv")
# Remove blacks from analysis
data<-subset(data, strain !="Blank")

allLong<-melt(data, 
            varnames = grep("X.*", colnames(data), value=T),
            id.vars = c("read","col", "row", "rep", "strain"),
            value.name = "val", variable.name = "time")
allLong$time<-as.numeric(gsub("X(.*)","\\1", allLong$time))


all<-ggplot(subset(allLong, read=="od"), aes(x=time, y=val))+
  stat_summary(fun.y="mean", geom="line")+
  stat_summary(data=subset(allLong, read=="lum"), aes(y=val/10000), fun.y="mean", geom="line", color="red")+
  facet_grid(~strain)

all+xlab("Hours") +
  ylab("OD600 (Black) and Lum/10000(red)") +
  ggtitle("Overnight growth Dual-read OD600 and Luminescence")


