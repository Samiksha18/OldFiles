#!/usr/bin/Rscript
####       Growth Curve Analysis        ####
####          Nick Waters 2015          ####
####              20150421              ####
# install.packages("XLConnect")
### 
argv<-commandArgs()

library(XLConnect)               # load XLConnect package 
library(ggplot2)
library(dplyr)
library(reshape2)
library(grid)

print(argv)
print(paste("length of argv:", length(argv)))
conff<-grep("CONFIG", argv, value=T)
source(conff)
# source("on0419CONFIG.R")
# source("201505089lucCONFIG.R" )
# source("tboxCONFIG.R")
source("growthcurve0820CONFIG.R")

# Remove blanks from analysis
#data<-subset(data, strain !=blank)
data<-subset(data, cond !=blank)
data<-subset(data, cond !="null")
media="TSB"
if(grepl("TSB", media)){
  data$rowMeans<-rowMeans(data[,grep("(X\\d.*)", colnames(data))])
  data$growth<-ifelse(data$rowMeans<=mediaMin, "no_growth", "growth")
  data<-data%>%group_by( col, row)%>% mutate(rm=ifelse(any(growth=="no_growth")&strain!=blank, "remove","keep"))
  print(paste("the following rows showed no perceptible change in", read_1, ":"))
  print(unique(data[data$rm=="remove", c("col","row","strain" )]))
  data<-data[data$rm=="keep",]
  print("they have been excluded")
  data$rm<-NULL
  data$growth<-NULL
  data$rowMeans<-NULL
}

allLong<-melt(data, 
              varnames = grep("X.*", colnames(data), value=T),
              id.vars = c(note, read, col, cond, row, rep, strain),
              value.name = "rawVal", variable.name = "time")
allLong$time<-as.numeric(gsub("(X)(\\d*)(\\D*)","\\2", allLong$time))
allLong$rawVal<-as.numeric(allLong$rawVal)

if(normMethod=="none"){allLong$val<-allLong$rawVal}
#Col
if(normMethod=="col"){
  allLong$key <-paste(allLong$time,allLong$read, allLong$col, sep="_")
  blanksdf<-allLong[allLong$strain==blank,c("key","rawVal") ]
  toMerge<-allLong
  allLong<-merge(toMerge, blanksdf, by="key", all.x=T)
  allLong$val<-ifelse((allLong$rawVal.x-allLong$rawVal.y)<=0,0, allLong$rawVal.x-allLong$rawVal.y)
}
#Row ----  make this work
if(normMethod=="row"){
  allLong$key <-paste(allLong$time,allLong$read, allLong$row, allLong$rep, sep="_")# sep="_"
  blanksdf<-allLong[allLong$strain==blank,c("key","rawVal") ]
  toMerge<-allLong
  allLong<-merge(toMerge, blanksdf, by="key", all.x=T)
  allLong$val<-allLong$rawVal.x-allLong$rawVal.y
}
#Mean by read type
#check
if(normMethod=="rawMean"){
  ifelse(length(table(allLong$read==read_1))==2,"Both read types accounted for!", "You have a problem, sir!")
  allLong$blankMean<-ifelse(allLong$read==read_1, mean(allLong[allLong$read==read_1,"rawVal"]), mean(allLong[allLong$read==read_2,"rawVal"]))
  allLong$val<-allLong$rawVal-allLong$blankMean
}  
###############  PLOTTING
allLong<-allLong[allLong$strain!=blank,]
all<-ggplot(subset(allLong, read==read_1), aes(x=time, y=val))+scale_y_log10()+
  stat_summary(fun.y="mean", geom="line")+
  facet_wrap(strain~cond)+
  xlab(time_unit)+
  ylab(paste(read_1," (Black)"))+
  ggtitle(paste(title))
####  dualreading?
if(dualRead==T){
  all<-ggplot(subset(allLong, read==read_1), aes(x=time, y=val))+scale_y_log10()+
    stat_summary(fun.y="mean", geom="line")+
    stat_summary(data=subset(allLong, read==read_2), aes(y=val/scale),
                   fun.y="mean", geom="line", color="red")+
    facet_wrap(strain~cond)+
    xlab(time_unit) +
    ylab(paste(read_1," (Black) and ", read_2,"/", scale, " (red)", sep="")) +
    ggtitle(paste(title, read_1,"vs.", read_2))}
all

      



p1<-ggplot(subset(allLong, read==read_1, color="black"), aes(x=time, y=val))+
  scale_y_log10()+
  stat_summary(fun.y="mean", geom="line")+facet_wrap(strain~cond)
  ####  dualreading?
p2<-ggplot(subset(allLong, read==read_2), aes(x=time, y=val))+ 
  scale_y_log10()+
  stat_summary(fun.y="mean", geom="line", color="red")+facet_wrap(strain~cond)
grid.newpage()
grid.draw(rbind(ggplotGrob(p1),ggplotGrob(p2), size="last"))

