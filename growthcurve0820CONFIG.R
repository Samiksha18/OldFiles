#   CONFIG for growth curve Analysis
suppressMessages(library(XLConnect))
library(dplyr)
#################   Experiment title        #####################
experiment_title<-"rsaD_codYrsaD_50%TSB_dualRead"
experiment_date<-"20150819"
title<-paste(experiment_date, "_",experiment_title, sep="")

####  Input from CSV
#data<-read.csv("../../Lab_Data/csvpath.csv")

####  Input from Excel
excel_data<- loadWorkbook("../../Lab_data/Maggie/ON730_731_337_20150820.xlsx")
data<-readWorksheet(excel_data, sheet="Magellan Sheet 1")

#####  Optional Subset    @@@@@@@@@@@
#data<-subset(data, !data$strain %in% grep("gua", data$strain, value = T))
###### Define Blank
blank<-  "blank"    ##  What is a blank well called?
normMethod<-"col"# "rawMean" row, , none, or col
mediaMin<-0.05
media<-"TSB"
#####   Nomenclature for column headings
note<-   "note"     ##  Where are the notes?
read<-   "read"     ##  Whats the col name for what readings are being taken?
col<-    "col"      ##  Column of plate (X)
row<-    "row"      ##  Row of plate    (Y)
cond<-    "cond"      ##  Use this to group
rep<-     "rep"     ##  unique per sample
strain<- "strain"   ##  What is delineating samples ?

#####  Sanitize as needed
colnames(data)<-gsub(read,"read", colnames(data))
colnames(data)<-gsub(col,"col", colnames(data))
colnames(data)<-gsub(row,"row", colnames(data))
colnames(data)<-gsub(cond,"cond", colnames(data))
if(length(grep("cond", colnames(data)))==0){
  data$cond<-"standard"}
####    Incidental filtering
data<-data[!data$strain=="WT",]
####

##  Missing ID Column?
ids<-grep("X\\d",colnames(data), invert = T, value = T)
print("ID columns:")
print(unlist(ids))
print("Everything look ok? Good!")
#data<-data%>%mutate(rep="1")
#print(colnames(data))
####
dualRead=T
read_1<-"OD600"
read_2<-"dsRed"
time_unit<-"Seconds"
scale<-10000    #######  reduce data? if no, than this is 1

