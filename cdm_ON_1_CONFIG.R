#   CONFIG for growth curve Analysis
library(XLConnect)
library(dplyr)
#################   Experiment title        #####################
experiment_title<-"CDM+-ilv"
experiment_date<-"20150521"
title<-paste(experiment_date, "_",experiment_title, sep="")

####  Input from CSV
#data<-read.csv("../../Lab_Data/csvpath.csv")

####  Input from Excel
excel_data<- loadWorkbook("../../Lab_data/Maggie/2105LuxONCDM.xlsx")
data<-readWorksheet(excel_data, sheet="Magellan Sheet 1")

#####  Optional Subset    @@@@@@@@@@@
#data<-subset(data, !data$strain %in% grep("gua", data$strain, value = T))
###### Define Blank
blank<-  "blank"    ##  What is a blank well called?

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
##  Missing ID Column?
ids<-grep("X\\d",colnames(data), invert = T, value = T)
print(paste("ID columns: ",unlist(ids)))
print("Everything look ok? Good!")
#data<-data%>%mutate(rep="1")
#print(colnames(data))
####
read_1<-"od"
read_2<-"lum"
time_unit<-"Hours"
scale<-10000    #######  reduce data? if no, than this is 1

