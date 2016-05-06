#   CONFIG for growth curve Analysis
#################   Experiment title        #####################
library(XLConnect)
library(dplyr)
experiment_title<-"luc_ON"
experiment_date<-"20150508"
title<-paste(experiment_date, "_",experiment_title, sep="")

####  Input from CSV
#data<-read.csv("../../Lab_Data/050815.xlsx")

####  Input from Excel
excel_data<- loadWorkbook("../../Lab_data/050815.xlsx")
data<-readWorksheet(excel_data, sheet="Magellan Sheet 1")

#####  Optional Subset    @@@@@@@@@@@
#data<-subset(data, !data$strain %in% grep("gua", data$strain, value = T))

###### Define Blank
blank<-  "x"    ##  What is a blank well called?

##### What do we scale the lum plots by to get them on the same plots?
scale<-5000

#####   Nomenclature for column headings
note<-   "note"     ##  for experiment notes
read<-   "read"     ##  Whats the col name for what readings are being taken?
col<-    "col"      ##  Column of plate (X)
row<-    "row"      ##  Row of plate    (Y)
cond<-    "cond"      ##  Use this to group
rep<-     "rep"     ##  unique per sample
strain<- "strain"   ##  What is delineating samples ?

##  Missing Column?
#print(colnames(data))
#data<-data%>%mutate(cond="1")
#print(colnames(data))

#####  Sanatize as needed
colnames(data)<-gsub(read,"read", colnames(data))
colnames(data)<-gsub(col,"col", colnames(data))
colnames(data)<-gsub(row,"row", colnames(data))
colnames(data)<-gsub(cond,"cond", colnames(data))
###
read_1<-"od"
read_2<-"lum"
time_unit<-"Hours"

