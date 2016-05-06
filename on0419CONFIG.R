#   CONFIG for growth curve Analysis
#################   Experiment title        #####################
experiment_title<-"luc_ON"
experiment_date<-"20150416"
title<-paste(experiment_date, "_",experiment_title, sep="")

####  Input from CSV
#data<-read.csv("../../Lab_Data/lucON20150416.csv")

####  Input from Excel
excel_data<- loadWorkbook("../../Lab_data/Maggie//on_20150417.xlsx")
data<-readWorksheet(excel_data, sheet="Magellan Sheet 1")

#####  Optional Subset    @@@@@@@@@@@
data<-subset(data, !data$strain %in% grep("gua", data$strain, value = T))
###### Define Blank
blank<-  "NA"    ##  What is a blank well called?

#####   Nomenclature for column headings
read<-   "reading"     ##  Whats the col name for what readings are being taken?
col<-    "col"      ##  Column of plate (X)
row<-    "row"      ##  Row of plate    (Y)
cond<-    "cond"      ##  Use this to group
rep<-     "rep"     ##  unique per sample
strain<- "strain"   ##  What is delineating samples ?
#####  Sanatize as needed
colnames(data)<-gsub(read,"read", colnames(data))
colnames(data)<-gsub(col,"col", colnames(data))
colnames(data)<-gsub(read,"row", colnames(data))
colnames(data)<-gsub(cond,"cond", colnames(data))
##  Missing Column?
print(colnames(data))
data<-data%>%mutate(rep="1")
print(colnames(data))


####
read_1<-"od"
read_2<-"lum"
time_unit<-"Hours"