#   CONFIG for qRT_PCR_Analysis.R
library(reshape2)
library(XLConnect)               # load XLConnect package 


#################   Experiment title        #####################
experiment_title<-"sRNA_Demo"  #@@@@@@@@@@@@@@@@@@@@@@@@@@
experiment_title<-"PCR-PULLDOWN-3"  #@@@@@@@@@@@@@@@@@@@@@@@@@@
#experiment_date<-"201504"               #@@@@@@@@@@@@@@@@@@@@@@@@@@

title<- paste(gsub("-","", Sys.Date()),experiment_title, sep="_")

### primers for normalization

foldChange0<-"337(37)"
#################   Input and Output Directories        #####################

inputDirStd<-"../../Lab_data/RT_PCR/stdCurvesSYBR20150408/"   #@@@@@@@@@@@@@@@@@@@@@@@@@@
inputDirSamp<-"../../Lab_data/RT_PCR/rpoC_palC/"   #@@@@@@@@@@@@@@@@@@@@@@@@@@
inputDirSamp2<-"../../Lab_data/RT_PCR/srn1020_rsaD/"   #@@@@@@@@@@@@@@@@@@@@@@@@@@
#inputDirSamp3<-"../../Lab_data/RT_PCR/brnq/"   #@@@@@@@@@@@@@@@@@@@@@@@@@@
#inputDirSamp4<-"../../Lab_data/RT_PCR/brnq/"   #@@@@@@@@@@@@@@@@@@@@@@@@@@


####  make a new directory
new_path<-paste("~/GitHub/R/qRT-PCR/", title, "/",sep="")
dir.create(paste(gsub("(.*)/", "\\1", new_path)), showWarnings = TRUE, recursive = FALSE)

######## What did the nanodrop of your genomic DNA read? #########
inp_dna_concentration<-1.82*10^8

###########revamp dilution notation 
#  define input genomic dna concentration (copies per ul), and dilution series used
H<-inp_dna_concentration/1
G<-inp_dna_concentration/10
F<-inp_dna_concentration/100
E<-inp_dna_concentration/1000
D<-inp_dna_concentration/10000
C<-inp_dna_concentration/100000
B<-inp_dna_concentration/1000000
A<-inp_dna_concentration/10000000



#################   Define Primer Set Names #####################
ps1<-"rpoC" #@@@@@@@@@@@@@@@@@@@@@@@@@@   Normalization Primer
ps2<-"polC" #@@@@@@@@@@@@@@@@@@@@@@@@@@   
ps3<-"srn1020"     #@@@@@@@@@@@@@@@@@@@@@@@@@@
ps4<-"rsaD"     #@@@@@@@@@@@@@@@@@@@@@@@@@@
ps5<-""     #@@@@@@@@@@@@@@@@@@@@@@@@@@
ps6<-""     #@@@@@@@@@@@@@@@@@@@@@@@@@@

norm_primer1<-ps1
#try(norm_primer2<-ps2, silent = T) ##  DEPRECIATED


plist<-c(ps1,ps2, ps3, ps4, ps5, ps6)
pexist<-plist[grep("\\w", plist)]
#################   Standard Curve Inputs  #####################
filesStd<-grep("^\\~", list.files(inputDirStd, all.files = F), value = T, invert = T)  ##  this should get rid of ~$ temps
standard_cqwk<- loadWorkbook(paste(inputDirStd,  as.character(grep("Quantification Cq Results.xlsx", filesStd, value = T)), sep=""))
standard_cqdf<- readWorksheet(standard_cqwk, sheet="0")
standard_ampwk<- loadWorkbook(paste(inputDirStd,  as.character(grep("Quantification Amplification Results.xlsx", filesStd, value = T)), sep=""))
######  quick clean
standard_cqdf$row<-gsub("(\\D)(.*)", "\\1", standard_cqdf$Well)
standard_cqdf$col<-as.numeric(gsub("(\\D)(.*)", "\\2", standard_cqdf$Well))
standard_cqdf$well<-paste(standard_cqdf$row, standard_cqdf$col, sep="")
standard_cqdf<-standard_cqdf[complete.cases(standard_cqdf),]

# ps1df<-readWorksheet(standard_ampwk, sheet=ps1)
# ps1df$ps<-ps1
# 
# merge_maestro<-function(x){
#   t<-readWorksheet(standard_ampwk, sheet=x)
#   return(t)
# }
# #Import cq and amp xlsx files
# allall<-pexist%>% merge_maestro(.)
# allall<-allall%>% mutate(ps<-a

ps1df<-readWorksheet(standard_ampwk, sheet=ps1)
ps1df$ps<-ps1
ps1df<-melt(ps1df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well")
ps1df$rep_group<-gsub("(\\D)(.+)", "\\1", ps1df$well)

try(ps2df<-readWorksheet(standard_ampwk, sheet=ps2), silent = T)
try(ps2df$ps<-ps2, silent = T)
try(ps2df<-melt(ps2df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well"), silent = T)
try(ps2df$rep_group<-gsub("(\\D)(.+)", "\\1", ps2df$well), silent = T)

try(ps3df<-readWorksheet(standard_ampwk, sheet=ps3), silent = T)
try(ps3df$ps<-ps3, silent = T)
try(ps3df<-melt(ps3df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well"), silent = T)
try(ps3df$rep_group<-gsub("(\\D)(.+)", "\\1", ps3df$well), silent = T)

try(ps4df<-readWorksheet(standard_ampwk, sheet=ps4), silent = T)
try(ps4df$ps<-ps4, silent = T)
try(ps4df<-melt(ps4df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well"), silent = T)
try(ps4df$rep_group<-gsub("(\\D)(.+)", "\\1", ps4df$well), silent = T)

try(ps5df<-readWorksheet(standard_ampwk, sheet=ps5), silent = T)
try(ps5df$ps<-ps5, silent = T)
try(ps5df<-melt(ps5df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well"), silent = T)
try(ps5df$rep_group<-gsub("(\\D)(.+)", "\\1", ps5df$well), silent = T)

try(ps6df<-readWorksheet(standard_ampwk, sheet=ps6), silent = T)
try(ps6df$ps<-ps6, silent = T)
try(ps6df<-melt(ps6df, na.rm = T, value.name = "val", id.vars = c("Cycle", "ps"), variable.name = "well"), silent = T)
try(ps6df$rep_group<-gsub("(\\D)(.+)", "\\1", ps6df$well), silent = T)
#######################

prebind<- grep("ps.df", ls(), value=T) 
primerBind<-suppressWarnings(rbind(get(prebind[1]), 
                                   try(get(prebind[2]), silent =T ), 
                                   try(get(prebind[3]), silent =T ),
                                   try(get(prebind[4]), silent =T ),
                                   try(get(prebind[5]), silent =T ),
                                   try(get(prebind[6]), silent =T )))
primerBind<-primerBind[complete.cases(primerBind),]


###################  Sample plate Layout ####################



######## Configure inputs of your test conditions         #########
filesSamp<-grep("^\\~", list.files(inputDirSamp, all.files = F), value = T, invert = T)  ##  this should get rid of ~$ temps
sample_cqwk<- loadWorkbook(paste(inputDirSamp,  as.character(grep("Quantification Cq Results.xlsx", filesSamp, value = T)), sep=""))
cqdfSamp<- readWorksheet(sample_cqwk, sheet="0")
ampwk1<- loadWorkbook(paste(inputDirSamp,  as.character(grep("Quantification Amplification Results.xlsx", filesSamp, value = T)), sep=""))
wk1<-(readWorksheet(ampwk1, sheet="SYBR")) #sheet="SYBR"))
wk1t <- as.data.frame(t(wk1[,2:ncol(wk1)]))# Transpose 
wk1t$well <-row.names(wk1t)
wk1t$strain<-cqdfSamp$Sample
wk1t$condition<-cqdfSamp$Content
wk1t$ps<-cqdfSamp$Target
# 
try(filesSamp2<-grep("^\\~", list.files(inputDirSamp2, all.files = F), value = T, invert = T), silent=T) 
try(sample_cqwk2<- loadWorkbook(paste(inputDirSamp2,  as.character(grep("Quantification Cq Results.xlsx", filesSamp2, value = T)), sep="")), silent=T) 
try(cqdfSamp2<- readWorksheet(sample_cqwk2, sheet="0"), silent=T) 
try(ampwk2<- loadWorkbook(paste(inputDirSamp2,  as.character(grep("Quantification Amplification Results.xlsx", filesSamp2, value = T)), sep="")), silent=T) 
try(wk2<-(readWorksheet(ampwk2, sheet="SYBR")), silent=T)  #sheet="SYBR"))
try(wk2t <- as.data.frame(t(wk2[,2:ncol(wk2)])), silent=T) # Transpose 
try(wk2t$well <-row.names(wk2t), silent=T) 
try(wk2t$strain<-cqdfSamp2$Sample, silent=T) 
try(wk2t$condition<-cqdfSamp2$Content, silent=T) 
try(wk2t$ps<-cqdfSamp2$Target, silent=T) 
#
try(filesSamp3<-grep("^\\~", list.files(inputDirSamp3, all.files = F), value = T, invert = T), silent=T) 
try(sample_cqwk3<- loadWorkbook(paste(inputDirSamp3,  as.character(grep("Quantification Cq Results.xlsx", filesSamp3, value = T)), sep="")), silent=T) 
try(cqdfSamp3<- readWorksheet(sample_cqwk3, sheet="0"), silent=T) 
try(ampwk3<- loadWorkbook(paste(inputDirSamp3,  as.character(grep("Quantification Amplification Results.xlsx", filesSamp3, value = T)), sep="")), silent=T) 
try(wk3<-(readWorksheet(ampwk3, sheet="SYBR")), silent=T)  #sheet="SYBR"))
try(wk3t <- as.data.frame(t(wk3[,3:ncol(wk3)])), silent=T) # Transpose 
try(wk3t$well <-row.names(wk3t), silent=T) 
try(wk3t$strain<-cqdfSamp3$Sample, silent=T) 
try(wk3t$condition<-cqdfSamp3$Content, silent=T) 
try(wk3t$ps<-cqdfSamp3$Target, silent=T) 
#
try(filesSamp4<-grep("^\\~", list.files(inputDirSamp4, all.files = F), value = T, invert = T), silent=T) 
try(sample_cqwk4<- loadWorkbook(paste(inputDirSamp4,  as.character(grep("Quantification Cq Results.xlsx", filesSamp4, value = T)), sep="")), silent=T) 
try(cqdfSamp4<- readWorksheet(sample_cqwk4, sheet="0"), silent=T) 
try(ampwk4<- loadWorkbook(paste(inputDirSamp4,  as.character(grep("Quantification Amplification Results.xlsx", filesSamp4, value = T)), sep="")), silent=T) 
try(wk4<-(readWorksheet(ampwk4, sheet="SYBR")), silent=T)  #sheet="SYBR"))
try(wk4t <- as.data.frame(t(wk4[,4:ncol(wk4)])), silent=T) # Transpose 
try(wk4t$well <-row.names(wk4t), silent=T) 
try(wk4t$strain<-cqdfSamp4$Sample, silent=T) 
try(wk4t$condition<-cqdfSamp4$Content, silent=T) 
try(wk4t$ps<-cqdfSamp4$Target, silent=T) 
#
try(filesSamp5<-grep("^\\~", list.files(inputDirSamp5, all.files = F), value = T, invert = T), silent=T) 
try(sample_cqwk5<- loadWorkbook(paste(inputDirSamp5,  as.character(grep("Quantification Cq Results.xlsx", filesSamp5, value = T)), sep="")), silent=T) 
try(cqdfSamp5<- readWorksheet(sample_cqwk5, sheet="0"), silent=T) 
try(ampwk5<- loadWorkbook(paste(inputDirSamp5,  as.character(grep("Quantification Amplification Results.xlsx", filesSamp5, value = T)), sep="")), silent=T) 
try(wk5<-(readWorksheet(ampwk5, sheet="SYBR")), silent=T)  #sheet="SYBR"))
try(wk5t <- as.data.frame(t(wk5[,5:ncol(wk5)])), silent=T) # Transpose 
try(wk5t$well <-row.names(wk5t), silent=T) 
try(wk5t$strain<-cqdfSamp5$Sample, silent=T) 
try(wk5t$condition<-cqdfSamp5$Content, silent=T) 
try(wk5t$ps<-cqdfSamp5$Target, silent=T) 
#
try(filesSamp6<-grep("^\\~", list.files(inputDirSamp6, all.files = F), value = T, invert = T), silent=T) 
try(sample_cqwk6<- loadWorkbook(paste(inputDirSamp6,  as.character(grep("Quantification Cq Results.xlsx", filesSamp6, value = T)), sep="")), silent=T) 
try(cqdfSamp6<- readWorksheet(sample_cqwk6, sheet="0"), silent=T) 
try(ampwk6<- loadWorkbook(paste(inputDirSamp6,  as.character(grep("Quantification Amplification Results.xlsx", filesSamp6, value = T)), sep="")), silent=T) 
try(wk6<-(readWorksheet(ampwk6, sheet="SYBR")), silent=T)  #sheet="SYBR"))
try(wk6t <- as.data.frame(t(wk6[,6:ncol(wk6)])), silent=T) # Transpose 
try(wk6t$well <-row.names(wk6t), silent=T) 
try(wk6t$strain<-cqdfSamp6$Sample, silent=T) 
try(wk6t$condition<-cqdfSamp6$Content, silent=T) 
try(wk6t$ps<-cqdfSamp6$Target, silent=T) 

workbookList<-grep("wk\\dt", ls(), value=T)
all<-suppressWarnings(rbind(get(workbookList[1]), 
                            try(get(workbookList[2]), silent =T ), 
                            try(get(workbookList[3]), silent =T ),
                            try(get(workbookList[4]), silent =T ),
                            try(get(workbookList[5]), silent =T ),
                            try(get(workbookList[6]), silent =T )))
all[grep("Error.*", all$V1),]<-NA
all<-na.omit(all)
