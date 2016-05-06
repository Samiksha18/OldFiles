#   CONFIG for qRT_PCR_Analysis.R
#################   Experiment title        #####################
experiment_title<-"sRNA_ColdShock"
experiment_date<-"20150417"

title<-paste(experiment_date, "_",experiment_title, sep="")
####  make a new directory
new_path<-paste("~/GitHub/R/qRT-PCR/", title, "/",sep="")
dir.create(paste(gsub("(.*)/", "\\1", new_path)), showWarnings = TRUE, recursive = FALSE)



#################   Define Primer Set Names #####################
ps1<-"rpoC"
ps2<-"polC"
ps3<-"srn1020"
ps4<-"rsaD"
ps5<-" "

primer_list<-grep("\\S", c(ps1,ps2, ps3, ps4, ps5), value = T)
#################   Standard Curve Inputs  #####################

standard_cqwk<- loadWorkbook("../../Lab_data/RT_PCR/stdCurvesSYBR20150408/standardCurvesRpocPalc1020Rsad -  Quantification Cq Results.xlsx") 
standard_cqdf<- readWorksheet(standard_cqwk, sheet="0")
standard_ampwk<- loadWorkbook("../../Lab_data/RT_PCR/stdCurvesSYBR20150408//standardCurvesRpocPalc1020Rsad -  Quantification Amplification Results.xlsx")
######  quick clean
standard_cqdf$row<-gsub("(\\D)(.*)", "\\1", standard_cqdf$Well)
standard_cqdf$col<-as.numeric(gsub("(\\D)(.*)", "\\2", standard_cqdf$Well))
standard_cqdf$well<-paste(standard_cqdf$row, standard_cqdf$col, sep="")


######## What did the nanodrop of your genomic DNA read? #########
inp_dna_concentration<-1.82*10^8

######## Configure inputs of your test conditions         #########

cqwk<- loadWorkbook("../../Lab_data/RT_PCR/rpoC_palC/rpoCpalC_samples -  Quantification Cq Results.xlsx")
cqdf<- readWorksheet(cqwk, sheet="0")
#    Load in one or more curve files
ampwk1<- loadWorkbook("../../Lab_data/RT_PCR/rpoC_palC/rpoCpalC_samples -  Quantification Amplification Results_monkey.xlsx")
wk1<-(readWorksheet(ampwk1, sheet="Wide")) #sheet="SYBR"))
ampwk2<- loadWorkbook("../../Lab_data/RT_PCR/srn1020_rsaD/rsaD_1020_samples -  Quantification Amplification Results.xlsx")
wk2<-(readWorksheet(ampwk2, sheet="Wide")) #sheet="SYBR"))

### Samples
s1<-"337(37)a"  #row A
s2<-"372(37)a"  #row B
s3<-"337(37)b"  #row C
s4<-"372(37)b"  #row D
s5<-"337(25)"   #row E
s6<-"372(25)"   #row F


