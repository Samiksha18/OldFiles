#!/usr/bin/Rscript

"
usage:
Rscript uams1_preseq.R counts_file, genome, manual_annotations

outputs a counts file with everything attached
"

args<-commandArgs(T)
# args<-c("../uams1rnaseq/counts_files/All_counts_Brinsmade2.txt",
#         "../../Py/uams1_with_sax_kegg_dbanno.csv", 
#         "../../Lab_data/alt_manual_annotation.xlsx"
#         )

# read in arguments  MAKE CHECKS HERE

counts_handle<-args[1]  
genome_handle<-args[2]
man_annno_bool<-F

if (length(args)>2){
  man_annno_bool<-T
  print(paste("reading in manual annotations from", args[3]))
  man_anno_handle<-args[3]
}


##############################################
input_genome<-read.csv(genome_handle, header=T,sep=",")
anno_column<-grep("db_anno", colnames(input_genome))
#input_genome_alt_ID<-read.csv(genome_handle_alt_ID, header=T,sep=",")
input_counts <-read.csv(counts_handle, header=T,sep="\t")
idRow<-1 # for counts

##############################################################################
# Get rid of NA cols
noNA<-input_counts[, colSums(is.na(input_counts)) != nrow(input_counts)]
print(paste("these are the Cols without NA's:", (colnames(noNA))))
print("got a problem with that? clean up your data, or write a program to handle your mess")
################################################################################
#        Get rid Rows with text or NA's
# Make dataframe with counts, make all characters
b<-as.data.frame(apply(noNA[,(idRow+1):length(colnames(noNA))], 2, as.character))
# from that,make all integers.  This forces non-integers to NA
c<-as.data.frame(apply(b, 2, as.numeric))
#Put back the row names before you do something foolish
row.names(c)<-input_counts[,idRow]
# Get rid of row will ALL NA's  (thanks to Wookai @SO)
noChar<- c[rowSums(is.na(c))!=ncol(c), ]

#noChar<-na.omit(c)
print(paste("removed", nrow(b)-nrow(noChar), "row(s) for containing characters"))
##  Still got any NA's lurking in there?  Clobber them with complete cases, and survey the damage
noBadRow<-noChar[complete.cases(noChar),] 
test1<-ifelse(nrow(noChar)-nrow(noBadRow)>0, "you got issues still", "looking pretty fine, there!")
print("Alright, say 'AWWWWWW', hmm , lets see now....")
print(test1)
#write.csv(noChar,'counts_files/Brinsmade2.csv')
################################################################################
#                         SPECIFIC TO QV15 DATASET
#    ID row contains IDs in triplicate separated by underscore, 
#        and there is an underscore in the ID, like QV15_00056_QV15_00056_QV15_00056
#       and there are also the wild genes like the first gryB and tuf 
#               so its like tuf_tuf_tuf
#                      oh joy.... 

#       release the REGEX!!!!!!
################################################################################

cleanCounts<-noChar
m<-gsub("(.*_)(.*_).*", "\\1", row.names(cleanCounts))
n<-gsub("(.*)_", "\\1", m)
o<-gsub("(.{10})(_.*)", "\\1", n)

row.names(cleanCounts)<-o

#shazam.  


################################################################################
# Write out simplekey with locus tags from both counts and genome
print("writing out key to unify disparate locus tags in the genome and counts file")
simpleKey<-data.frame(row.names(cleanCounts), input_genome$locus_tag )
colnames(simpleKey)<-c("counts", "genome")
print(paste("writing Key to ", getwd(), "and heres a preview:"))
print(head(simpleKey))
write.csv(simpleKey, paste(getwd(),"/", gsub("(.*\\/)(.*?)(\\.csv)","\\2", genome_handle), "_simpleKey.csv", sep=""))

################################################################################
input_counts<-cleanCounts
################################################################################
if (!all( row.names(input_counts) == input_genome$locus_tag)){
  
#####  Replace Gene code with informative Gene name
input_counts$locus_tag<-row.names(input_counts)
input_countsUnif<-merge(simpleKey,input_counts, by.y ="locus_tag", by.x ="counts", sort = F )
input_countsUnif$oldschool<-paste(input_countsUnif$genome, "|",input_countsUnif$counts, sep="")
input_counts$both_locus_tags<-input_countsUnif$oldschool
input_counts$locus_tag<-gsub("(^.+)\\|(.*)", "\\1", input_counts$both_locus_tags)
}



##### Add manual annotation 
combined<-merge(input_counts, input_genome, by = "locus_tag", all.y=T)

if (man_annno_bool){
  if(grepl("\\.xlsx", args[3])){
    suppressMessages(library(XLConnect))
    man_anno_data<- loadWorkbook(man_anno_handle)
    manualAnnotation<-readWorksheet(man_anno_data, sheet="Sheet1")
  } else if (grepl("\\.csv", args[3])){
    manualAnnotation<-read.csv(args[3], header=T,sep=",")
  } else{
    stop("your manual input must be either a excel file on a tab called 'Sheet1' or a csv file ")
  }
  manualAnnotation<-manualAnnotation[!duplicated(manualAnnotation$locus_tag), ]
  manualAnnotation<-manualAnnotation[complete.cases(manualAnnotation$locus_tag),]
  #virGenes$vir<-"V"
  combined<-merge(combined, manualAnnotation, by = "locus_tag", all.x=T)
  
}
################################################################################
####################  Prepare incoming datasets  ################################
#  for using the uniprot mapping:
col_to_extract<-grep("db_anno", colnames(combined), value=T) #  Raw description col
#col_to_extract<-"uniprot300_anno"
desc1<-paste("altAnnotation_path", sep="")               # secondary description or accession 
#desc1<-"usa300_uni"
desc2<-paste("altAnnotation_locus_tag", sep="")              # index name
#desc2<-"usa300gene"
desc3<-paste("altAnnotation_desc", sep="")              # primary description name
#desc3<-"usa300_desc"
# desc4<-"AltID"                     # alternate genome id
# desc5<-"Other"                     # alternate genome desc

combined[desc1]<-gsub("(.+?)\\|(.+)\\|(.+)\\|(.+)", "\\4", combined[,col_to_extract])
combined[desc2]<-gsub("(.+?)\\|(.+)\\|(.+)\\|(.+)", "\\1", combined[,col_to_extract])
combined[desc3]<-gsub("(\\(RefSeq\\)\\s.*?)","", gsub("(.+?)\\|(.+)\\|(.+)\\|(.+)", "\\3", combined[,col_to_extract]))
#combined[desc4]<-gsub("(.+?)\\|(.+)\\|(.+)\\|(.+)", "\\1",combined$altID)
combined$note<-as.character(combined$note)

# write out combinedd data  NOW or later?
write.csv(combined, paste(getwd(),"/",gsub("(.*)\\/(.*)\\..*","\\2", args[1]),"_groomed", sep=''))

