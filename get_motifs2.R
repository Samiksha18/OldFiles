###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
library(seqinr)
input_path<-paste("~/GitHub/R/gbparse_output/N315.csv")

targetdf_gbparse<-read.csv(input_path, stringsAsFactors = F)
rot_targets_doc <-read.csv("rot_tagets.csv", stringsAsFactors = F)
reannotate<-read.csv("~/20160209_reannotate_output/uams1_with_CP000253_8325_dbanno.csv", stringsAsFactors = F)
protein<-"rot"
genome<-"N315"
genes_of_interest<-rot_targets_doc[rot_targets_doc$target==protein &
                                     rot_targets_doc$genome==genome, "old_locus_tag"]
homologues<-reannotate[reannotate$target_locus_tag %in% genes_of_interest, "locus_tag"]


subsetdf<- targetdf_gbparse[targetdf_gbparse$locus_tag %in% homologues,]
if (all(is.na(subsetdf$old_locus_tag))){subsetdf$old_locus_tag<-subsetdf$locus_tag}
subsetdf$full_seq<-""
subsetdf$full_seq<-subsetdf$prom_region


subsetdfnames<-as.list(paste(subsetdf$old_locus_tag,"|",subsetdf$locus_tag, "|",subsetdf$name, "|",subsetdf$direction, sep=""))
subsetdfseq<-as.list(subsetdf$full_seq)

output_name<- gsub("(.*)(\\.csv)","\\1", gsub("(.*\\/)(.*\\.csv)","\\2", input_path))
output_path<-paste("~/Desktop/MEME/","MEME_", gsub("-|:|\\D", "", Sys.time()),"_",output_name,"_",protein, "_input.fasta", sep="")
write.fasta(sequences=subsetdfseq, 
            names=subsetdfnames, nbchar = 60, 
            file.out= output_path, open = "w")

###############################################################################
###############################################################################
############################### CodY ##########################################
###############################################################################
###############################################################################
library(seqinr)
input_path<-paste("~/GitHub/R/gbparse_output/uams1.csv")
protein<-"codY"

targetdf_gbparse<-read.csv(input_path, stringsAsFactors = F)
reannotate<-read.csv("../uams1rnaseq/20151217_Kegg_UAMS-1_RNA-seq/20151217_Kegg_UAMS-1_RNA-seq_foldChangesWTvsNull.csv", stringsAsFactors = F)
threshold<-5
thresholds<-c(1/threshold,threshold)
#genome<-"CP000253_8325"
genes_of_interest<-reannotate[reannotate$fc<=thresholds[1] |
                                reannotate$fc>=thresholds[2], "qv_anno"]
genes_of_interest<-gsub("(.{10}) .*", "\\1", genes_of_interest)
min_length<-50

subsetdf<- targetdf_gbparse[targetdf_gbparse$locus_tag %in% genes_of_interest &
                              nchar(targetdf_gbparse$preceeding_intergenic_region)>min_length,]
if (all(is.na(subsetdf$old_locus_tag))){subsetdf$old_locus_tag<-subsetdf$locus_tag}
subsetdf$full_seq<-""
subsetdf$full_seq<-subsetdf$preceeding_intergenic_region


subsetdfnames<-as.list(paste(subsetdf$old_locus_tag,"|",subsetdf$locus_tag, "|",subsetdf$name, "|",subsetdf$direction, sep=""))
subsetdfseq<-as.list(subsetdf$full_seq)

output_name<- gsub("(.*)(\\.csv)","\\1", gsub("(.*\\/)(.*\\.csv)","\\2", input_path))
output_path<-paste("~/Desktop/MEME/","MEME_", gsub("-|:|\\D", "", Sys.time()),"_",output_name,"_",protein, "_input.fasta", sep="")
write.fasta(sequences=subsetdfseq, 
            names=subsetdfnames, nbchar = 60, 
            file.out= output_path, open = "w")

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
library(seqinr)
input_path<-paste("~/GitHub/R/gbparse_output/N315.csv")
protein<-"rot"

targetdf_gbparse<-read.csv(input_path, stringsAsFactors = F)
reannotate<-read.csv("rot_tagets.csv", stringsAsFactors = F)
threshold<-7
#thresholds<-c(1/threshold,threshold)
#genome<-"CP000253_8325"
genes_of_interest<-reannotate[reannotate$target==protein &
                                reannotate$genome=="N315" &
                                reannotate$fold_change>threshold, "old_locus_tag"]
#genes_of_interest<-gsub("(.{10}) .*", "\\1", genes_of_interest)
min_length<-40

subsetdf<- targetdf_gbparse[targetdf_gbparse$old_locus_tag %in% genes_of_interest &
                              nchar(targetdf_gbparse$preceeding_intergenic_region)>min_length,]
if (all(is.na(subsetdf$old_locus_tag))){subsetdf$old_locus_tag<-subsetdf$locus_tag}
subsetdf$full_seq<-""
subsetdf$full_seq<-subsetdf$preceeding_intergenic_region


subsetdfnames<-as.list(paste(subsetdf$old_locus_tag,"|",subsetdf$locus_tag, "|",subsetdf$name, "|",subsetdf$direction, sep=""))
subsetdfseq<-as.list(subsetdf$full_seq)

output_name<- gsub("(.*)(\\.csv)","\\1", gsub("(.*\\/)(.*\\.csv)","\\2", input_path))
output_path<-paste("~/Desktop/MEME/","MEME_", gsub("-|:|\\D", "", Sys.time()),"_",output_name,"_",protein, "_input.fasta", sep="")
write.fasta(sequences=subsetdfseq, 
            names=subsetdfnames, nbchar = 60, 
            file.out= output_path, open = "w")

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################