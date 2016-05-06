#  make fasta from gbparseR.r
library(seqinr)
input_path<-"~/GitHub/R/gbparse_output/N315.csv"
a<-read.csv(input_path)
a<-a[a$type=="CDS",]
names<-as.list(paste(a$old_locus_tag,"|",a$locus_tag, "|",a$name, "|",a$direction, sep=""))
seqs <-as.list(a$prom_region)

write.fasta(sequences=seqs, 
            names=names, nbchar = 60, 
            file.out= gsub("\\.csv","_upstream\\.fasta",input_path), open = "w")
