#!/usr/bin/Rscript
library(shiny)
library(Biostrings)
setwd("~/GitHub/R/shinyApps/shinyAppV6/")

genome<-read.csv("genome_flatfile.csv", stringsAsFactors = F)
genome[is.na(genome)]<-""


shiny::runApp(host="0.0.0.0",port=3168)
