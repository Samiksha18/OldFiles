#!/usr/bin/Rscript
##  SHiny VieweR for MyE Genbank Sequences

#  This file underwent a major revision on 20151208 with the completion of version
#  0.2.0 of the gbparse project.  this will now solely house the Shiny bits for 
#  browsing the GenomeInfoDb
library(shiny)
#dir.create("~/GitHub/R/shinyApps/shiVr/")
setwd("~/GitHub/R/shinyApps/shiVr/")
shiny::runApp(host="0.0.0.0",port=3169)
