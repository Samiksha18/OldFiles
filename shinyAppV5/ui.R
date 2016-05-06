#!/usr/bin/env Rscript
library(shiny)
library(ggplot2)
#source("server.R")
load("shiny_data.rda")

shinyUI(navbarPage("Brinsmade 'Omics Apps", 
                   tabPanel("Bagel",
fluidPage(    
  titlePanel("BAGEL: BrinsmAde Gene Expression visuaLizer"),
  
  sidebarLayout(          
    # Define the sidebar with one input
    sidebarPanel(
      selectInput("dataset","Subset", setList, selected = setList[1]),
      selectInput("gene", "Gene of Interest:", ""),
      #        paste(shinyDf$qv_anno, shinyDf$gene,shinyDf$desc, shinyDf$AltID)),#choices=fclist)
      hr(),
      helpText("Gene:"),
      helpText(verbatimTextOutput("geneInfo")),

      helpText("Fold Change (1 = no change):" ),
      helpText(verbatimTextOutput("sigLev")),
      plotOutput(("slicePlot")),
      downloadButton('downloadData', "Download Selected Subset"),
      helpText("(for the staunch of heart only)"),
      
    width=6),
    # Create a spot for the barplot
    mainPanel(
      plotOutput("volcanoPlot",click = "plot_click"),
      verbatimTextOutput("volcano_info"),
      width=6
      
    )
  )
)
)))


#library(shinyapps)
#shinyapps::setAccountInfo(name='srbdata', token='BE10526EB5667E21AEE2379BF9869296', secret='7/2aZj3by2hi1kjvL6aNrTkO9U0S7pQ/Y9qvEDxl')
#deployApp()
