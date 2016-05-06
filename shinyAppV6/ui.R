#!/usr/bin/env Rscript
library(shiny)
library(ggplot2)
#source("server.R")
load("shiny_data.rda")

shinyUI(navbarPage("Brinsmade 'Omics Apps", 
                   tabPanel("BAGEL",
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
),tabPanel("shivR",
           fluidPage(
             title = 'ShiVrMeGenbanks',
             
             sidebarLayout(
               sidebarPanel(
                 selectInput(
                   'locus_tag', 'locus_tag',  choices = unlist(genome$locus_tag),
                   selected = genome$locus_tag[2]),
                 helpText("Coords:"),
                 helpText(verbatimTextOutput("range")),
                 helpText("Length:"),
                 helpText(verbatimTextOutput("length")),
                 helpText("product:"),
                 helpText(verbatimTextOutput("product")),
                 helpText("direction:"),
                 helpText(verbatimTextOutput("direction")),
                 helpText("protein_id:"),
                 helpText(verbatimTextOutput("protein_id")),
#                 helpText("revcomp:"),
#                 helpText(verbatimTextOutput("danseq_revcomp")),
                 helpText("WARNING:  DNA SEQUENCE DISPLAYED IS ON THE LEADING STRAND; 
               IF SELECTED GENE IS ON COMPLEMENT STRAND, it will be displayed in CAPS")
                 
               ),
               mainPanel(
                 titlePanel("Genome Browser"),
                 mainPanel(htmlOutput("aaseq"),
                           htmlOutput("dnaseq"),
                           htmlOutput("dnaseq_revcomp")),
              #            ifelse(genome[genome$locus_tag==input$locus_tag, "direction"]=="complement",
              #                   htmlOutput("aaseq"), htmlOutput("dnaseq_revcomp"))),
              #                    htmlOutput()
                 
                 width = 6)
               
             )
           )
           )
))
           

