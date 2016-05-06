library(shiny)
library(ggplot2)
library(dplyr)
# shinyDf<-read.csv("../shinyApps/v5/shinyDF.csv")
# contrastdf<-read.csv("../shinyApps/v5/contrastdf.csv")
# resdf<-read.csv("../shinyApps/v5/resdf.csv")
# key<-read.csv("../shinyApps/v5/key.csv")
#attach("../shinyApps/v5/shiny_data.rda")



#Define a server for the Shiny app
server<-function(input, output) {
  attach("shiny_data.rda")
  ig<-reactive(input$gene)
  output$slicePlot <- renderPlot({
    slice<-unlist(contrastdf[gsub("(.{10})(.*)","\\1", input$gene),cols2])
    g<-ggplot(data.frame(slice))+
      geom_point(aes(x=cols2, y=slice), size=5)+
      scale_x_discrete(limits= ordering_val)+
      labs(title = paste( "Expression Overlay"), 
           y="% expression",
           x="CodY variant")
    
    g+stat_smooth(aes(x=cols2, y=slice, group=1))
  })
  output$volcanoPlot<-renderPlot({
    ###Source:  http://bioinformatics.knowledgeblog.org/2011/06/21/volcano-plots-of-microarray-data/
    ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    #     no_of_genes = length(rownames(dds))
    xlm<-c(-5, 5)
    ylm<-c(0, 10)
    #     shinyDf<-data.frame(aRes)
    #     shinyDf$threshold = as.factor((abs(shinyDf$log2FoldChange) > fcl2) & (shinyDf$padj < pval_thresh/no_of_genes))
    #     shinyDf$qv_anno<-gsub("(.{10})(.*)", "\\1", row.names(shinyDf))
    #     shinyDf<-merge(shinyDf, data.frame(qv_anno=combined$qv_anno, desc=combined$sax_kegg_anno), by="qv_anno", all.x=T)
    #     shinyDf<-na.omit(shinyDf)
    h<-ggplot(data=shinyDf, aes(x=volx, y=voly))+#+, colour=threshold)) +
      geom_point(alpha=0.4, size=4, pch=21 ,aes(fill=threshold, color=threshold)) +
      geom_point(data=shinyDf[shinyDf$peak=="peaked",], size=5, pch=21 ,alpha=1, aes(fill=NULL), color="chartreuse1") +
      theme(legend.position = "none")
    #+xlim(xlm) +ylim(ylm)
    
    h+geom_point(data=shinyDf[shinyDf$qv_anno==as.character(gsub("(.{10})(.*)","\\1", input$gene)),], aes(x=volx, y=voly), size=5.5, colour="black", vjust=0)+
      geom_hline(yintercept=.01, color="green")+
      geom_vline(xintercept=c(-1, 1), color="red")+
      scale_y_log10()+
      #scale_x_log10()+
      labs(title="Volcano Plot",
           x="Log2 Fold Change; peaked genes in green",
           y=" -log10 BH adj. p-value ")
  })
  output$geneInfo <- renderText({
    paste(shinyDf[shinyDf$qv_anno==as.character(gsub("(.{10})(.*)","\\1", input$gene)), "desc"])
  })
  output$sigLev <- renderText({
    paste(shinyDf[shinyDf$qv_anno==gsub("(.{10})(.*)","\\1", input$gene), "fc"])
  })
  output$volcano_info <- renderPrint({
    # With base graphics, need to tell it what the x and y variables are.
    nearPoints(shinyDf, input$plot_click,threshold = 10, maxpoints = 1,
               addDist = TRUE)})
}




ui<-fluidPage(    
  titlePanel("Brinsmade Lab Gene Expression Browser"),
  sidebarLayout(          
    # Define the sidebar with one input
    sidebarPanel(
      selectInput("gene", "Gene of Interest:", 
                  paste(shinyDf$qv_anno, shinyDf$desc)),#choices=fclist)
      hr(),
      helpText("Gene:"),
      helpText(verbatimTextOutput("geneInfo")),
      helpText("Maximum Fold Change Threshold:" ),
      helpText(verbatimTextOutput("sigLev"))
    ),
    # Create a spot for the barplot
    mainPanel(
      plotOutput("slicePlot"),
      plotOutput("volcanoPlot",click = "plot_click"),
      verbatimTextOutput("volcano_info")
      
    )
  )
)

shinyAppV5<-shinyApp(ui=ui, server=server, )
# runApp(Version_5, launch.browser=T)
#########
#install.packages('devtools')


