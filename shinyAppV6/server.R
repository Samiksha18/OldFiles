library(shiny)
library(ggplot2)
load("shiny_data.rda",.GlobalEnv)

function(input, output, session) {
  ig<-reactive(input$gene)
  outVar = reactive({
    shinyDf = get(input$dataset)
    paste(shinyDf$qv_anno, shinyDf$gene,shinyDf$desc, shinyDf$AltID)})
  observe({
    updateSelectInput(session, "gene",
                      choices = outVar() )})
  
  output$slicePlot <- renderPlot({
    #slice<-unlist(contrastdf[gsub("(.{10})(.*)","\\1", "QV15_06075 USA300HOU_1218"),cols2])
    slice<-unlist(contrastdf[gsub("(.{10})(.*)","\\1", input$gene),cols2])
    g<-ggplot(data.frame(slice))+
      geom_point(aes(x=cols2, y=slice), size=5)+
      scale_x_discrete(limits= ordering_val)+
      labs(title = paste(title_of_experiment, "Expression Overlay"), 
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
    
    h<-h+geom_point(data=shinyDf[shinyDf$qv_anno==as.character(gsub("(.{10})(.*)","\\1", input$gene)),], aes(x=volx, y=voly), size=5.5, colour="black", vjust=0)+
      scale_y_log10()+
      geom_hline(yintercept=-log(pval_thresh), color="green")+
      geom_vline(xintercept=c(-fcl2, fcl2), color="red")+
      #scale_x_log10()+
      labs(title="Volcano Plot",
           x="Log2 Fold Change; peaked genes in green",
           y=" -log10 BH adj. p-value ")+
      geom_point(data=shinyDf[shinyDf$qv_anno==as.character(gsub("(.{10})(.*)","\\1", input$gene)),], aes(x=peakvolx, y=peakvoly), size=5.5, colour="black", alpha = I(0.2), vjust=0)+
      geom_text(data=shinyDf[shinyDf$qv_anno==as.character(gsub("(.{10})(.*)","\\1", input$gene)),], aes(x=peakvolx, y=peakvoly+0.5, label=text), colour="black")
    h
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
  output$downloadData<-downloadHandler(
    filename=function(){
      paste("BAGEL_", input$dataset, ".csv",sep="")
    },
    content=function(con){
      write.csv(shinyDf,con)
    }  
  )
  ##############################  SHIVR###################################
  output$text1 <- renderText({paste("You have selected", input$locus_tag)})
  output$dnaseq <- renderText({genome[genome$locus_tag==input$locus_tag, "dnaseq"]})
  output$loc_start <- renderText({genome[genome$locus_tag==input$locus_tag, "loc_start"]})
  output$loc_end <- renderText({genome[genome$locus_tag==input$locus_tag, "loc_end"]})
  output$geneshiv <- renderText({genome[genome$locus_tag==input$locus_tag, "gene"]})
  output$direction <- renderText({genome[genome$locus_tag==input$locus_tag, "direction"]})
  output$note <- renderText({genome[genome$locus_tag==input$locus_tag, "note"]})
  output$product <- renderText({genome[genome$locus_tag==input$locus_tag, "product"]})
  output$protein_id <- renderText({genome[genome$locus_tag==input$locus_tag, "protein_id"]})
  output$db_xref <- renderText({genome[genome$locus_tag==input$locus_tag, "db_xref"]})
  output$translation <- renderText({genome[genome$locus_tag==input$locus_tag, "translation"]})
  output$dnaseq_upstream <- renderText({genome[genome$locus_tag==input$locus_tag, "dnaseq_upstream"]})
  output$dnaseq_downstream <- renderText({genome[genome$locus_tag==input$locus_tag, "dnaseq_downstream"]})
  output$scaffold <- renderText({genome[genome$locus_tag==input$locus_tag, "scaffold"]})
  output$type <- renderText({genome[genome$locus_tag==input$locus_tag, "type"]})
  output$range<- renderText(paste(genome[genome$locus_tag==input$locus_tag, "loc_start"], "to",
                                  genome[genome$locus_tag==input$locus_tag, "loc_end"], 
                                  "on",
                                  gsub("(.*)(\\.txt)","\\1",genome[genome$locus_tag==input$locus_tag, "scaffold"]),
                                  sep=" "))
  output$length<- renderText(as.numeric(genome[genome$locus_tag==input$locus_tag, "loc_end"])-
                               as.numeric(genome[genome$locus_tag==input$locus_tag, "loc_start"]))
  output$aaseq<-renderUI({
    HTML(paste(paste("Amino acid sequence for", input$locus_tag),
               genome[genome$locus_tag==input$locus_tag, "translation"],
               sep = '<br/>'))
  })
  output$dnaseq<-renderUI({
    HTML(paste(paste("DNA sequence for", input$locus_tag),
               genome[genome$locus_tag==input$locus_tag, "dnaseq"],
               sep = '<br/>'))
  })
  output$dnaseq_revcomp<-renderUI({
    (ifelse(genome[genome$locus_tag==input$locus_tag, "direction"]=="compliment",
    HTML(paste(paste("Reverse Complement:", input$locus_tag),"\n",
               as.character(reverseComplement(DNAString(genome[genome$locus_tag==input$locus_tag, "dnaseq"]))),
               sep = '')),
    HTML(paste(""))))
  })
}



