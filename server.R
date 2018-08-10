library("shiny")
library("limma")
library(preprocessCore)
library(gplots)
library(ggplot2)
library(RColorBrewer)

options(shiny.maxRequestSize=300*1024^2) 
server = function(input, output) {

  
  df <- reactive({ 
    req(input$file)
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   row.names = 1,
                   check.names=F)
    
    return(df)
  })
  pheno <- reactive({ 
    req(input$file1)
    pheno <- read.csv(input$file1$datapath,
                      header = TRUE,
                      sep = input$sep1,
                      
                      check.names=F)
    
    return(pheno)
  })
  
  
  
  output$table <- renderTable({
    req(input$file)
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   
                   check.names=F)
    
    return(head(df))
  })
  output$mytable3 <- DT::renderDataTable({
    DT::datatable(iris[,-5], colnames = c('Genes','Sample1', 'Sample2', 'Sample3', 'Sample4'), options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  })
  
  norm <- reactive({ 
    if(input$norm == "0"){
      norm <- df()
    }
    if(input$norm == "1"){
      norm <- log2(df())
    }
    if(input$norm == "2"){
      norm <- normalize.quantiles(as.matrix(df()))
      colnames(norm) <- colnames(df())
    }
    if(input$norm == "3"){
      norm <- log2(df()+1)
    }
    rownames(norm) <- rownames(df())
    return(norm)
  }) 
  plotBox <- function(){
    
    boxplot(norm(), main=input$boxtitle, xlab= input$xbox, ylab= input$ybox )+
      geom_point()
  }
  output$boxplot <- renderPlot({
    plotBox()
    
  }) 
  volcano <- function(){
    group <- factor(pheno()[,2])
    design <- model.matrix(~ 0 + group)
    ## Make the column names of the design matrix a bit nicer
    colnames(design) <- levels(group)
    fit <- lmFit(norm(),design)
    fit.cont <- contrasts.fit(fit, c(-1,1))
    fit.cont <- eBayes(fit.cont)
    
    return(fit.cont)
  }
  output$volcanoplot <- renderPlot({
    
    volcanoplot(volcano(),coef=1, main=input$voltitle)
  })
  
  output$voldata <- downloadHandler(
    
    filename = function() {
      paste("Results.csv", sep = "")
    },
    content = function(file) {

      write.csv(volcano(), file, quote = F, row.names = TRUE)
    }
  )
  
  output$volheatmap <- renderPlot({
    probeset.list1 <- topTable(volcano(), coef=1, n="Inf" , lfc =input$lfc)
    genes <- row.names(probeset.list1)
    #genes <- probeset.list1[,1]
    norm.data <- norm()[genes,]
    heatmap.2(as.matrix(norm.data),col=rev(color()(50)),trace="none", main=input$heattitle1,
              scale="row" ,Colv="NA")
  })
  
  output$normdata <- downloadHandler(
    
    filename = function() {
      paste("normdata.csv", sep = "")
    },
    content = function(file) {
      write.csv(norm(), file, row.names = TRUE)
    }
  )
  
  
  heat <- reactive({ 
    
    logcounts <- norm()
    rownames(logcounts) <- rownames(df())
    var_genes <- apply(logcounts, 1, var)
    
    # Get the gene names for the top 500 most variable genes
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$vargenes]
    
    # Subset logcounts matrix
    highly_variable_lcpm <- logcounts[select_var,]
    
    
    
    return(highly_variable_lcpm)
  })
  color <- reactive({ 
    if(input$color == "RedYellowGreen"){
      mypalette <- brewer.pal(11,"RdYlGn")
      morecols <- colorRampPalette(mypalette)
    }
    if(input$color == "GreenBlue"){
      morecols <- colorRampPalette(c("blue", "white", "green"))
    }
    if(input$color == "RedYellowBlue"){
      mypalette <- brewer.pal(11,"RdYlBu")
      morecols <- colorRampPalette(mypalette)
    }
    if(input$color == "RedBlue"){
      
      morecols <- colorRampPalette(c("blue", "white", "red"))
    }
    return(morecols)
  }) 
  output$heatmap <- renderPlot({
    
    heatmap.2(as.matrix(heat()),col=rev(color()(50)),trace="none", main=input$heattitle,
              scale="row" ,Colv="NA")
  }
  )
  output$vardata <- downloadHandler(
    
    filename = function() {
      paste("variance.csv", sep = "")
    },
    content = function(file) {
      write.csv(heat(), file, row.names = TRUE)
    }
  )
  
  
}
