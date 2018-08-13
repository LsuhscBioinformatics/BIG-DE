# Packages required
library("shiny")
library("limma")
library(preprocessCore)
library(gplots)
library(ggplot2)
library(RColorBrewer)

#Size of upload file up to 300MB
options(shiny.maxRequestSize=300*1024^2) 

#server function
server = function(input, output) {

#Read Expresssion input from user  
  df <- reactive({ 
    req(input$file)
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   row.names = 1,
                   check.names=F)
    
    return(df)
  })
  
#Read Pheno input from user  
  pheno <- reactive({ 
    req(input$file1)
    pheno <- read.csv(input$file1$datapath,
                      header = TRUE,
                      sep = input$sep1,
                      
                      check.names=F)
    
    return(pheno)
  })
  
  
 #Ouput the head of Expression file for display 
  output$table <- renderTable({
    req(input$file)
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   
                   check.names=F)
    
    return(head(df))
  })
  
  #Output example dataset to display on home page
  #Iris data from library(DT)
  output$mytable3 <- DT::renderDataTable({
    DT::datatable(iris[,-5], colnames = c('Genes','Sample1', 'Sample2', 'Sample3', 'Sample4'), options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  })
  
  #select normalization method based on user selection
  norm <- reactive({ 
    if(input$norm == "0"){
      #raw data
      norm <- df()
    }
    if(input$norm == "1"){
      #log2 transform
      norm <- log2(df())
    }
    if(input$norm == "2"){
      #Quantile normalization
      norm <- normalize.quantiles(as.matrix(df()))
      colnames(norm) <- colnames(df())
    }
    if(input$norm == "3"){
      # 1+log2 transform
      norm <- log2(df()+1)
    }
    rownames(norm) <- rownames(df())
    return(norm)
  }) 
  
  #Print the normalized data to a file named "normdata.csv" 
    output$normdata <- downloadHandler(
    filename = function() {
      paste("normdata.csv", sep = "")
    },
    content = function(file) {
      write.csv(norm(), file, row.names = TRUE)
    }
  )
  
  #Function to output boxplot based on the normalization selected
  plotBox <- function(){
      boxplot(norm(), main=input$boxtitle, xlab= input$xbox, ylab= input$ybox )+
      geom_point()
  }
  
  #Output the boxplot generated above
  output$boxplot <- renderPlot({
    plotBox()
  }) 
  
  #Differential Expression using Limma package
  #Function to output the model fit for volcano plot
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
  
  #Output volcano plot from the result generated above
  output$volcanoplot <- renderPlot({
    volcanoplot(volcano(),coef=1, main=input$voltitle)
  })
  
  #print the results from Differential expression analysis to a file named "Results.csv"
  output$voldata <- downloadHandler(
    filename = function() {
      paste("Results.csv", sep = "")
    },
    content = function(file) {
      write.csv(volcano(), file, quote = F, row.names = TRUE)
    }
  )
  
  #Take the LogFC filter from user and generate heatmap based on the D.E. genes
  output$volheatmap <- renderPlot({
    probeset.list1 <- topTable(volcano(), coef=1, n="Inf" , lfc =input$lfc)
    genes <- row.names(probeset.list1)
    norm.data <- norm()[genes,]
    heatmap.2(as.matrix(norm.data),col=rev(color()(50)),trace="none", main=input$heattitle1,
              scale="row" ,Colv="NA")
  })
  
#Function to select number of genes with most variation
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
  
  #Variable to select different colors by user
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
  
  #Generate heatmap for the most variable genes generated by function heat
  output$heatmap <- renderPlot({
    heatmap.2(as.matrix(heat()),col=rev(color()(50)),trace="none", main=input$heattitle,
              scale="row" ,Colv="NA")
  })
  
  #Print the gene expression of most variable genes to a file named "variance.csv"
  output$vardata <- downloadHandler(
    filename = function() {
      paste("variance.csv", sep = "")
    },
    content = function(file) {
      write.csv(heat(), file, row.names = TRUE)
    }
  )
}
