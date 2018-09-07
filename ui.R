library("shiny")
library(DT)

#User interface 
ui = tagList(
  
  #Navigation bar for different tabs
  navbarPage(
    # theme = "cerulean",  # <--- To use a theme, uncomment this
    "BIG-DE",
    tabPanel("About", 
             sidebarPanel(
               tags$h4("BIG-DE (BioInformatics and Genomics - Differential Expression) is an open source, Shiny-based application for exploring, visualizing, and analyzing gene expression data with minimal coding required. Users can use the application in their own RNA-seq or microarray expression data, and download the results. Built-in visualization tools include box plots, heatmaps and volcano plots."),
               tags$hr(),
               tags$h4("Please look at the sample data"),
               tags$hr(),
               tags$h4("Please make a separate Pheno table for Differential expression as shown below (.csv file)"),
               tags$h4("This application supports only 2 classes for now. We'll improve this application in near future."),
               tags$br(),
               tags$h4("Sample     ,     Class/Time"),
               tags$h4("sample1    ,     Tumor/Day0"),
               tags$h4("sample2    ,     Tumor/Day0"),
               tags$h4("sample3    ,     Normal/Day1"),
               tags$h4("sample4    ,     Normal/Day1"),
               tags$hr(),
               tags$div(
                 "Please contact", 
                 tags$a(href="https://www.medschool.lsuhsc.edu/bioinformatics/contact.aspx", "Bioinformatics group"),
                 "for any queries "
               )),
               
             #output example data 
             mainPanel(
               DT::dataTableOutput("mytable3")
             )
    ),
    
    #Tab for data upload
    tabPanel("Data Upload",
             sidebarPanel(
               
               fileInput("file", "Upload Expression File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               # Input: Checkbox if file has header ----
               checkboxInput("header", "Header", TRUE),
               
               # Input: Select separator ----
               radioButtons("sep", "Separator",
                            choices = c(Comma = ",",
                                        Semicolon = ";",
                                        Tab = "\t"),
                            selected = ",")
               
               
             ),
             
             #Output the uploaded file 
             mainPanel(
               
               h4("Table Header"),
               tableOutput("table")
               )
             
    ),
    
    #Tab for normalization
    tabPanel("Normalization", 
             sidebarPanel(
               
               #Button to select the type of normalization method used
               radioButtons("norm", "Normalization Method",
                            choices = c("Raw data" = "0",
                                        "Log2" = "1",
                                        "Quantile" = "2",
                                        "1+Log2" = "3"),
                            selected = "0"),
               textInput("boxtitle", "Title:", "Boxplot"),
               textInput("xbox", "X-axis:", "Samples"),
               textInput("ybox", "Y-axis:", "Expression value"),
               downloadButton("normdata", "Download normalized data")
             ),
             
             #Output boxplot for the normalized method used
             mainPanel(
               plotOutput("boxplot", width = "100%", height = 600)
             )
    ),
    
    #Tab for D.E. analysis
    tabPanel("Differential Expression", 
             sidebarPanel(
               conditionalPanel(condition="input.tabselected==1",
               tags$h4("Differential Expression using Limma"),
               tags$hr(),
                                
               #Upload Pheno file for class/time variable
               fileInput("file1", "Upload Pheno File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
                                
               # Input: Checkbox if file has header ----
               checkboxInput("header1", "Header", TRUE),
               
               # Input: Select separator ----
               radioButtons("sep1", "Separator",
                            choices = c(Comma = ",",
                                        Semicolon = ";",
                                        Tab = "\t"),
                            selected = ","),
              #Title for the plot
               textInput("voltitle", "Title:", "Volcano plot"),
              #Download button for normalized data
               downloadButton("voldata", "Download D.E. results")
               ),
               
               conditionalPanel(condition="input.tabselected==2",
               tags$h4("Heatmap"),
               tags$hr(),
               textInput("lfc", "LogFC cutoff", "2"),

               textInput("heattitle1", "Title of Heatmap:", "Tumor vs Normal"),
               
               selectInput("color", "Choose color:",
                           choices = c("RedYellowGreen", "GreenBlue", "RedYellowBlue", "RedBlue"))
               )
             ),
             
             #Output Volcano plot and Heatmap
             mainPanel(
               #tableOutput("contrast"),
               tabsetPanel(type = "tabs",
               tabPanel("Volcano", value = 1, plotOutput("volcanoplot", width = "100%", height = 600)),
               tabPanel("Heatmap", value = 2,  plotOutput("volheatmap", width = "100%", height = 600)),
               id= "tabselected"
               )
             )
    ),
    
    #Tab for variance and Heatmap
    tabPanel("Variance", 
             sidebarPanel(
               textInput("heattitle", "Title of Heatmap:", "Top 500 most variable genes across samples"),
               textInput("vargenes", "Top __ variable genes:", "500"),
               selectInput("color", "Choose color:",
                           choices = c("RedYellowGreen", "GreenBlue", "RedYellowBlue", "RedBlue")),
               downloadButton("vardata", "Download top variance data")

             ),
             #Output heatmap
             mainPanel(
               plotOutput("heatmap", width = "100%", height = 600)
             )
    )
  )
)
