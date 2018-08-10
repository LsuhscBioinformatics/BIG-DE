BIG-DE (BioInformatics and Genomics - Differential Expression) is an open source, Shiny-based application for exploring, visualizing, and analyzing gene expression data with minimal coding required. Users can use the application in their own RNA-seq or microarray expression data, and download the results. Built-in visualization tools include box plots, heatmaps and volcano plots. Please download the sample data provided to test the application.

## Access the application [here] (http://big.bio.lsuhsc.edu:3838/)


How to install the required packages before using the application?

```
source("https://bioconductor.org/biocLite.R")

biocLite("limma")

biocLite("preprocessCore")

install.packages("DT")

install.packages("ggplot2")

install.packages("gplots")
```

## Contributors
Main contributors: Tarun karthik kumar Mamidi, M.S.
Additional contributors: Jiande Wu PhD, Chindo Hicks PhD


Please report any issues or bugs to bioinfo.services@lsuhsc.edu

