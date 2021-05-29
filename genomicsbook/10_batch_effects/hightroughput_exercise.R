.libPaths('C:/Users/BrBl1834/R/win-library')
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

# 1
library(tidyverse)
sampleInfo$date <- as.Date(sampleInfo$date)
sampleInfo%>%filter(date==as.Date("2005-06-27", "%Y-%m-%d"))

# 2
geneAnnotation%>%group_by(CHR)%>%summarize(n())
geneAnnotation%>%filter(CHR=="chrY")%>%count()

# 3
sampleInfo$date <- as.Date(sampleInfo$date)
sampleInfo%>%filter(date==as.Date("2005-06-27", "%Y-%m-%d"))
subject<-sampleInfo%>%filter(date==as.Date("2005-06-10", "%Y-%m-%d"))

dim(geneExpression)

subject_expression <- geneExpression[,subject$filename]
gene<-"ARPC1A"
gene_an<-geneAnnotation%>%filter(SYMBOL==gene)
subject_expression[which(names(subject_expression)==gene_an$PROBEID)]


