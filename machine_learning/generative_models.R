.libPaths('C:/Users/BrBl1834/R/win-library')
library(dslabs)
library(caret)
library(tidyverse)
data("tissue_gene_expression")

# set.seed(1993) #if using R 3.5 or earlier
set.seed(1993, sample.kind="Rounding") # if using R 3.6 or later
ind <- which(tissue_gene_expression$y %in% c("cerebellum", "hippocampus"))
y <- droplevels(tissue_gene_expression$y[ind])
x <- tissue_gene_expression$x[ind, ]
x <- x[, sample(ncol(x), 10)]

# Q1 train on all data (overfitting!)
model<-train(x,y, method = "lda")

# Q2
model
model$finalModel$means
plotmodel<-data.frame(t(model$finalModel$means))
ggplot(data = plotmodel,aes(x=plotmodel$cerebellum,y=plotmodel$hippocampus)) + 
  geom_point() + geom_label(aes(label=rownames(plotmodel)))


# Q3

library(dslabs)      
library(caret)
data("tissue_gene_expression")

set.seed(1993) #set.seed(1993, sample.kind="Rounding") if using R 3.6 or later
ind <- which(tissue_gene_expression$y %in% c("cerebellum", "hippocampus"))
y <- droplevels(tissue_gene_expression$y[ind])
x <- tissue_gene_expression$x[ind, ]
x <- x[, sample(ncol(x), 10)]

model<-train(x,y, method = "qda")
model

# Q4

model
model$finalModel$means
plotmodel<-data.frame(t(model$finalModel$means))
ggplot(data = plotmodel,aes(x=plotmodel$cerebellum,y=plotmodel$hippocampus)) + 
  geom_point() + geom_label(aes(label=rownames(plotmodel)))

# Q5 preprocess: some predictors are low/high in both groups
model<-train(x,y, method = "lda", preProcess = 'center')
model
plotmodel<-data.frame(t(model$finalModel$means))
ggplot(data = plotmodel,aes(x=plotmodel$cerebellum,y=plotmodel$hippocampus)) + 
  geom_point() + geom_label(aes(label=rownames(plotmodel)))

# Q6
library(dslabs)      
library(caret)
data("tissue_gene_expression")

# set.seed(1993) # if using R 3.5 or earlier
set.seed(1993, sample.kind="Rounding") # if using R 3.6 or later
y <- tissue_gene_expression$y
x <- tissue_gene_expression$x
x <- x[, sample(ncol(x), 10)]
model<-train(x,y, method = "lda", preProcess = 'center')
model
plotmodel<-data.frame(t(model$finalModel$means))
ggplot(data = plotmodel,aes(x=plotmodel$cerebellum,y=plotmodel$hippocampus)) + 
  geom_point() + geom_label(aes(label=rownames(plotmodel)))