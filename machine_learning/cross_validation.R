.libPaths('C:/Users/BrBl1834/R/win-library')

library(tidyverse)
library(caret)

# set.seed(1996) #if you are using R 3.5 or earlier
set.seed(1996, sample.kind="Rounding") #if you are using R 3.6 or later
n <- 1000
p <- 10000
x <- matrix(rnorm(n*p), n, p)
colnames(x) <- paste("x", 1:ncol(x), sep = "_")
y <- rbinom(n, 1, 0.5) %>% factor()

# take 100 random columns, use these as predictors in Q1
x_subset <- x[ ,sample(p, 100)]


# Q1
fit <- train(x_subset, y, method = "glm")
fit$results

# Q2 use t test for each predictor with sample groups being y=1 and y=0
install.packages("BiocManager")
BiocManager::install("genefilter")
library(genefilter)
tt <- colttests(x, y)

pvals <- tt$p.value

# Q3 which predictors are statistically significant
ind<-which(pvals<0.01)
length(ind)

# Q4 rerun cross-validation using only the significant columns defined in Q3
x_subset <- x[ ,ind]

fit <- train(x_subset, y, method = "glm")
fit$results

# Q5 this time, use knn for modelling
k <- seq(101, 301, 25)

fit <- train(x_subset, y, method = "knn", tuneGrid = data.frame(k = seq(101, 301, 25)))
ggplot(fit)

# Q6 - non-code question
# Q7
library(dslabs)
data(tissue_gene_expression)

fit <- train(tissue_gene_expression$x, tissue_gene_expression$y, method = "knn", tuneGrid = data.frame(k = seq(1,7,2)))
ggplot(fit)
fit$results
