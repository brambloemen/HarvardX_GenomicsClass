.libPaths('C:/Users/BrBl1834/R/win-library')
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

g <- sampleInfo$group
g

e <- geneExpression[25,]

library(rafalib)
mypar(1,2)

qqnorm(e[g==1])
qqline(e[g==1])

qqnorm(e[g==0])
qqline(e[g==0])

t.test(e[g==1],e[g==0])$p.value


myttest <- function(x) t.test(x[g==1],x[g==0],var.equal=TRUE)$p.value
pvals <- apply(geneExpression,1,myttest)

sum(pvals<0.05)

# same amount of t-tests on random data (where null hypothesis is always true) also results in lots of p<0.05
# --> to be expected when performing many t tests
# !!!! P value isn't everything
set.seed(1)
m <- nrow(geneExpression)
n <- ncol(geneExpression)
randomData <- matrix(rnorm(n*m),m,n)
nullpvals <- apply(randomData,1,myttest)
sum(nullpvals<0.05)


# bioc genefilter package includes rowttest function to more efficiently perform the t tests for each gene
library(genefilter)
results <- rowttests(geneExpression,factor(g))
max(abs(pvals-results$p))