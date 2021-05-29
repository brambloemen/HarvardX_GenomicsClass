.libPaths("C:/Users/BrBl1834/R/win-library")

library(dslabs)
data(tissue_gene_expression)

dim(tissue_gene_expression$x)

table(tissue_gene_expression$y)

# Q1
d <- dist(tissue_gene_expression$x)

# Q2
d <- as.matrix(d)

d[1:2,1:2]
d[39:40,39:40]
d[73:74,73:74]

d[1,]
plot(d[1,])

