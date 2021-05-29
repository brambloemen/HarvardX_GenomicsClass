.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(dslabs)

# Q1
d <- dist(tissue_gene_expression$x - rowMeans(tissue_gene_expression$x))

# Q2
h <- hclust(d)

plot(h, cex = 0.5, main = "", xlab = "")

# Q3
library(RColorBrewer)
sds <- matrixStats::colSds(tissue_gene_expression$x)
ind <- order(sds, decreasing = TRUE)[1:50]
colors <- brewer.pal(7, "Dark2")[as.numeric(tissue_gene_expression$y)]

heatmap(t(tissue_gene_expression$x[,ind]), col = brewer.pal(11, "RdBu"), scale = "row", ColSideColors = colors)
