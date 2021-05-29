.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(dslabs)

data("tissue_gene_expression")
dim(tissue_gene_expression$x)

# Q1

t<-rbind(tissue_gene_expression$y, tissue_gene_expression$x)
t<-t[2:190,]
s<-svd(t)

# pc2<-s$u[,1:2] * s$d[1:2] %*% t(s$v[,1:2])
pc2 <- sweep(s$u[, 1:2], 2, s$d[1:2], FUN="*")
plot(pc2)

pc2 <- data.frame(tissue=tissue_gene_expression$y,pc2)
ggplot(pc2, aes(x=pc2[,2],y=pc2[,3],col=tissue)) + geom_point()


# solution: use function PRcomp
pc <- prcomp(tissue_gene_expression$x)
data.frame(pc_1 = pc$x[,1], pc_2 = pc$x[,2], 
           tissue = tissue_gene_expression$y) %>%
  ggplot(aes(pc_1, pc_2, color = tissue)) +
  geom_point()

# Q2
obs_means <- rowMeans(tissue_gene_expression$x)

plot(pc2[,2],obs_means)
cor(pc2[,2],obs_means)

# solution
avgs <- rowMeans(tissue_gene_expression$x)
data.frame(pc_1 = pc$x[,1], avg = avgs, 
           tissue = tissue_gene_expression$y) %>%
  ggplot(aes(avgs, pc_1, color = tissue)) +
  geom_point()
cor(avgs, pc$x[,1])

# Q3
x <- with(tissue_gene_expression, sweep(x, 1, rowMeans(x)))
pc <- prcomp(x)
data.frame(pc_1 = pc$x[,1], pc_2 = pc$x[,2], 
           tissue = tissue_gene_expression$y) %>%
  ggplot(aes(pc_1, pc_2, color = tissue)) +
  geom_point()

pc$x[,7]

p_pc7<-data.frame(tissue=tissue_gene_expression$y,pc$x[,7])

ggplot(p_pc7, aes(tissue, p_pc7[,2])) + geom_boxplot()

summary(pc)
