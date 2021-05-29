# when running on local
.libPaths('C:/Users/BrBl1834/R/win-library')

############# example: student scores across subjects correlate 
library(MASS)
library(rafalib)
# generate random data representing student scores
n <- 250
p <- 6
set.seed(1, sample.kind = "Rounding")
g <- mvrnorm(n,c(0,0),matrix(c(1,0.5,0.5,1),2,2))
Ystem <- g[,1] + matrix(rnorm(n*p/2,0,0.65),n,p/2)
Yhum <- g[,2] + matrix(rnorm(n*p/2,0,0.65),n,p/2)
Y <- cbind(Ystem,Yhum)
colnames(Y) <- c("Math","Science","CS","Eng","Hist","Classics")

# high correlations between subjects
round(cor(Y),2)

# plot correlation of data vs independent data
library(RColorBrewer)
mypar(1,2)
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
eps = matrix(rnorm(n*p),n,p)
par(mar = c(8.1, 8.1, 3.5, 2.1))
image(1:ncol(Y),1:ncol(Y),cor(Y)[,6:1],xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1),main="Actual Data")
axis(1,1:ncol(Y),colnames(Y),las=2)
axis(2,1:ncol(Y),rev(colnames(Y)),las=2)
image(1:ncol(Y),1:ncol(Y),cor(eps)[,6:1],xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1),main="Independent Data")
axis(1,1:ncol(Y),colnames(Y),las=2)
axis(2,1:ncol(Y),rev(colnames(Y)),las=2)

# image shows that correlation is high across all subjects, but even more so amon STEM/ humanities respectively
# --> there are hidden factors explaining these correlations (e.g. student academic performance, interest in STEM vs humanities)

# apparently, under certain assumptions, first two PC's can explain the student effect + interest effect

s <- svd(Y)
What <- t(s$v[,1:2])
colnames(What)<-colnames(Y)
round(What,2)

# first two PC's explain lot of the variation
fit = s$u[,1:2]%*% (s$d[1:2]*What)
var(as.vector(fit))/var(as.vector(Y))

# Correlated units --> linear models aren't appropriate anymore


############# example 2: cor structures in high throughput gene expression

library(Biobase)

library(GSE5859) #deprecated --> manually load the GSE5859.rda file
data(GSE5859)
wd<-getwd()
load(paste0(wd,"/GSE5859.rda"))

n <- nrow(pData(e))
o <- order(pData(e)$date)
Y=exprs(e)[,o]
cors=cor(Y-rowMeans(Y))
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)

mypar()
image(1:n,1:n,cors,xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1))

# --> several correlation structures, explained by 1,2,....,k factors
# how many k's should we pick?

############# exercises
rm(list = ls())
library(Biobase)
library(GSE5859Subset) #deprecated --> manually load the GSE5859Subset.rda file
data(GSE5859Subset)
wd<-getwd()
load(paste0(wd,"/GSE5859Subset.rda"))

y = geneExpression[,1:2]

# compare correlation structure between samples: no order vs ordered by date
y = geneExpression - rowMeans(geneExpression)
library(RColorBrewer)
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
image(1:ncol(y),1:ncol(y),cor(y)[],xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1),main="not ordered")

dateorder <- order(sampleInfo$date)
y_orderbydate <- y[,dateorder]
image(1:ncol(y_orderbydate),1:ncol(y_orderbydate),cor(y_orderbydate)[],xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1),main="ordered by date")

# using SVD to extract first two PCs
pcs = svd(y)$v[,1:2]

library(tidyverse)
p1dat<-data.frame(pcs[,1],sampleInfo$date,months=months(sampleInfo$date))
ggplot(data=p1dat,
       aes(x=pcs[,1], y=sampleInfo$date, col=months)) + 
  geom_point()
p1dat

p2dat <- data.frame(pcs[,2],sampleInfo$date,months=months(sampleInfo$date))
ggplot(data=p2dat,
       aes(x=pcs[,2], y=sampleInfo$date, col=months)) + 
  geom_point()
p2dat

o = order(sampleInfo$date)
cols = as.numeric(month)[o]
mypar(2,1)
for(i in 1:2){
  plot(pcs[o,i],col=cols,xaxt="n",xlab="")
  label = gsub("2005-","",sampleInfo$date[o])
  axis(1,1:ncol(y),label,las=2)
}

# variance explained per pc
s = svd(y)
varexplained = s$d^2/ sum(s$d^2)
plot(varexplained)
sum(varexplained>0.10)

# cor of pcs with months
month = factor( format(sampleInfo$date,"%m"))
cors = cor( as.numeric(month),s$v)
plot(t(cors))
which.max(abs(cors))
max(abs(cors))

# corr of pcs with sex
sex = factor(sampleInfo$group)
cors = cor( as.numeric(sex),s$v)
plot(t(cors))
which.max(abs(cors))
max(abs(cors))

# add first two pc's to remove batch effect
X <- model.matrix(~sex+s$v[,1:2])

y <- geneExpression[5,] #test
fit <- lm(y~X)
summary(fit)$coef[2,c(1,4)]

res <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit <- lm(y~X)
  summary(fit)$coef[2,c(1,4)] #keep only sex effect=> month effect is removed
} ) )
res<-data.frame(res)
names(res)<-c('est','pval')

library(qvalue)
qvals <- qvalue(res$pval)$qvalue

sum(qvals<0.1)

index <- which(qvals<0.1)

goi<-geneAnnotation[index,]
mean(goi$CHR %in% c('chrX','chrY'))
