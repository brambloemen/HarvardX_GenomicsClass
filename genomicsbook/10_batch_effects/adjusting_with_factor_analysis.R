# when running on local
# .libPaths('C:/Users/BrBl1834/R/win-library')

# load GSE5859subset
wd<-getwd()
load(paste0(wd,"/GSE5859Subset.rda"))


# course rmd code ---------------------------------------------------------
# plot images
library(rafalib)
library(RColorBrewer)
library(genefilter)
sex <- sampleInfo$group
batch <- factor(format(sampleInfo$date,"%m"))
chr <- geneAnnotation$CHR
tt<-rowttests(geneExpression,batch)
ind1 <- which(chr=="chrY") #real differences
ind2 <- setdiff(c(order(tt$dm)[1:25],order(-tt$dm)[1:25]),ind1)
set.seed(1)
ind0 <- setdiff(sample(seq(along=tt$dm),50),c(ind2,ind1))
geneindex<-c(ind2,ind0,ind1)
mat<-geneExpression[geneindex,]
mat <- mat -rowMeans(mat)
icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
mypar(1,2)




# plot month and sex effects
image(t(mat),xaxt="n",yaxt="n",col=icolors)
y <- geneExpression - rowMeans(geneExpression)

# plot sample to sample correlation
image(1:ncol(y),1:ncol(y),cor(y),col=icolors,zlim=c(-1,1),
      xaxt="n",xlab="",yaxt="n",ylab="")
axis(2,1:ncol(y),sex,las=2)
axis(1,1:ncol(y),sex,las=2)


# months
times <-sampleInfo$date 
mypar(1,1)
o=order(times)
plot(times[o],pch=21,bg=as.numeric(batch)[o],ylab="date")
o=order(times)
plot(times[o],pch=21,bg=as.numeric(batch)[o],ylab="date")

# PCA to see if day has effect
s <- svd(y)
mypar(1,1)
o<-order(times)
cols <- as.numeric( batch)
plot(s$v[o,1],pch=21,cex=1.25,bg=cols[o],ylab="First PC",xaxt="n",xlab="")
legend("topleft",c("Month 1","Month 2"),col=1:2,pch=16,box.lwd=0)

# how much variability is explained by the first PC
mypar(1,1)
plot(s$d^2/sum(s$d^2),ylab="% variance explained",xlab="Principal component")

# first 6 PC's are at least partially driven by date: they vary significantly between dates
mypar(3,4)
for(i in 1:12){
  days <- gsub("2005-","",times)  
  boxplot(split(s$v[,i],gsub("2005-","",days)))
}

# remove first 6 PC's, perform rowttests
D <- s$d; D[1:4]<-0 #take out first 2
cleandat <- sweep(s$u,2,D,"*")%*%t(s$v)
res <-rowttests(cleandat,factor(sex))
library(qvalue)
mypar(1,2)
hist(res$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))

plot(res$dm,-log10(res$p.value))
points(res$dm[which(chr=="chrX")],-log10(res$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res$dm[which(chr=="chrY")],-log10(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)

# removing first 6 PC's also removes a lot of the biological variation
qvals <- qvalue(res$p.value)$qvalue
index <- which(qvals<0.1)

cat("Total genes with q-value < 0.1: ",length(index),"\n",
    "Number of selected genes on chrY: ", sum(chr[index]=="chrY",na.rm=TRUE),"\n",
    "Number of selected genes on chrX: ", sum(chr[index]=="chrX",na.rm=TRUE),sep="")



# SVA: fit model with covariate of interest, as well as with batch covariates
library(sva)
library(limma)
mod <- model.matrix(~sex)
cind <- order( as.Date(sampleInfo$date) )
dates <- gsub("2005-","",sampleInfo$date)
weights=rep(1,nrow(y))
par(mar = c(4.1, 2.1, 3.5, 2.1), 
    mgp = c(1.5, 0.5, 0))
layout(matrix(c(1:6),nrow=2,byrow=TRUE),widths=c(5,1.5,5))
for(b in 1:2){
  image(1:ncol(mat),1:nrow(mat),t(mat[,cind]*weights[geneindex]),xaxt="n",yaxt="n",col=icolors,xlab="",ylab="")
  axis(side=1,seq(along=dates),dates[cind],las=2)
  abline(v=12.5)
  
  svafit <- sva(y,mod,B=b,n.sv=5)
  weights = svafit$pprob.gam*(1-svafit$pprob.b)
  
  surrogate <- svd( y*weights)$v[,1]#Weighted SVD
  
  image(matrix(weights[geneindex],nrow=1),xaxt="n",yaxt="n",col=brewer.pal(9,"Blues"))
  plot(surrogate[cind],bg=sex[cind]+1,pch=21,xlab="",xaxt="n",ylab="Surrogate variable",ylim=c(-.5,.5),cex=1.5)
  axis(side=1,seq(along=dates),dates[cind],las=2)
  abline(v=12.5)
  text(1,0.5,"June")
  text(13.5,0.5,"Oct")
  legend("bottomright",c("0","1"),col=c(1,2),pch=16)
}

# let SVA function find significant surrogate vars
library(limma)
svafit <- sva(geneExpression,mod)

# make model matrix of outcome + surr vars, use to fit lm
svaX<-model.matrix(~sex+svafit$sv)
lmfit <- lmFit(geneExpression,svaX)
tt<- lmfit$coef[,2]*sqrt(lmfit$df.residual)/(2*lmfit$sigma)

res <- data.frame(dm= -lmfit$coef[,2],
                  p.value=2*(1-pt(abs(tt),lmfit$df.residual[1]) ) )
mypar(1,2)
hist(res$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))

plot(res$dm,-log10(res$p.value))
points(res$dm[which(chr=="chrX")],-log10(res$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res$dm[which(chr=="chrY")],-log10(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)

qvals <- qvalue(res$p.value)$qvalue
index <- which(qvals<0.1)

cat("Total genes with q-value < 0.1: ",length(index),"\n",
    "Number of selected genes on chrY: ", sum(chr[index]=="chrY",na.rm=TRUE),"\n",
    "Number of selected genes on chrX: ", sum(chr[index]=="chrX",na.rm=TRUE),sep="")


Batch<- lmfit$coef[geneindex,3:7]%*%t(svaX[,3:7])
Signal<-lmfit$coef[geneindex,1:2]%*%t(svaX[,1:2])
error <- geneExpression[geneindex,]-Signal-Batch
##demean for plot
Signal <-Signal-rowMeans(Signal)
mat <- geneExpression[geneindex,]-rowMeans(geneExpression[geneindex,])
mypar(1,4,mar = c(2.75, 4.5, 2.6, 1.1))
image(t(mat),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")
image(t(Signal),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")
image(t(Batch),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")
image(t(error),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")


# exercises ---------------------------------------------------------
library(Biobase)

# SVD: first factor highly correlates with sex
s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])

# SVA: estimate factors, but downweighs outcome of interest
sex = sampleInfo$group
mod = model.matrix(~sex)
svafit = sva(geneExpression,mod)
head(svafit$sv)

# compare SVA factors with PC from SVD
for(i in 1:ncol(svafit$sv)){
  print( cor(s$v[,i],svafit$sv[,i]) )
}

# fit linear model based on sex (outcome of interest) + 5 factors from SVA

svaX<-model.matrix(~sex + svafit$sv)
lmfit <- lmFit(geneExpression,svaX)
tt<- lmfit$coef[,2]*sqrt(lmfit$df.residual)/(2*lmfit$sigma)

res <- data.frame(dm= -lmfit$coef[,2],
                  p.value=2*(1-pt(abs(tt),lmfit$df.residual[1]) ) )
mypar(1,2)

qvals <- qvalue(res$p.value)

# how many have qval<0.1
sum(qvals$qvalues<0.1)

# on which chromosomes?
genes_qval_below_0.1<-geneAnnotation[which(qvals$qvalues<0.1),]
genes_qval_below_0.1
