.libPaths('C:/Users/BrBl1834/R/win-library')

library(devtools)
library(Biobase) ##available from Bioconductor
library(genefilter) 
devtools::install_github(GSE5859)
library(GSE5859) ##available from github
GSE5859
data(GSE5859)

geneExpression = exprs(e)
sampleInfo = pData(e)
head(sampleInfo$date)

year = factor( format(sampleInfo$date,"%y") )
tab = table(year,sampleInfo$ethnicity)
print(tab)

library(genefilter)

##remove control genes
out <- grep("AFFX",rownames(geneExpression))

eth <- sampleInfo$ethnicity
ind<- which(eth%in%c("CEU","ASN"))
res1 <- rowttests(geneExpression[-out,ind],droplevels(eth[ind]))
ind <- which(year%in%c("02","03") & eth=="CEU")
res2 <- rowttests(geneExpression[-out,ind],droplevels(year[ind]))

XLIM <- max(abs(c(res1$dm,res2$dm)))*c(-1,1)
YLIM <- range(-log10(c(res1$p,res2$p)))
mypar(1,2)
plot(res1$dm,-log10(res1$p),xlim=XLIM,ylim=YLIM,
     xlab="Effect size",ylab="-log10(p-value)",main="Populations")
plot(res2$dm,-log10(res2$p),xlim=XLIM,ylim=YLIM,
     xlab="Effect size",ylab="-log10(p-value)",main="2003 v 2002")


###############
# exercises
library(rafalib)
library(Biobase)
library(GSE5859) ##Available from GitHub
data(GSE5859)

cors <- cor(exprs(e))
Pairs=which(abs(cors)>0.9999,arr.ind=TRUE)
out = Pairs[which(Pairs[,1]<Pairs[,2]),,drop=FALSE]
if(length(out[,2])>0) e=e[,-out[2]]

out <- grep("AFFX",featureNames(e))
e <- e[-out,]


# detrend
y <- exprs(e)-rowMeans(exprs(e))
dates <- pData(e)$date
eth <- pData(e)$ethnicity

annotation(e)

BiocManager::install("hgfocus.db")
library(hgfocus.db)

annot <- select(hgfocus.db, keys=featureNames(e), keytype="PROBEID",
                columns=c("CHR"))
##for genes with multiples, pick one
annot <-annot[match(featureNames(e),annot$PROBEID),]
annot$CHR <- ifelse(is.na(annot$CHR),NA,paste0("chr",annot$CHR))
##compute median expression on chromosome Y
chryexp<- colMeans(y[which(annot$CHR=="chrY"),])

mypar()
hist(chryexp)

sex <- factor(ifelse(chryexp<0,"F","M"))

# dimension reduction
s <- svd(y)
dim(s$v)

pc<-prcomp(y)

library(RColorBrewer)
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
image ( cor(y) ,col=cols,zlim=c(-1,1))

# variance explained plot for independent
y0 <- matrix( rnorm( nrow(y)*ncol(y) ) , nrow(y), ncol(y) )
d0 <- svd(y0)$d
plot(d0^2/sum(d0^2),ylim=c(0,.25))

# for our data
plot(s$d^2/sum(s$d^2))

# plot pc1 vs pc2, color by ethnicity
cols = as.numeric(eth)
mypar()
plot(s$v[,1],s$v[,2],col=cols,pch=16,
     xlab="PC1",ylab="PC2")
legend("bottomleft",levels(eth),col=seq(along=levels(eth)),pch=16)

year = factor(format(dates,"%y"))
table(year,eth)


# plot pc1 vs pc2, but now with colors=year
cols = as.numeric(year)
mypar()
plot(s$v[,1],s$v[,2],col=cols,pch=16,
     xlab="PC1",ylab="PC2")
legend("bottomleft",levels(year),col=seq(along=levels(year)),pch=16)

# boxplots of pc's

month <- format(dates,"%y%m")
length( unique(month))

variable <- as.numeric(month)
mypar(2,2)
for(i in 1:4){
  boxplot(split(s$v[,i],variable),las=2,range=0)
  stripchart(split(s$v[,i],variable),add=TRUE,vertical=TRUE,pch=1,cex=.5,col=1)
  
  
}

# linear regression (anova)
corr <- sapply(1:ncol(s$v),function(i){
  fit <- lm(s$v[,i]~as.factor(month))
  return( summary(fit)$adj.r.squared  )
})
mypar()
plot(seq(along=corr), corr, xlab="PC")

# f stat to compare within month with month to month variabiliyt
Fstats<- sapply(1:ncol(s$v),function(i){
  fit <- lm(s$v[,i]~as.factor(month))
  Fstat <- summary(aov(fit))[[1]][1,4]
  return(Fstat)
})
mypar()
plot(seq(along=Fstats),sqrt(Fstats))
p <- length(unique(month))
abline(h=sqrt(qf(0.995,p-1,ncol(s$v)-1)))

