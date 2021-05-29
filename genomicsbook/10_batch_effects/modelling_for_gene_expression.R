.libPaths('C:/Users/BrBl1834/R/win-library')

library(parathyroidSE) ##available from Bioconductor
data(parathyroidGenesSE)
se <- parathyroidGenesSE

# interpretation: each cell in parathyroidGenesSE represents #reads for a given gene,sample
# here, we plotthe average log2 of the genes across two samples
# the lower the average expression, the larger the variation across samples
x <- assay(se)[,23]
y <- assay(se)[,24]
ind=which(y>0 & x>0)##make sure no 0s due to ratio and log
plot((log2(x)+log2(y))/2,log(x/y),subset=ind)


# four samples
vars=rowVars(assay(se)[,c(2,8,16,21)]) ##we now these four are 4
means=rowMeans(assay(se)[,c(2,8,16,21)]) ##different individulsa

plot(means,vars,log="xy",subset=which(means>0&vars>0)) ##plot a subset of data
abline(0,1,col=2,lwd=2)