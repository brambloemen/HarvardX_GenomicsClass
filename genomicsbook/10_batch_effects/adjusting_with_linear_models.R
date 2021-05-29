
# open GSE58589subset.rda before running script
library(tidyverse)
library(genefilter)
library(qvalue)
# correcting for confounding with linear models

# exercises

sex = sampleInfo$group
month = factor( format(sampleInfo$date,"%m"))
table( sampleInfo$group, month)

# t test to find differences across sex
pvals <- rowttests(geneExpression, fac = factor(sex))

# qval to account for False positives due to many t tests being performed
qvals <- qvalue(pvals$p.value)
signif_diff_expressed <- geneAnnotation[qvals$qvalues<0.1,]

# which of the sign diff expressed are on X or Y chromosomes
sign_diff_allosome <- signif_diff_expressed%>%filter(CHR %in% c("chrX","chrY"))
sign_diff_autosome <- signif_diff_expressed%>%filter(!(CHR %in% c("chrX","chrY")))

nrow(sign_diff_allosome)/nrow(signif_diff_expressed)

# which of sign diff are on Y
sign_diff_Y <- signif_diff_expressed%>%filter(CHR=="chrY")


# june vs october: batch effects?

expr_sign_diff_autosome <- geneExpression[sign_diff_autosome$PROBEID,]

pvals <- rowttests(expr_sign_diff_autosome, fac = factor(month))

mean(pvals$p.value<0.05)


#the autosomal expr. diff. are not a (always) result of the sex, 
#but also a result of processing month -> confounding factor




# correct for batch effect
X = model.matrix(~sex+month)

res <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,c(1,4)] #keep only sex effect=> month effect is removed
} ) )
res <- data.frame(res)
names(res)<-c('estimate','pvalue')
qvals <- qvalue(res$pvalue)$qvalue
index <- which(qvals<0.1)

# compare corrected with data not corrected for batch (=month) effect
signif_diff_batch_corrected <- geneAnnotation[index,]

sign_diff_corr_allo <- signif_diff_batch_corrected%>%filter(CHR %in% c("chrX","chrY"))
sign_diff_corr_auto <- signif_diff_batch_corrected%>%filter(!(CHR %in% c("chrX","chrY")))
nrow(sign_diff_corr_allo)/nrow(signif_diff_batch_corrected)

sign_diff_corr_Y<- signif_diff_batch_corrected%>%filter(CHR =="chrY")


# pvals and qvals for batch effect

res_month <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit <- lm(y~X-1)
  summary(fit)$coef[3,c(1,4)] #keep only month effect=> sex effect is removed
} ) )
res_month <- data.frame(res_month)
names(res_month)<-c('estimate','pvalue')
qvals <- qvalue(res_month$pvalue)$qvalue
index <- which(qvals<0.1)

signif_diff_batch_effect <- geneAnnotation[index,]
