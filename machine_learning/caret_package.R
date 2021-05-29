.libPaths('C:/Users/BrBl1834/R/win-library')
library(rpart)
library(dslabs)
library(caret)
data("tissue_gene_expression")


# Q1 regression tree to predict tissue type
set.seed(1991, sample.kind = "Rounding") # if using R 3.6 or later
cp<-data.frame(cp=seq(0, 0.1, 0.01))
fit<-train(tissue_gene_expression$x,tissue_gene_expression$y,method = "rpart",tuneGrid = cp)
plot(fit)


# Q2 try multiple minsplit arguments (because low number of placenta tissue types)

set.seed(1991, sample.kind = "Rounding") # if using R 3.6 or later
cp<-data.frame(cp=seq(0, 0.1, 0.01))
fit_rpart<-train(tissue_gene_expression$x,tissue_gene_expression$y,method = "rpart",tuneGrid = cp,control = rpart.control(minsplit = 0))
plot(fit_rpart)
# don't do this--> model is massively overtrained because all data was used to build it
y_hat<-predict(fit_rpart,tissue_gene_expression$x)
confusionMatrix(y_hat,tissue_gene_expression$y)

fit_rpart$results$Accuracy
fit_rpart$results%>%ggplot(aes(x=cp,y=Accuracy)) + geom_line() + geom_point() + 
  geom_errorbar(aes(x=cp, ymin=Accuracy-AccuracySD, ymax=Accuracy+AccuracySD))


# Q3 plot best tree

plot(fit$finalModel)
text(fit$finalModel)

# Q4 use random forest with diff # of trees (mtry)
set.seed(1991, sample.kind = "Rounding") # if using R 3.6 or later
mtry<-data.frame(mtry=seq(50, 200, 25))
fit<-train(tissue_gene_expression$x,tissue_gene_expression$y,method = "rf",tuneGrid = mtry, nodesize=1)

fit$results

# Q5
imp <- varImp(fit)
  imp

# Q6
tree_terms <- as.character(unique(fit_rpart$finalModel$frame$var[!(fit_rpart$finalModel$frame$var == "<leaf>")]))
tree_terms



q6<-imp$importance%>%mutate(gene=rownames(imp$importance))%>%filter(gene %in% tree_terms)
q6

q6b<-imp$importance%>%mutate(gene=rownames(imp$importance))
q6b_order<-q6b$gene[order(q6b$Overall)]
q6b_order<-rev(q6b_order)
     
which(q6b_order=="CFHR4")
