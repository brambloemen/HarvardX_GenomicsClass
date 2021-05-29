.libPaths('C:/Users/BrBl1834/R/win-library')
library(caret)
library(tidyverse)
library(dslabs)


# Q1
# my answer: contains errors
data(heights)

# recoding not necessary
# heights$sex<-as.factor(recode(heights$sex, Male=1, Female=0))

set.seed(1,sample.kind="Rounding")

test_index<-createDataPartition(heights$sex, times=1, p=0.5, list=FALSE)

train<-heights[-test_index,]
test<-heights[test_index,]

k<-seq(1, 101, 3)

generate_knnfits<-function(k,data=heights){
  kfit<-knn3(sex~height,train,k=k)
  sex_hat<-predict(kfit, test, type="class")
  F_meas(data = sex_hat, reference = test$sex)
}
knns<-sapply(k, generate_knnfits)

kfit<-knn3(sex~height,heights,1)
sex_hat<-predict(kfit, heights, type="class")

# solution
########################
library(dslabs)
library(tidyverse)
library(caret)
data("heights")

# set.seed(1) # if using R 3.5 or earlier
set.seed(1, sample.kind = "Rounding") # if using R 3.6 or later
test_index <- createDataPartition(heights$sex, times = 1, p = 0.5, list = FALSE)
test_set <- heights[test_index, ]
train_set <- heights[-test_index, ]     

ks <- seq(1, 101, 3)
F_1 <- sapply(ks, function(k){
  fit <- knn3(sex ~ height, data = train_set, k = k)
  y_hat <- predict(fit, test_set, type = "class") %>% 
    factor(levels = levels(train_set$sex))
  F_meas(data = y_hat, reference = test_set$sex)
})
plot(ks, F_1)
max(F_1)
ks[which.max(F_1)]
########################

# Q2

library(dslabs)
library(caret)
data("tissue_gene_expression")

set.seed(1, sample.kind = "Rounding") # if using R 3.6 or later
test_index <- createDataPartition(tissue_gene_expression$y, times = 1, p = 0.5, list = FALSE)
test_set <- list(y=tissue_gene_expression$y[test_index], x=tissue_gene_expression$x[test_index, ])
train_set <- list(y=tissue_gene_expression$y[-test_index], x=tissue_gene_expression$x[-test_index, ])  

k <- seq(1, 11, 2)

knn_accs<-sapply(k,function(k){
  fit<-knn3(train_set$x, train_set$y,k=k)
  y_hat<-predict(fit, test_set$x,type="class")
  mean(y_hat==test_set$y)
})

plot(k, knn_accs)

knn_accs[which(k==1)]
knn_accs[which(k==3)]
knn_accs[which(k==5)]
knn_accs[which(k==7)]
knn_accs[which(k==9)]
knn_accs[which(k==11)]
