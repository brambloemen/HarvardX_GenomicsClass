.libPaths('C:/Users/BrBl1834/R/win-library')
library(caret)
library(dslabs)
library(tidyverse)


# Q1

models <- c("glm", "lda", "naive_bayes", "svmLinear", "knn", "gamLoess", "multinom", "qda", "rf", "adaboost")

# set.seed(1) # if using R 3.5 or earlier
set.seed(1, sample.kind = "Rounding") # if using R 3.6 or later
data("mnist_27")

fits <- lapply(models, function(model){ 
  print(model)
  train(y ~ ., method = model, data = mnist_27$train)
}) 

names(fits) <- models

# Q2 predicions for all models

y_hat_all <- sapply(fits, FUN=predict, mnist_27$test)

# Q3
compare <- y_hat_all==mnist_27$test$y
mean(compare)

accs <- colMeans(y_hat_all == mnist_27$test$y)
accs
mean(accs)

# Q4
ensemble <- y_hat_all==7
ensemble[,1:10] <- as.numeric(ensemble[,1:10])
ensemble <- rowSums(ensemble[,1:10])
ensemble <- ifelse(ensemble>5,7,2)

mean(ensemble==mnist_27$test$y)

votes <- rowMeans(y_hat_all == "7")
y_hat <- ifelse(votes > 0.5, "7", "2")
mean(y_hat == mnist_27$test$y)

# Q5
acc<-numeric()
for(i in 1:10){
  acc[i]<-mean(compare[,i])
}

betterthanensemble<- models[which(acc>0.81)]

ind <- accs > mean(y_hat == mnist_27$test$y)
sum(ind)
models[ind]

# Q6
mean(acc)

# Q7

y_hat_best <- y_hat_all[,ind]
votes_best <- rowMeans(y_hat_best == "7")
y_hat_best <- ifelse(votes_best > 0.5, "7", "2")

mean(y_hat_best == mnist_27$test$y)
