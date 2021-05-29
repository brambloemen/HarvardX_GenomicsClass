.libPaths('C:/Users/BrBl1834/R/win-library')
library(dslabs)
library(dplyr)
library(lubridate)
library(caret)

data(iris)
iris <- iris[-which(iris$Species=='setosa'),]
y <- iris$Species

set.seed(2, sample.kind="Rounding")
# line of code
test_index<-createDataPartition(y,times=1,p=0.5,list = FALSE)
test <- iris[test_index,]
train <- iris[-test_index,]

# for different cutoffs, give accuracy ->repeat for all vars
#  for one var 

cutoff <- seq(min(train$Sepal.Length), max(train$Sepal.Length), 0.1)
accuracy <- sapply(cutoff, function(x){
  y_hat <- ifelse(train$Sepal.Length > x, "virginica", "versicolor") %>% 
    factor(levels = levels(train$Species))
    mean(y_hat == train$Species)
})


# function by var

plot_predict_species<-function(var){
cutoff <- seq(min(var), max(var), 0.1)
accuracy <- sapply(cutoff, function(x){
  y_hat <- ifelse(var > x, "virginica", "versicolor") %>% 
    factor(levels = levels(train$Species))
    mean(y_hat == train$Species)
})
plot(cutoff, accuracy)
}

# plot accuracies
plot_predict_species(train$Sepal.Length)
plot_predict_species(train$Sepal.Width)
plot_predict_species(train$Petal.Length)
plot_predict_species(train$Petal.Width)


# return vectors only
predict_species<-function(var){
  cutoff <- seq(min(var), max(var), 0.1)
  accuracy <- sapply(cutoff, function(x){
    y_hat <- ifelse(var > x, "virginica", "versicolor") %>% 
      factor(levels = levels(train$Species))
    mean(y_hat == train$Species)
  })
}

a<-predict_species(train$Sepal.Length)
b<-predict_species(train$Sepal.Width)
c<-predict_species(train$Petal.Length)
d<-predict_species(train$Petal.Width)


# now use best cutoff for best var to predict test
y_hat <- ifelse(test$Petal.Length > 4.7, "virginica", "versicolor") %>% 
  factor(levels = levels(test$Species))
mean(y_hat == test$Species)

# see if for test data, there is better predictor
predict_species<-function(var){
  cutoff <- seq(min(var), max(var), 0.1)
  accuracy <- sapply(cutoff, function(x){
    y_hat <- ifelse(var > x, "virginica", "versicolor") %>% 
      factor(levels = levels(test$Species))
    mean(y_hat == test$Species)
  })
}

a<-predict_species(test$Sepal.Length)
b<-predict_species(test$Sepal.Width)
c<-predict_species(test$Petal.Length)
d<-predict_species(test$Petal.Width)


# use best cutoffs from train for both petal length and width: 4.7 and 1.5, respectively

  y_hat <- ifelse(test$Petal.Length > 4.7 | test$Petal.Width > 1.5, "virginica", "versicolor") %>% 
    factor(levels = levels(test$Species))
  mean(y_hat == test$Species)

