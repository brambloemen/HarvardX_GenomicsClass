.libPaths('C:/Users/BrBl1834/R/win-library')

library(dslabs)
library(caret)
data(mnist_27)
# set.seed(1995) # if R 3.5 or earlier
set.seed(1995, sample.kind="Rounding") # if R 3.6 or later
indexes <- createResample(mnist_27$train$y, 10)

# Q1 how often do certain indexes occur
sum(indexes[[1]]==3)
sum(indexes[[1]]==4)
sum(indexes[[1]]==7)

# Q2 how often does 3 occur in all samples
for (i in 2:10){
  test[i]<-sum(indexes[[i]]==3)
}
sum(test)

# Q3
y <- rnorm(100, 0, 1)
qnorm(0.75)
quantile(y, 0.75)

set.seed(1,sample.kind = "Rounding")
q3<-replicate(10000,{
  y <- rnorm(100, 0, 1)
  quantile(y, 0.75)
})

mean(q3)
sd(q3) #Remember that the standard error of a sample statistic is the standard deviation of its sampling distribution

# Q4

# set.seed(1) # if R 3.5 or earlier
set.seed(1, sample.kind = "Rounding") # if R 3.6 or later
y <- rnorm(100, 0, 1)

set.seed(1, sample.kind = "Rounding") # if R 3.6 or later

indexes <- createResample(y, 10)
q4<-numeric()
for(i in 1:10){
  q4[i]<-quantile(y[indexes[[i]]],0.75)
}
mean(q4)
sd(q4)

# Q5: cfr Q4, now 10000 bootstrap samples

# set.seed(1) # if R 3.5 or earlier
set.seed(1, sample.kind = "Rounding") # if R 3.6 or later
y <- rnorm(100, 0, 1)

set.seed(1, sample.kind = "Rounding") # if R 3.6 or later

indexes <- createResample(y, 10000)
q4<-numeric()
for(i in 1:10000){
  q4[i]<-quantile(y[indexes[[i]]],0.75)
}
mean(q4)
sd(q4)