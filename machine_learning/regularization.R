.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
options(digits=7)

# set.seed(1986) # if using R 3.5 or earlier
set.seed(1986, sample.kind="Rounding") # if using R 3.6 or later
n <- round(2^rnorm(1000, 8, 1))

# set.seed(1) # if using R 3.5 or earlier
set.seed(1, sample.kind="Rounding") # if using R 3.6 or later
mu <- round(80 + 2*rt(1000, 5))
range(mu)
schools <- data.frame(id = paste("PS",1:1000),
                      size = n,
                      quality = mu,
                      rank = rank(-mu))

schools %>% top_n(10, quality) %>% arrange(desc(quality))

# set.seed(1) # if using R 3.5 or earlier
set.seed(1, sample.kind="Rounding") # if using R 3.6 or later
mu <- round(80 + 2*rt(1000, 5))

scores <- sapply(1:nrow(schools), function(i){
  scores <- rnorm(schools$size[i], schools$quality[i], 30)
  scores
})
schools <- schools %>% mutate(score = sapply(scores, mean))

# Q1 top schools
top10 <- head(schools[rev(order(schools$score)),],10)
top10

# Q2
median(schools$size)
median(top10$size)

# Q3
worst10 <- tail(schools[rev(order(schools$score)),],10)
median(worst10$size)

# Q4
top10qual<-head(schools[rev(order(schools$quality)),],10)
ggplot(schools,aes(size, score)) + geom_point() + geom_smooth() + geom_point(data=top10qual,aes(size, score,col="red"))

# Q5 regularization
overall <- mean(sapply(scores, mean))

alpha <- 25
score_reg <- sapply(scores, function(x){
  overall + sum(x-overall)/(length(x)+alpha)
  })
schools %>% mutate(score_reg = score_reg) %>%
  top_n(10, score_reg) %>% arrange(desc(score_reg))



# Q6 better alpha
alphas <- seq(10,250,1)

# my code: not correct -> using sapply with function(i)-> i is not the index of "alphas", 
# but the value for each element of alpha -->first i is 10 --> alpha=alphas[10]=19, and not 10 as I thought it would
RMSEs <- sapply(alphas, function(i){
  alpha<-alphas[i]
  print(alpha)
  score_reg <- sapply(scores, function(x){
  overall + sum(x-overall)/(length(x)+alpha)
  })
  RMSE <- sqrt(mean((score_reg-schools$quality)^2))
  RMSE
})
  
alphas[which.min(RMSEs)]


# solution
alphas <- seq(10,250)
rmse <- sapply(alphas, function(alpha){
  score_reg <- sapply(scores, function(x) overall+sum(x-overall)/(length(x)+alpha))
  sqrt(mean((score_reg - schools$quality)^2))
})
plot(alphas, rmse)
alphas[which.min(rmse)]


# Q7 apply best alpha
overall <- mean(sapply(scores, mean))

alpha <- alphas[which.min(rmse)]
score_reg <- sapply(scores, function(x){
  overall + sum(x-overall)/(length(x)+alpha)
})
schools %>% mutate(score_reg = score_reg) %>%
  top_n(10, score_reg) %>% arrange(desc(score_reg))

# Q8: cfr Q6, but without centering the values around 0 (i.e. without substracting overall mean)
alphas <- seq(10,250)
rmse <- sapply(alphas, function(alpha){
  score_reg <- sapply(scores, function(x) sum(x)/(length(x)+alpha))
  sqrt(mean((score_reg - schools$quality)^2))
})
plot(alphas, rmse)
alphas[which.min(rmse)]
