.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)

# set.seed(1) # if using R 3.5 or earlier
set.seed(1, sample.kind = "Rounding") # if using R 3.6 or later
disease <- sample(c(0,1), size=1e6, replace=TRUE, prob=c(0.98,0.02))
test <- rep(NA, 1e6)
test[disease==0] <- sample(c(0,1), size=sum(disease==0), replace=TRUE, prob=c(0.90,0.10))
test[disease==1] <- sample(c(0,1), size=sum(disease==1), replace=TRUE, prob=c(0.15, 0.85))

# prob of positive test
p_pos<-mean(test==1)

# prob of having disease when negative test
mean(disease[test==0]==1)

# prob of having disease when test positive
mean(disease[test==1]==1)

# prev among pos tested vs total prev
mean(disease[test==1]==1)/0.02


# part 2

library(dslabs)
data("heights")

# Q6
heights %>%
  mutate(height=round(height))%>%
  group_by(height)%>%
  summarize(p=mean(sex=="Male"))%>%
qplot(height, p, data =.)

# Q7
ps <- seq(0, 1, 0.1)
heights %>% 
  mutate(g=cut(height, quantile(height,ps), include.lowest = TRUE))%>%
  group_by(g) %>%
  summarize(p = mean(sex == "Male"), height = mean(height)) %>%
  qplot(height, p, data =.)

# Q8
library(MASS)
Sigma <- 9*matrix(c(1,0.5,0.5,1), 2, 2)
dat <- MASS::mvrnorm(n = 10000, c(69, 69), Sigma) %>%
  data.frame() %>% setNames(c("x", "y"))
plot(dat)

ps <- seq(0, 1, 0.1)
dat %>% 
  mutate(g = cut(x, quantile(x, ps), include.lowest = TRUE)) %>%
  group_by(g) %>%
  summarize(y = mean(y), x = mean(x)) %>%
  qplot(x, y, data =.)