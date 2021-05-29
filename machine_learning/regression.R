.libPaths('C:/Users/BrBl1834/R/win-library')
library(caret)
library(tidyverse)

# Q1

# set.seed(2) #if you are using R 3.5 or earlier
set.seed(2, sample.kind="Rounding") #if you are using R 3.6 or later
make_data <- function(n = 1000, p = 0.5, 
                      mu_0 = 0, mu_1 = 2, 
                      sigma_0 = 1,  sigma_1 = 1){
  
  y <- rbinom(n, 1, p)
  f_0 <- rnorm(n, mu_0, sigma_0)
  f_1 <- rnorm(n, mu_1, sigma_1)
  x <- ifelse(y == 1, f_1, f_0)
  
  test_index <- createDataPartition(y, times = 1, p = 0.5, list = FALSE)
  
  list(train = data.frame(x = x, y = as.factor(y)) %>% slice(-test_index),
       test = data.frame(x = x, y = as.factor(y)) %>% slice(test_index))
}
dat <- make_data()


set.seed(1, sample.kind="Rounding")

mu_1 <- seq(0, 3, len=25)
my_dat<-list()

for (i in 1:25){
my_dat[[i]] <- make_data(mu_1=mu_1[i])
}

accuracies<-numeric()
for (i in 1:25){
glm_log_fit <- glm(y~x,data = my_dat[[i]]$train,family = "binomial")
y_hat <- ifelse(
  predict.glm(glm_log_fit,newdata = my_dat[[i]]$test, type = "response")>0.5,
  1,0
)
acc<- mean(y_hat==my_dat[[i]]$test$y)
accuracies[i]<-acc
}


#smoothing
# Q1
library(tidyverse)
library(lubridate)
library(purrr)
library(pdftools)

fn <- system.file("extdata", "RD-Mortality-Report_2015-18-180531.pdf", package="dslabs")
dat <- map_df(str_split(pdf_text(fn), "\n"), function(s){
  s <- str_trim(s)
  header_index <- str_which(s, "2015")[1]
  tmp <- str_split(s[header_index], "\\s+", simplify = TRUE)
  month <- tmp[1]
  header <- tmp[-1]
  tail_index  <- str_which(s, "Total")
  n <- str_count(s, "\\d+")
  out <- c(1:header_index, which(n==1), which(n>=28), tail_index:length(s))
  s[-out] %>%
    str_remove_all("[^\\d\\s]") %>%
    str_trim() %>%
    str_split_fixed("\\s+", n = 6) %>%
    .[,1:5] %>%
    as_data_frame() %>% 
    setNames(c("day", header)) %>%
    mutate(month = month,
           day = as.numeric(day)) %>%
    gather(year, deaths, -c(day, month)) %>%
    mutate(deaths = as.numeric(deaths))
}) %>%
  mutate(month = recode(month, "JAN" = 1, "FEB" = 2, "MAR" = 3, "APR" = 4, "MAY" = 5, "JUN" = 6, 
                        "JUL" = 7, "AGO" = 8, "SEP" = 9, "OCT" = 10, "NOV" = 11, "DEC" = 12)) %>%
  mutate(date = make_date(year, month, day)) %>%
  dplyr::filter(date <= "2018-05-01")

total_days <- diff(as.numeric(range(dat$date)))
span<-61/total_days
fit <- loess(deaths ~ as.numeric(date), degree=1, span = span, data=dat)
dat %>% ggplot(aes(date, deaths)) +
  geom_point() +
  geom_smooth(color="red", span = span, method = "loess", method.args = list(degree=1))


# Q2
dat %>% 
  mutate(smooth = predict(fit, as.numeric(date)), day = yday(date), year = as.character(year(date))) %>%
  ggplot(aes(day, smooth, col = year)) +
  geom_line(lwd = 2)

# Q3
library(broom)
library(dslabs)
data(mnist_27)
mnist_27$train %>% glm(y ~ x_2, family = "binomial", data = .) %>% tidy()
qplot(x_2, y, data = mnist_27$train)
## my solution
fit <- loess(as.numeric(y) ~ x_2, degree=1, span = 0.05, data=mnist_27$train)
plot(mnist_27$train$x_2,predict(fit,mnist_27$train$x_2))

## suggested method
mnist_27$train %>% 
  mutate(y = ifelse(y=="7", 1, 0)) %>%
  ggplot(aes(x_2, y)) + 
  geom_smooth(method = "loess")

data("mnist")
