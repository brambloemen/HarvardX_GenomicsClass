.libPaths('C:/Users/BrBl1834/R/win-library')


library(lubridate)
library(tidyverse)
library(dslabs)
data("movielens")

# Q1

ratings_by_year<-movielens%>%count(title,movieId, year, rating)%>%filter(year>1980)
ratings_by_year$n<-sqrt(ratings_by_year$n)
ratings_by_year$year<-as.factor(ratings_by_year$year)

ratings_by_year%>%ggplot(aes(x=year,y=n)) +
  geom_boxplot()

# Q2
ratings_by_year<-movielens%>%count(title,movieId, year)%>%filter(year>=1993)
ratings_by_year<-ratings_by_year%>%mutate(avg_yearly_ratings=n/(2018-year))
ratings_by_year <- ratings_by_year[rev(order(ratings_by_year$n)),]

# Q3
ratings_by_year<-merge(ratings_by_year, movielens, by=c("movieId","title","year"),all.x = TRUE,all.y = TRUE)

ratings_by_year<-ratings_by_year%>%group_by(movieId,title, avg_yearly_ratings)%>%summarize(avg_rating=mean(rating))


ggplot(ratings_by_year, aes(x=avg_yearly_ratings, y=avg_rating)) + geom_line()

movielens %>% 
  filter(year >= 1993) %>%
  group_by(movieId) %>%
  summarize(n = n(), years = 2018 - first(year),
            title = title[1],
            rating = mean(rating)) %>%
  mutate(rate = n/years) %>%
  ggplot(aes(rate, rating)) +
  geom_point() +
  geom_smooth()

# Q5
movielens <- mutate(movielens, date = as_datetime(timestamp))

# Q6
movielens$date<-as.Date(movielens$date)
 movielens%>%
  mutate(week=round_date(date, unit = "week"))%>%
  group_by(week)%>%
  summarize(week_mean=mean(rating))%>%
  ggplot(aes(week,week_mean))+
  geom_point() +
  geom_smooth()
 
 movielens %>% mutate(date = round_date(date, unit = "week")) %>%
   group_by(date) %>%
   summarize(rating = mean(rating)) %>%
   ggplot(aes(date, rating)) +
   geom_point() +
   geom_smooth()

 # Q7
 
 movielens%>%
   group_by(genres)%>%
   summarize(n=n(), rating=mean(rating), sd=sd(rating))%>%filter(n>1000)%>%
   ggplot(aes(genres,rating)) +
   geom_bar(stat="identity")
 