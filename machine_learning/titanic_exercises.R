.libPaths('C:/Users/BrBl1834/R/win-library')
library(titanic)    # loads titanic_train data frame
library(caret)
library(tidyverse)
library(rpart)

# 3 significant digits
options(digits = 3)

# clean the data - `titanic_train` is loaded with the titanic package
titanic_clean <- titanic_train %>%
  mutate(Survived = factor(Survived),
         Embarked = factor(Embarked),
         Age = ifelse(is.na(Age), median(Age, na.rm = TRUE), Age), # NA age to median age
         FamilySize = SibSp + Parch + 1) %>%    # count family members
  select(Survived,  Sex, Pclass, Age, Fare, SibSp, Parch, FamilySize, Embarked)

# Q1 
set.seed(42,sample.kind = "Rounding")

test_index <- createDataPartition(titanic_clean$Survived,times=1,p=0.2, list=FALSE)

trainset<-titanic_clean[-test_index,]
testset<-titanic_clean[test_index,]

mean(trainset$Survived==1)

# Q2
set.seed(3,sample.kind = "Rounding")

y_hat<-sample(c(0,1),nrow(testset),replace=TRUE)

mean(y_hat==testset$Survived)

# Q3
# a

male <- trainset%>%filter(Sex=="male")
mean(male$Survived==1)

female <- trainset%>%filter(Sex=="female")
mean(female$Survived==1)


# b
y_hat_sex<-ifelse(testset$Sex=="male", 0,1)
mean(y_hat_sex==testset$Survived)

# Q4
# a survival by class
trainset%>%group_by(Pclass)%>%summarise(prob_surv=mean(Survived==1))

# b predict survival by class
y_hat_class <- ifelse(testset$Pclass==1,1,0)
mean(y_hat_class==testset$Survived)

# c both sex and class
trainset%>%group_by(Pclass,Sex)%>%summarise(prob_surv=mean(Survived==1))

# d 
probs<-trainset%>%group_by(Pclass,Sex)%>%summarise(prob_surv=mean(Survived==1))
probs<-probs%>%mutate(prob_surv=ifelse(
  prob_surv>0.5, 1, 0
))

y_hat_sexclass <- merge(testset,probs,by=c("Pclass","Sex"))
mean(y_hat_sexclass$Survived==y_hat_sexclass$prob_surv)

# Q5
# a) confusion matrix
y_hat_sex<-as.factor(y_hat_sex)
y_hat_class<-as.factor(y_hat_class)
y_hat_sexclass$prob_surv<-as.factor(y_hat_sexclass$prob_surv)

testset$Survived<-as.factor(testset$Survived)

confusionMatrix(y_hat_sex,testset$Survived)
confusionMatrix(y_hat_class,testset$Survived)
confusionMatrix(y_hat_sexclass$prob_surv,testset$Survived)

# Q6 F1 scores
print("Sex")
F_meas(y_hat_sex,testset$Survived)
print("Class")
F_meas(y_hat_class,testset$Survived)
print("sex + class")
F_meas(y_hat_sexclass$prob_surv,testset$Survived)


# Q7 LDA + QDA on fares

fit_lda<-train(Survived~Fare, data=trainset,method="lda")
y_hat_lda <- predict(fit_lda,testset) 
mean(y_hat_lda==testset$Survived)


set.seed(1, sample.kind = "Rounding")
fit_qda<-train(Survived~Fare, data=trainset,method="qda")
y_hat_qda <- predict(fit_qda,testset) 
mean(y_hat_qda==testset$Survived)

# Q8 glm
fit_glm_age<-train(Survived~Age, data=trainset,method="glm")
y_hat_glm_age <- predict(fit_glm_age,testset) 
mean(y_hat_glm_age==testset$Survived)

fit_glm<-train(Survived~ Sex + Pclass + Fare + Age, data=trainset,method="glm")
y_hat_glm <- predict(fit_glm,testset) 
mean(y_hat_glm==testset$Survived)

fit_glm_all<-train(Survived~. , data=trainset,method="glm")
y_hat_glm_all <- predict(fit_glm_all,testset) 
mean(y_hat_glm_all==testset$Survived)


# Q9 Knn

set.seed(6,sample.kind = "Rounding")
ks<-data.frame(k=seq(3, 51, 2))

fit_knn<-train(Survived~., data=trainset,method="knn",tuneGrid=ks)
plot(fit_knn)

y_hat_knn <- predict(fit_knn, testset)
mean(y_hat_knn==testset$Survived)

# Q10 crossvalidation

set.seed(8,sample.kind = "Rounding")
ks<-data.frame(k=seq(3, 51, 2))
train_control <- trainControl('cv',10,p=0.1)

fit_knn_cv<-train(Survived~., data=trainset,method="knn",tuneGrid=ks, trControl=train_control)
plot(fit_knn_cv)

y_hat_knn_cv <- predict(fit_knn_cv, testset)
mean(y_hat_knn_cv==testset$Survived)

# Q11 CART
set.seed(10,sample.kind = "Rounding")

cp <- data.frame(cp=seq(0, 0.05, 0.002))
fit_CART <- train(Survived~., data=trainset,method="rpart",tuneGrid=cp)
plot(fit_CART)

y_hat_CART <- predict(fit_CART, testset)
mean(y_hat_CART==testset$Survived)

plot(fit_CART$finalModel)
text(fit_CART$finalModel)

# Q12 random forest
set.seed(14, sample.kind = "Rounding") # if using R 3.6 or later
mtry<-data.frame(mtry=seq(1:7))
fit_rf<-train(Survived~., data=trainset,method = "rf",tuneGrid = mtry, ntree=100)

plot(fit_rf)

y_hat_rf <- predict(fit_rf, testset)
mean(y_hat_rf==testset$Survived)
y_hat_rf <- as.factor(y_hat_rf)
confusionMatrix(y_hat_rf,testset$Survived)

imp <- varImp(fit_rf)
