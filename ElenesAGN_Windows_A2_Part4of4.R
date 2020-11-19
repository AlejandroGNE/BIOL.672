# Alejandro G.N. Elenes
# Code written in Windows 10 Home
# R packages needed to run this script: ----
library(readr)
library(class)
library(gmodels)
library(ggplot2)
library(caret)
library(GGally)
library(corrplot)
library(Rmisc)
library(e1071)
library(caret)
library(MASS)
library(dplyr)
library(randomForest)
library(ada)
# Data files necessary to run this script:
# breastcancer.csv - included with script, can also be downloaded from https://www.kaggle.com/uciml/breast-cancer-wisconsin-data

# Read in data ----
bcd <- read.csv("breastcancer.csv")

str(bcd)
bcd <- subset(bcd,select=-c(id,X))
table(bcd$diagnosis)
round(prop.table(table(bcd$diagnosis)) * 100, digits = 1)
sum(is.na(bcd))
head(bcd)

bcd_predictors_only <- subset(bcd,select=-diagnosis)
cor(bcd_predictors_only)
summary(bcd)
knn_scatter_full <- ggplot(bcd, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point()

# normalize predictors ----
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
normalized_bcd <- as.data.frame(lapply(bcd_predictors_only, normalize))
normalized_bcd_predictors_only <- normalized_bcd
normalized_bcd$diagnosis <- bcd$diagnosis

summary(normalized_bcd$radius_mean,normalized_bcd$smoothness_mean)

predictors <- as.data.frame(normalized_bcd_predictors_only)
full <- as.data.frame(normalized_bcd)
full$diagnosis <- factor(full$diagnosis)
diagnosis_only <- full$diagnosis
# corrs <- predictors[, sapply(predictors, is.numeric)] %>% na.omit() %>%
#   ggpairs(aes(col = diagnosis_only, alpha=.4))

# 5-fold cross-validation ----
set.seed(123)
Partitions <- createDataPartition(bcd$diagnosis,5,p=0.8)

# 1 of 5 ----
train1 <- normalized_bcd[Partitions$Resample1,]
validation1 <- normalized_bcd[-Partitions$Resample1,]
train <- subset(train1, select=-diagnosis)
validation <- subset(validation1, select=-diagnosis)

rndfor_train1 <- randomForest(as.factor(diagnosis)~., train1, ntree=500)
rndfor_pred1 <- predict(rndfor_train1, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm1 <- CrossTable(x = validation1$diagnosis , y = rndfor_pred1, prop.chisq = FALSE)
rndfor_accuracy1 <- ((rndfor_cm1$t[1]+rndfor_cm1$t[4])/sum(rndfor_cm1$t))*100

ada_train1 <- ada(as.factor(diagnosis)~., train1, 500)
ada_pred1 <- predict(ada_train1, validation, probability = FALSE, decision.values = TRUE)
ada_cm1 <- CrossTable(x = validation1$diagnosis , y = ada_pred1, prop.chisq = FALSE)
ada_accuracy1 <- ((ada_cm1$t[1]+ada_cm1$t[4])/sum(ada_cm1$t))*100

# 2 of 5 ----
train2 <- normalized_bcd[Partitions$Resample2,]
validation2 <- normalized_bcd[-Partitions$Resample2,]
train <- subset(train2, select=-diagnosis)
validation <- subset(validation2, select=-diagnosis)

rndfor_train2 <- randomForest(as.factor(diagnosis)~., train2, ntree=500)
rndfor_pred2 <- predict(rndfor_train2, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm2 <- CrossTable(x = validation2$diagnosis , y = rndfor_pred2, prop.chisq = FALSE)
rndfor_accuracy2 <- ((rndfor_cm2$t[1]+rndfor_cm2$t[4])/sum(rndfor_cm2$t))*100

ada_train2 <- ada(as.factor(diagnosis)~., train2, 500)
ada_pred2 <- predict(ada_train2, validation, probability = FALSE, decision.values = TRUE)
ada_cm2 <- CrossTable(x = validation2$diagnosis , y = ada_pred2, prop.chisq = FALSE)
ada_accuracy2 <- ((ada_cm2$t[1]+ada_cm2$t[4])/sum(ada_cm2$t))*100

# 3 of 5 ----
train3 <- normalized_bcd[Partitions$Resample3,]
validation3 <- normalized_bcd[-Partitions$Resample3,]
train <- subset(train3, select=-diagnosis)
validation <- subset(validation3, select=-diagnosis)

rndfor_train3 <- randomForest(as.factor(diagnosis)~., train3, ntree=500)
rndfor_pred3 <- predict(rndfor_train3, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm3 <- CrossTable(x = validation3$diagnosis , y = rndfor_pred3, prop.chisq = FALSE)
rndfor_accuracy3 <- ((rndfor_cm3$t[1]+rndfor_cm3$t[4])/sum(rndfor_cm3$t))*100

ada_train3 <- ada(as.factor(diagnosis)~., train3, 500)
ada_pred3 <- predict(ada_train3, validation, probability = FALSE, decision.values = TRUE)
ada_cm3 <- CrossTable(x = validation3$diagnosis , y = ada_pred3, prop.chisq = FALSE)
ada_accuracy3 <- ((ada_cm3$t[1]+ada_cm3$t[4])/sum(ada_cm3$t))*100

# 4 of 5 ----
train4 <- normalized_bcd[Partitions$Resample4,]
validation4 <- normalized_bcd[-Partitions$Resample4,]
train <- subset(train4, select=-diagnosis)
validation <- subset(validation4, select=-diagnosis)

rndfor_train4 <- randomForest(as.factor(diagnosis)~., train4, ntree=500)
rndfor_pred4 <- predict(rndfor_train4, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm4 <- CrossTable(x = validation4$diagnosis , y = rndfor_pred4, prop.chisq = FALSE)
rndfor_accuracy4 <- ((rndfor_cm4$t[1]+rndfor_cm4$t[4])/sum(rndfor_cm4$t))*100

ada_train4 <- ada(as.factor(diagnosis)~., train4, 500)
ada_pred4 <- predict(ada_train4, validation, probability = FALSE, decision.values = TRUE)
ada_cm4 <- CrossTable(x = validation4$diagnosis , y = ada_pred4, prop.chisq = FALSE)
ada_accuracy4 <- ((ada_cm4$t[1]+ada_cm4$t[4])/sum(ada_cm4$t))*100

# 5 of 5 ----
train5 <- normalized_bcd[Partitions$Resample5,]
validation5 <- normalized_bcd[-Partitions$Resample5,]
train <- subset(train5, select=-diagnosis)
validation <- subset(validation5, select=-diagnosis)

rndfor_train5 <- randomForest(as.factor(diagnosis)~., train5, ntree=500)
rndfor_pred5 <- predict(rndfor_train5, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm5 <- CrossTable(x = validation5$diagnosis , y = rndfor_pred5, prop.chisq = FALSE)
rndfor_accuracy5 <- ((rndfor_cm5$t[1]+rndfor_cm5$t[4])/sum(rndfor_cm5$t))*100

ada_train5 <- ada(as.factor(diagnosis)~., train5, 500)
ada_pred5 <- predict(ada_train5, validation, probability = FALSE, decision.values = TRUE)
ada_cm5 <- CrossTable(x = validation5$diagnosis , y = ada_pred5, prop.chisq = FALSE)
ada_accuracy5 <- ((ada_cm5$t[1]+ada_cm5$t[4])/sum(ada_cm5$t))*100

# Accuracies ----
rndfor_accuracies <- c(rndfor_accuracy1,rndfor_accuracy2,rndfor_accuracy3,rndfor_accuracy4,rndfor_accuracy5)
ada_accuracies <- c(ada_accuracy1,ada_accuracy2,ada_accuracy3,ada_accuracy4,ada_accuracy5)
rndfor_accuracies
ada_accuracies

# Validation scatterplots ----
scatter_validation1 <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 1
scatter_validation2 <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 2
scatter_validation3 <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 3
scatter_validation4 <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 4
scatter_validation5 <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 5

# Prediction scatterplots rndfor
rndfor_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(rndfor_pred1))) + geom_point()
rndfor_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(rndfor_pred2))) + geom_point()
rndfor_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(rndfor_pred3))) + geom_point()
rndfor_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(rndfor_pred4))) + geom_point()
rndfor_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(rndfor_pred5))) + geom_point()

# Prediction scatterplots ada
ada_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(ada_pred1))) + geom_point()
ada_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(ada_pred2))) + geom_point()
ada_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(ada_pred3))) + geom_point()
ada_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(ada_pred4))) + geom_point()
ada_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(ada_pred5))) + geom_point()

# export plots ----

pdf(file = 'ElenesAGN_A2_Problem4_GraphicOutputs.pdf')

multiplot(scatter_validation1,rndfor_scatter_1_pred)
multiplot(scatter_validation2,rndfor_scatter_2_pred)
multiplot(scatter_validation3,rndfor_scatter_3_pred)
multiplot(scatter_validation4,rndfor_scatter_4_pred)
multiplot(scatter_validation5,rndfor_scatter_5_pred)

multiplot(scatter_validation1,ada_scatter_1_pred)
multiplot(scatter_validation2,ada_scatter_2_pred)
multiplot(scatter_validation3,ada_scatter_3_pred)
multiplot(scatter_validation4,ada_scatter_4_pred)
multiplot(scatter_validation5,ada_scatter_5_pred)

dev.off()

# png(filename="Pairs_correlograms.png",width=3840,height=2160)
# multiplot(corrs)
# dev.off()

# confusion matrices ----

rndfor_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(rndfor_pred1))
rndfor_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(rndfor_pred2))
rndfor_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(rndfor_pred3))
rndfor_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(rndfor_pred4))
rndfor_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(rndfor_pred5))

ada_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(ada_pred1))
ada_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(ada_pred2))
ada_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(ada_pred3))
ada_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(ada_pred4))
ada_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(ada_pred5))

# export text outputs ----

sink(file = 'ElenesAGN_A2_Problem4_TextOutputs.txt')

writeLines(" \n\n Confusion matrices for 5-fold cross validation of random forest \n")
print(rndfor_confusionmatrix1)
print(rndfor_confusionmatrix2)
print(rndfor_confusionmatrix3)
print(rndfor_confusionmatrix4)
print(rndfor_confusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of boosted gradient \n")
print(ada_confusionmatrix1)
print(ada_confusionmatrix2)
print(ada_confusionmatrix3)
print(ada_confusionmatrix4)
print(ada_confusionmatrix5)

writeLines(" \n\n Accuracies for random forest and booster gradient \n")
print(rndfor_accuracies)
print(ada_accuracies)

writeLines(" \n\n Mean accuracies for random forest and booster gradient \n")
print(mean(rndfor_accuracies))
print(mean(ada_accuracies))

writeLines(" \n\n Boosted gradient performed better than random forest in this case. But none of these overperformed the predictions using SVMs. \n")

sink()

