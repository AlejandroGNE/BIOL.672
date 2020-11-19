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
library(neuralnet)

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

nn_train1 <- neuralnet(as.factor(diagnosis)~., train1, hidden=3, linear.output = FALSE)
nn_pred1 <- predict(nn_train1, validation, rep=1, all.units=FALSE)
nn_cm1 <- CrossTable(x = validation1$diagnosis , y = max.col(nn_pred1), prop.chisq = FALSE)
nn_accuracy1 <- ((nn_cm1$t[1]+nn_cm1$t[4])/sum(nn_cm1$t))*100
print(plot(nn_train1))

# 2 of 5 ----
train2 <- normalized_bcd[Partitions$Resample2,]
validation2 <- normalized_bcd[-Partitions$Resample2,]
train <- subset(train2, select=-diagnosis)
validation <- subset(validation2, select=-diagnosis)

nn_train2 <- neuralnet(as.factor(diagnosis)~., train2, hidden=3, linear.output = FALSE)
nn_pred2 <- predict(nn_train2, validation, rep=1, all.units=FALSE)
nn_cm2 <- CrossTable(x = validation2$diagnosis , y = max.col(nn_pred2), prop.chisq = FALSE)
nn_accuracy2 <- ((nn_cm2$t[1]+nn_cm2$t[4])/sum(nn_cm2$t))*100
print(plot(nn_train2))

# 3 of 5 ----
train3 <- normalized_bcd[Partitions$Resample3,]
validation3 <- normalized_bcd[-Partitions$Resample3,]
train <- subset(train3, select=-diagnosis)
validation <- subset(validation3, select=-diagnosis)

nn_train3 <- neuralnet(as.factor(diagnosis)~., train3, hidden=3, linear.output = FALSE)
nn_pred3 <- predict(nn_train3, validation, rep=1, all.units=FALSE)
nn_cm3 <- CrossTable(x = validation3$diagnosis , y = max.col(nn_pred3), prop.chisq = FALSE)
nn_accuracy3 <- ((nn_cm3$t[1]+nn_cm3$t[4])/sum(nn_cm3$t))*100
print(plot(nn_train3))

# 4 of 5 ----
train4 <- normalized_bcd[Partitions$Resample4,]
validation4 <- normalized_bcd[-Partitions$Resample4,]
train <- subset(train4, select=-diagnosis)
validation <- subset(validation4, select=-diagnosis)

nn_train4 <- neuralnet(as.factor(diagnosis)~., train4, hidden=3, linear.output = FALSE)
nn_pred4 <- predict(nn_train4, validation, rep=1, all.units=FALSE)
nn_cm4 <- CrossTable(x = validation4$diagnosis , y = max.col(nn_pred4), prop.chisq = FALSE)
nn_accuracy4 <- ((nn_cm4$t[1]+nn_cm4$t[4])/sum(nn_cm4$t))*100
print(plot(nn_train4))

# 5 of 5 ----
train5 <- normalized_bcd[Partitions$Resample5,]
validation5 <- normalized_bcd[-Partitions$Resample5,]
train <- subset(train5, select=-diagnosis)
validation <- subset(validation5, select=-diagnosis)

nn_train5 <- neuralnet(as.factor(diagnosis)~., train5, hidden=3, linear.output = FALSE)
nn_pred5 <- predict(nn_train5, validation, rep=1, all.units=FALSE)
nn_cm5 <- CrossTable(x = validation5$diagnosis , y = max.col(nn_pred5), prop.chisq = FALSE)
nn_accuracy5 <- ((nn_cm5$t[1]+nn_cm5$t[4])/sum(nn_cm5$t))*100
print(plot(nn_train5))

# Accuracies ----

nn_accuracies <- c(nn_accuracy1,nn_accuracy2,nn_accuracy3,nn_accuracy4,nn_accuracy5)
nn_accuracies

# Validation scatterplots ----
scatter_validation1 <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 1
scatter_validation2 <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 2
scatter_validation3 <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 3
scatter_validation4 <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 4
scatter_validation5 <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 5

# Prediction scatterplots nn
nn_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred1)))) + geom_point()
nn_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred2)))) + geom_point()
nn_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred3)))) + geom_point()
nn_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred4)))) + geom_point()
nn_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred5)))) + geom_point()

# export plots ----

pdf(file = 'ElenesAGN_A2_Problem3_GraphicOutputs.pdf')

multiplot(scatter_validation1,nn_scatter_1_pred)
multiplot(scatter_validation2,nn_scatter_2_pred)
multiplot(scatter_validation3,nn_scatter_3_pred)
multiplot(scatter_validation4,nn_scatter_4_pred)
multiplot(scatter_validation5,nn_scatter_5_pred)
dev.off()

# png(filename="Pairs_correlograms.png",width=3840,height=2160)
# multiplot(corrs)
# dev.off()

# confusion matrices ----

nn_confusionmatrix1 <- confusionMatrix(as.factor(as.integer(as.factor(validation1$diagnosis))),as.factor(max.col(nn_pred1)))
nn_confusionmatrix2 <- confusionMatrix(as.factor(as.integer(as.factor(validation2$diagnosis))),as.factor(max.col(nn_pred2)))
nn_confusionmatrix3 <- confusionMatrix(as.factor(as.integer(as.factor(validation3$diagnosis))),as.factor(max.col(nn_pred3)))
nn_confusionmatrix4 <- confusionMatrix(as.factor(as.integer(as.factor(validation4$diagnosis))),as.factor(max.col(nn_pred4)))
nn_confusionmatrix5 <- confusionMatrix(as.factor(as.integer(as.factor(validation5$diagnosis))),as.factor(max.col(nn_pred5)))

# export text outputs ----

sink(file = 'ElenesAGN_A2_Problem3_TextOutputs.txt')

writeLines(" \n\n Confusion matrices for 5-fold cross validation of neural net model \n")
print(nn_confusionmatrix1)
print(nn_confusionmatrix2)
print(nn_confusionmatrix3)
print(nn_confusionmatrix4)
print(nn_confusionmatrix5)

writeLines(" \n\n Mean accuracy for neural net \n")
print(mean(nn_accuracies))

writeLines(" \n\n Prediction power did not improve with neural net model. SVMs predicted better. \n")

sink()

