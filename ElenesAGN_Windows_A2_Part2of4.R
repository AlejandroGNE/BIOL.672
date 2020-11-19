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
library(mda)
library(dplyr)
library(kernlab)

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
corrs <- predictors[, sapply(predictors, is.numeric)] %>% na.omit() %>%
  ggpairs(aes(col = diagnosis_only, alpha=.4))

# 5-fold cross-validation ----
set.seed(123)
Partitions <- createDataPartition(bcd$diagnosis,5,p=0.8)

# 1 of 5 ----
train1 <- normalized_bcd[Partitions$Resample1,]
validation1 <- normalized_bcd[-Partitions$Resample1,]
train <- subset(train1, select=-diagnosis)
validation <- subset(validation1, select=-diagnosis)

svme_train1 <-  svm(as.factor(diagnosis)~., train1, probability = TRUE, type = "C-classification", kernel = "linear")
svme_pred1 <- predict(svme_train1, validation, probability = FALSE, decision.values = TRUE)
svme_cm1 =CrossTable(x = validation1$diagnosis , y = svme_pred1, prop.chisq = FALSE)
svme_accuracy1 <- ((svme_cm1$t[1]+svme_cm1$t[4])/sum(svme_cm1$t))*100

svmp_train1 <- svm(as.factor(diagnosis)~., train1, probability = TRUE, type = "C-classification", kernel = "polynomial", degree = 3)
svmp_pred1 <- predict(svmp_train1, validation, probability = FALSE, decision.values = TRUE)
svmp_cm1 =CrossTable(x = validation1$diagnosis , y = svmp_pred1, prop.chisq = FALSE)
svmp_accuracy1 <- ((svmp_cm1$t[1]+svmp_cm1$t[4])/sum(svmp_cm1$t))*100

svmr_train1 <- svm(as.factor(diagnosis)~., train1, probability = TRUE, type = "C-classification", kernel = "radial", gamma = 0.1)
svmr_pred1 <- predict(svmr_train1, validation, probability = FALSE, decision.values = TRUE)
svmr_cm1 =CrossTable(x = validation1$diagnosis , y = svmr_pred1, prop.chisq = FALSE)
svmr_accuracy1 <- ((svmr_cm1$t[1]+svmr_cm1$t[4])/sum(svmr_cm1$t))*100

svmke_train1 <- ksvm(as.matrix(train), as.factor(train1$diagnosis), kernel='vanilladot')
svmke_pred1 <- predict(svmke_train1, validation, type='response')
svmke_cm1 =CrossTable(x = validation1$diagnosis , y = svmke_pred1, prop.chisq = FALSE)
svmke_accuracy1 <- ((svmke_cm1$t[1]+svmke_cm1$t[4])/sum(svmke_cm1$t))*100

svmkp_train1 <- ksvm(as.matrix(train), as.factor(train1$diagnosis), kernel='polydot')
svmkp_pred1 <- predict(svmkp_train1, validation, type='response')
svmkp_cm1 =CrossTable(x = validation1$diagnosis , y = svmkp_pred1, prop.chisq = FALSE)
svmkp_accuracy1 <- ((svmkp_cm1$t[1]+svmkp_cm1$t[4])/sum(svmkp_cm1$t))*100

svmkr_train1 <- ksvm(as.matrix(train), as.factor(train1$diagnosis), kernel='rbfdot')
svmkr_pred1 <- predict(svmkr_train1, validation, type='response')
svmkr_cm1 =CrossTable(x = validation1$diagnosis , y = svmkr_pred1, prop.chisq = FALSE)
svmkr_accuracy1 <- ((svmkr_cm1$t[1]+svmkr_cm1$t[4])/sum(svmkr_cm1$t))*100

# 2 of 5 ----
train2 <- normalized_bcd[Partitions$Resample2,]
validation2 <- normalized_bcd[-Partitions$Resample2,]
train <- subset(train2, select=-diagnosis)
validation <- subset(validation2, select=-diagnosis)

svme_train2 <-  svm(as.factor(diagnosis)~., train2, probability = TRUE, type = "C-classification", kernel = "linear")
svme_pred2 <- predict(svme_train2, validation, probability = FALSE, decision.values = TRUE)
svme_cm2 =CrossTable(x = validation2$diagnosis , y = svme_pred2, prop.chisq = FALSE)
svme_accuracy2 <- ((svme_cm2$t[1]+svme_cm2$t[4])/sum(svme_cm2$t))*100

svmp_train2 <- svm(as.factor(diagnosis)~., train2, probability = TRUE, type = "C-classification", kernel = "polynomial", degree = 3)
svmp_pred2 <- predict(svmp_train2, validation, probability = FALSE, decision.values = TRUE)
svmp_cm2 =CrossTable(x = validation2$diagnosis , y = svmp_pred2, prop.chisq = FALSE)
svmp_accuracy2 <- ((svmp_cm2$t[1]+svmp_cm2$t[4])/sum(svmp_cm2$t))*100

svmr_train2 <- svm(as.factor(diagnosis)~., train2, probability = TRUE, type = "C-classification", kernel = "radial", gamma = 0.1)
svmr_pred2 <- predict(svmr_train2, validation, probability = FALSE, decision.values = TRUE)
svmr_cm2 =CrossTable(x = validation2$diagnosis , y = svmr_pred2, prop.chisq = FALSE)
svmr_accuracy2 <- ((svmr_cm2$t[1]+svmr_cm2$t[4])/sum(svmr_cm2$t))*100

svmke_train2 <- ksvm(as.matrix(train), as.factor(train2$diagnosis), kernel='vanilladot')
svmke_pred2 <- predict(svmke_train2, validation, type='response')
svmke_cm2 =CrossTable(x = validation2$diagnosis , y = svmke_pred2, prop.chisq = FALSE)
svmke_accuracy2 <- ((svmke_cm2$t[1]+svmke_cm2$t[4])/sum(svmke_cm2$t))*100

svmkp_train2 <- ksvm(as.matrix(train), as.factor(train2$diagnosis), kernel='polydot')
svmkp_pred2 <- predict(svmkp_train2, validation, type='response')
svmkp_cm2 =CrossTable(x = validation2$diagnosis , y = svmkp_pred2, prop.chisq = FALSE)
svmkp_accuracy2 <- ((svmkp_cm2$t[1]+svmkp_cm2$t[4])/sum(svmkp_cm2$t))*100

svmkr_train2 <- ksvm(as.matrix(train), as.factor(train2$diagnosis), kernel='rbfdot')
svmkr_pred2 <- predict(svmkr_train2, validation, type='response')
svmkr_cm2 =CrossTable(x = validation2$diagnosis , y = svmkr_pred2, prop.chisq = FALSE)
svmkr_accuracy2 <- ((svmkr_cm2$t[1]+svmkr_cm2$t[4])/sum(svmkr_cm2$t))*100

# 3 of 5 ----
train3 <- normalized_bcd[Partitions$Resample3,]
validation3 <- normalized_bcd[-Partitions$Resample3,]
train <- subset(train3, select=-diagnosis)
validation <- subset(validation3, select=-diagnosis)

svme_train3 <-  svm(as.factor(diagnosis)~., train3, probability = TRUE, type = "C-classification", kernel = "linear")
svme_pred3 <- predict(svme_train3, validation, probability = FALSE, decision.values = TRUE)
svme_cm3 =CrossTable(x = validation3$diagnosis , y = svme_pred3, prop.chisq = FALSE)
svme_accuracy3 <- ((svme_cm3$t[1]+svme_cm3$t[4])/sum(svme_cm3$t))*100

svmp_train3 <- svm(as.factor(diagnosis)~., train3, probability = TRUE, type = "C-classification", kernel = "polynomial", degree = 3)
svmp_pred3 <- predict(svmp_train3, validation, probability = FALSE, decision.values = TRUE)
svmp_cm3 =CrossTable(x = validation3$diagnosis , y = svmp_pred3, prop.chisq = FALSE)
svmp_accuracy3 <- ((svmp_cm3$t[1]+svmp_cm3$t[4])/sum(svmp_cm3$t))*100

svmr_train3 <- svm(as.factor(diagnosis)~., train3, probability = TRUE, type = "C-classification", kernel = "radial", gamma = 0.1)
svmr_pred3 <- predict(svmr_train3, validation, probability = FALSE, decision.values = TRUE)
svmr_cm3 =CrossTable(x = validation3$diagnosis , y = svmr_pred3, prop.chisq = FALSE)
svmr_accuracy3 <- ((svmr_cm3$t[1]+svmr_cm3$t[4])/sum(svmr_cm3$t))*100

svmke_train3 <- ksvm(as.matrix(train), as.factor(train3$diagnosis), kernel='vanilladot')
svmke_pred3 <- predict(svmke_train3, validation, type='response')
svmke_cm3 =CrossTable(x = validation3$diagnosis , y = svmke_pred3, prop.chisq = FALSE)
svmke_accuracy3 <- ((svmke_cm3$t[1]+svmke_cm3$t[4])/sum(svmke_cm3$t))*100

svmkp_train3 <- ksvm(as.matrix(train), as.factor(train3$diagnosis), kernel='polydot')
svmkp_pred3 <- predict(svmkp_train3, validation, type='response')
svmkp_cm3 =CrossTable(x = validation3$diagnosis , y = svmkp_pred3, prop.chisq = FALSE)
svmkp_accuracy3 <- ((svmkp_cm3$t[1]+svmkp_cm3$t[4])/sum(svmkp_cm3$t))*100

svmkr_train3 <- ksvm(as.matrix(train), as.factor(train3$diagnosis), kernel='rbfdot')
svmkr_pred3 <- predict(svmkr_train3, validation, type='response')
svmkr_cm3 =CrossTable(x = validation3$diagnosis , y = svmkr_pred3, prop.chisq = FALSE)
svmkr_accuracy3 <- ((svmkr_cm3$t[1]+svmkr_cm3$t[4])/sum(svmkr_cm3$t))*100

# 4 of 5 ----
train4 <- normalized_bcd[Partitions$Resample4,]
validation4 <- normalized_bcd[-Partitions$Resample4,]
train <- subset(train4, select=-diagnosis)
validation <- subset(validation4, select=-diagnosis)

svme_train4 <-  svm(as.factor(diagnosis)~., train4, probability = TRUE, type = "C-classification", kernel = "linear")
svme_pred4 <- predict(svme_train4, validation, probability = FALSE, decision.values = TRUE)
svme_cm4 =CrossTable(x = validation4$diagnosis , y = svme_pred4, prop.chisq = FALSE)
svme_accuracy4 <- ((svme_cm4$t[1]+svme_cm4$t[4])/sum(svme_cm4$t))*100

svmp_train4 <- svm(as.factor(diagnosis)~., train4, probability = TRUE, type = "C-classification", kernel = "polynomial", degree = 3)
svmp_pred4 <- predict(svmp_train4, validation, probability = FALSE, decision.values = TRUE)
svmp_cm4 =CrossTable(x = validation4$diagnosis , y = svmp_pred4, prop.chisq = FALSE)
svmp_accuracy4 <- ((svmp_cm4$t[1]+svmp_cm4$t[4])/sum(svmp_cm4$t))*100

svmr_train4 <- svm(as.factor(diagnosis)~., train4, probability = TRUE, type = "C-classification", kernel = "radial", gamma = 0.1)
svmr_pred4 <- predict(svmr_train4, validation, probability = FALSE, decision.values = TRUE)
svmr_cm4 =CrossTable(x = validation4$diagnosis , y = svmr_pred4, prop.chisq = FALSE)
svmr_accuracy4 <- ((svmr_cm4$t[1]+svmr_cm4$t[4])/sum(svmr_cm4$t))*100

svmke_train4 <- ksvm(as.matrix(train), as.factor(train4$diagnosis), kernel='vanilladot')
svmke_pred4 <- predict(svmke_train4, validation, type='response')
svmke_cm4 =CrossTable(x = validation4$diagnosis , y = svmke_pred4, prop.chisq = FALSE)
svmke_accuracy4 <- ((svmke_cm4$t[1]+svmke_cm4$t[4])/sum(svmke_cm4$t))*100

svmkp_train4 <- ksvm(as.matrix(train), as.factor(train4$diagnosis), kernel='polydot')
svmkp_pred4 <- predict(svmkp_train4, validation, type='response')
svmkp_cm4 =CrossTable(x = validation4$diagnosis , y = svmkp_pred4, prop.chisq = FALSE)
svmkp_accuracy4 <- ((svmkp_cm4$t[1]+svmkp_cm4$t[4])/sum(svmkp_cm4$t))*100

svmkr_train4 <- ksvm(as.matrix(train), as.factor(train4$diagnosis), kernel='rbfdot')
svmkr_pred4 <- predict(svmkr_train4, validation, type='response')
svmkr_cm4 =CrossTable(x = validation4$diagnosis , y = svmkr_pred4, prop.chisq = FALSE)
svmkr_accuracy4 <- ((svmkr_cm4$t[1]+svmkr_cm4$t[4])/sum(svmkr_cm4$t))*100

# 5 of 5 ----
train5 <- normalized_bcd[Partitions$Resample5,]
validation5 <- normalized_bcd[-Partitions$Resample5,]
train <- subset(train5, select=-diagnosis)
validation <- subset(validation5, select=-diagnosis)

svme_train5 <-  svm(as.factor(diagnosis)~., train5, probability = TRUE, type = "C-classification", kernel = "linear")
svme_pred5 <- predict(svme_train5, validation, probability = FALSE, decision.values = TRUE)
svme_cm5 =CrossTable(x = validation5$diagnosis , y = svme_pred5, prop.chisq = FALSE)
svme_accuracy5 <- ((svme_cm5$t[1]+svme_cm5$t[4])/sum(svme_cm5$t))*100

svmp_train5 <- svm(as.factor(diagnosis)~., train5, probability = TRUE, type = "C-classification", kernel = "polynomial", degree = 3)
svmp_pred5 <- predict(svmp_train5, validation, probability = FALSE, decision.values = TRUE)
svmp_cm5 =CrossTable(x = validation5$diagnosis , y = svmp_pred5, prop.chisq = FALSE)
svmp_accuracy5 <- ((svmp_cm5$t[1]+svmp_cm5$t[4])/sum(svmp_cm5$t))*100

svmr_train5 <- svm(as.factor(diagnosis)~., train5, probability = TRUE, type = "C-classification", kernel = "radial", gamma = 0.1)
svmr_pred5 <- predict(svmr_train5, validation, probability = FALSE, decision.values = TRUE)
svmr_cm5 =CrossTable(x = validation5$diagnosis , y = svmr_pred5, prop.chisq = FALSE)
svmr_accuracy5 <- ((svmr_cm5$t[1]+svmr_cm5$t[4])/sum(svmr_cm5$t))*100

svmke_train5 <- ksvm(as.matrix(train), as.factor(train5$diagnosis), kernel='vanilladot')
svmke_pred5 <- predict(svmke_train5, validation, type='response')
svmke_cm5 =CrossTable(x = validation5$diagnosis , y = svmke_pred5, prop.chisq = FALSE)
svmke_accuracy5 <- ((svmke_cm5$t[1]+svmke_cm5$t[4])/sum(svmke_cm5$t))*100

svmkp_train5 <- ksvm(as.matrix(train), as.factor(train5$diagnosis), kernel='polydot')
svmkp_pred5 <- predict(svmkp_train5, validation, type='response')
svmkp_cm5 =CrossTable(x = validation5$diagnosis , y = svmkp_pred5, prop.chisq = FALSE)
svmkp_accuracy5 <- ((svmkp_cm5$t[1]+svmkp_cm5$t[4])/sum(svmkp_cm5$t))*100

svmkr_train5 <- ksvm(as.matrix(train), as.factor(train5$diagnosis), kernel='rbfdot')
svmkr_pred5 <- predict(svmkr_train5, validation, type='response')
svmkr_cm5 =CrossTable(x = validation5$diagnosis , y = svmkr_pred5, prop.chisq = FALSE)
svmkr_accuracy5 <- ((svmkr_cm5$t[1]+svmkr_cm5$t[4])/sum(svmkr_cm5$t))*100

# Accuracies ----

svme_accuracies <- c(svme_accuracy1,svme_accuracy2,svme_accuracy3,svme_accuracy4,svme_accuracy5)
svmp_accuracies <- c(svmp_accuracy1,svmp_accuracy2,svmp_accuracy3,svmp_accuracy4,svmp_accuracy5)
svmr_accuracies <- c(svmr_accuracy1,svmr_accuracy2,svmr_accuracy3,svmr_accuracy4,svmr_accuracy5)
svmke_accuracies <- c(svmke_accuracy1,svmke_accuracy2,svmke_accuracy3,svmke_accuracy4,svmke_accuracy5)
svmkp_accuracies <- c(svmkp_accuracy1,svmkp_accuracy2,svmkp_accuracy3,svmkp_accuracy4,svmkp_accuracy5)
svmkr_accuracies <- c(svmkr_accuracy1,svmkr_accuracy2,svmkr_accuracy3,svmkr_accuracy4,svmkr_accuracy5)

svme_accuracies
svmp_accuracies
svmr_accuracies
svmke_accuracies
svmkp_accuracies
svmkr_accuracies

# Validation scatterplots ----
scatter_validation1 <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 1
scatter_validation2 <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 2
scatter_validation3 <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 3
scatter_validation4 <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 4
scatter_validation5 <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 5

# Prediction scatterplots svme
svme_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(svme_pred1))) + geom_point()
svme_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(svme_pred2))) + geom_point()
svme_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(svme_pred3))) + geom_point()
svme_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(svme_pred4))) + geom_point()
svme_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(svme_pred5))) + geom_point()

# Prediction scatterplots svmp
svmp_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(svmp_pred1))) + geom_point()
svmp_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(svmp_pred2))) + geom_point()
svmp_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(svmp_pred3))) + geom_point()
svmp_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(svmp_pred4))) + geom_point()
svmp_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(svmp_pred5))) + geom_point()

# Prediction scatterplots svmr
svmr_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(svmr_pred1))) + geom_point()
svmr_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(svmr_pred2))) + geom_point()
svmr_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(svmr_pred3))) + geom_point()
svmr_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(svmr_pred4))) + geom_point()
svmr_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(svmr_pred5))) + geom_point()

# Prediction scatterplots svmke
svmke_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(svmke_pred1))) + geom_point()
svmke_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(svmke_pred2))) + geom_point()
svmke_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(svmke_pred3))) + geom_point()
svmke_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(svmke_pred4))) + geom_point()
svmke_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(svmke_pred5))) + geom_point()

# Prediction scatterplots svmkp
svmkp_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(svmkp_pred1))) + geom_point()
svmkp_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(svmkp_pred2))) + geom_point()
svmkp_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(svmkp_pred3))) + geom_point()
svmkp_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(svmkp_pred4))) + geom_point()
svmkp_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(svmkp_pred5))) + geom_point()

# Prediction scatterplots svmkr
svmkr_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(svmkr_pred1))) + geom_point()
svmkr_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(svmkr_pred2))) + geom_point()
svmkr_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(svmkr_pred3))) + geom_point()
svmkr_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(svmkr_pred4))) + geom_point()
svmkr_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(svmkr_pred5))) + geom_point()


# export plots ----

pdf(file = 'ElenesAGN_A2_Problem2_GraphicOutputs.pdf')

multiplot(scatter_validation1,svme_scatter_1_pred)
multiplot(scatter_validation2,svme_scatter_2_pred)
multiplot(scatter_validation3,svme_scatter_3_pred)
multiplot(scatter_validation4,svme_scatter_4_pred)
multiplot(scatter_validation5,svme_scatter_5_pred)

multiplot(scatter_validation1,svmp_scatter_1_pred)
multiplot(scatter_validation2,svmp_scatter_2_pred)
multiplot(scatter_validation3,svmp_scatter_3_pred)
multiplot(scatter_validation4,svmp_scatter_4_pred)
multiplot(scatter_validation5,svmp_scatter_5_pred)

multiplot(scatter_validation1,svmr_scatter_1_pred)
multiplot(scatter_validation2,svmr_scatter_2_pred)
multiplot(scatter_validation3,svmr_scatter_3_pred)
multiplot(scatter_validation4,svmr_scatter_4_pred)
multiplot(scatter_validation5,svmr_scatter_5_pred)

multiplot(scatter_validation1,svmke_scatter_1_pred)
multiplot(scatter_validation2,svmke_scatter_2_pred)
multiplot(scatter_validation3,svmke_scatter_3_pred)
multiplot(scatter_validation4,svmke_scatter_4_pred)
multiplot(scatter_validation5,svmke_scatter_5_pred)

multiplot(scatter_validation1,svmkp_scatter_1_pred)
multiplot(scatter_validation2,svmkp_scatter_2_pred)
multiplot(scatter_validation3,svmkp_scatter_3_pred)
multiplot(scatter_validation4,svmkp_scatter_4_pred)
multiplot(scatter_validation5,svmkp_scatter_5_pred)

multiplot(scatter_validation1,svmkr_scatter_1_pred)
multiplot(scatter_validation2,svmkr_scatter_2_pred)
multiplot(scatter_validation3,svmkr_scatter_3_pred)
multiplot(scatter_validation4,svmkr_scatter_4_pred)
multiplot(scatter_validation5,svmkr_scatter_5_pred)

dev.off()

# png(filename="Pairs_correlograms.png",width=3840,height=2160)
# multiplot(corrs)
# dev.off()

# export confusion matrices ----

svme_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(svme_pred1))
svme_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(svme_pred2))
svme_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(svme_pred3))
svme_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(svme_pred4))
svme_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(svme_pred5))

svmp_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(svmp_pred1))
svmp_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(svmp_pred2))
svmp_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(svmp_pred3))
svmp_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(svmp_pred4))
svmp_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(svmp_pred5))

svmr_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(svmr_pred1))
svmr_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(svmr_pred2))
svmr_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(svmr_pred3))
svmr_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(svmr_pred4))
svmr_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(svmr_pred5))

svmke_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(svmke_pred1))
svmke_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(svmke_pred2))
svmke_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(svmke_pred3))
svmke_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(svmke_pred4))
svmke_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(svmke_pred5))

svmkp_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(svmkp_pred1))
svmkp_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(svmkp_pred2))
svmkp_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(svmkp_pred3))
svmkp_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(svmkp_pred4))
svmkp_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(svmkp_pred5))

svmkr_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(svmkr_pred1))
svmkr_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(svmkr_pred2))
svmkr_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(svmkr_pred3))
svmkr_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(svmkr_pred4))
svmkr_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(svmkr_pred5))

sink(file = 'ElenesAGN_A2_Problem2_TextOutputs.txt')

print(svme_confusionmatrix1)
print(svme_confusionmatrix2)
print(svme_confusionmatrix3)
print(svme_confusionmatrix4)
print(svme_confusionmatrix5)

print(svmp_confusionmatrix1)
print(svmp_confusionmatrix2)
print(svmp_confusionmatrix3)
print(svmp_confusionmatrix4)
print(svmp_confusionmatrix5)

print(svmr_confusionmatrix1)
print(svmr_confusionmatrix2)
print(svmr_confusionmatrix3)
print(svmr_confusionmatrix4)
print(svmr_confusionmatrix5)

print(svmke_confusionmatrix1)
print(svmke_confusionmatrix2)
print(svmke_confusionmatrix3)
print(svmke_confusionmatrix4)
print(svmke_confusionmatrix5)

print(svmkp_confusionmatrix1)
print(svmkp_confusionmatrix2)
print(svmkp_confusionmatrix3)
print(svmkp_confusionmatrix4)
print(svmkp_confusionmatrix5)

print(svmkr_confusionmatrix1)
print(svmkr_confusionmatrix2)
print(svmkr_confusionmatrix3)
print(svmkr_confusionmatrix4)
print(svmkr_confusionmatrix5)

sink()

