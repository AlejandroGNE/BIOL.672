# Alejandro G.N. Elenes
# Code written in Windows 10 Home
# R packages needed to run this script: ----
library(readr)
library(class)
library(gmodels)
library(ggplot2)
library(GGally)
library(corrplot)
library(Rmisc)
library(e1071)
library(caret)
library(MASS)
library(mda)
library(dplyr)
library(kernlab)
library(randomForest)
library(neuralnet)
library(ada)
library(scales)
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

knn_pred1 <- knn(train,validation,train1$diagnosis, k=sqrt(nrow(bcd)), l=0, prob=FALSE, use.all=TRUE)
summary(knn_pred1)
knn_cm1 =CrossTable(x = validation1$diagnosis , y = knn_pred1, prop.chisq = FALSE)
knn_accuracy1 <- ((knn_cm1$t[1]+knn_cm1$t[4])/sum(knn_cm1$t))*100
knn_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(knn_pred1))

NBtrain1 <-  naiveBayes(as.factor(diagnosis)~., train1, laplace = 0)
NBpred1 <- predict(NBtrain1, validation, probability = FALSE, decision.values = TRUE)
NB_cm1 =CrossTable(x = validation1$diagnosis , y = NBpred1, prop.chisq = FALSE)
NB_accuracy1 <- ((NB_cm1$t[1]+NB_cm1$t[4])/sum(NB_cm1$t))*100

ldatrain1 <- lda(as.factor(diagnosis)~., train1)
ldapred1 <- predict(ldatrain1, validation, prior=ldatrain1$prior, method=c("plug-in", "predictive", "debiased"))
lda_cm1 =CrossTable(x = validation1$diagnosis , y = ldapred1$class, prop.chisq = FALSE)
lda_accuracy1 <- ((lda_cm1$t[1]+lda_cm1$t[4])/sum(lda_cm1$t))*100

qdatrain1 <- qda(as.factor(diagnosis)~., train1)
qdapred1 <- predict(qdatrain1, validation, prior=qdatrain1$prior)
qda_cm1 =CrossTable(x = validation1$diagnosis , y = qdapred1$class, prop.chisq = FALSE)
qda_accuracy1 <- ((qda_cm1$t[1]+qda_cm1$t[4])/sum(qda_cm1$t))*100

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

rndfor_train1 <- randomForest(as.factor(diagnosis)~., train1, ntree=500)
rndfor_pred1 <- predict(rndfor_train1, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm1 <- CrossTable(x = validation1$diagnosis , y = rndfor_pred1, prop.chisq = FALSE)
rndfor_accuracy1 <- ((rndfor_cm1$t[1]+rndfor_cm1$t[4])/sum(rndfor_cm1$t))*100

ada_train1 <- ada(as.factor(diagnosis)~., train1, 500)
ada_pred1 <- predict(ada_train1, validation, probability = FALSE, decision.values = TRUE)
ada_cm1 <- CrossTable(x = validation1$diagnosis , y = ada_pred1, prop.chisq = FALSE)
ada_accuracy1 <- ((ada_cm1$t[1]+ada_cm1$t[4])/sum(ada_cm1$t))*100

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

knn_pred2 <- knn(train,validation,train2$diagnosis, k=sqrt(nrow(bcd)), l=0, prob=FALSE, use.all=TRUE)
summary(knn_pred2)
knn_cm2 =CrossTable(x = validation2$diagnosis , y = knn_pred2, prop.chisq = FALSE)
knn_accuracy2 <- ((knn_cm2$t[1]+knn_cm2$t[4])/sum(knn_cm2$t))*100

NBtrain2 <-  naiveBayes(as.factor(diagnosis)~., train2, laplace = 0)
NBpred2 <- predict(NBtrain2, validation, probability = FALSE, decision.values = TRUE)
NB_cm2 =CrossTable(x = validation2$diagnosis , y = NBpred2, prop.chisq = FALSE)
NB_accuracy2 <- ((NB_cm2$t[1]+NB_cm2$t[4])/sum(NB_cm2$t))*100

ldatrain2 <- lda(as.factor(diagnosis)~., train2)
ldapred2 <- predict(ldatrain2, validation, prior=ldatrain2$prior, method=c("plug-in", "predictive", "debiased"))
lda_cm2 =CrossTable(x = validation2$diagnosis , y = ldapred2$class, prop.chisq = FALSE)
lda_accuracy2 <- ((lda_cm2$t[1]+lda_cm2$t[4])/sum(lda_cm2$t))*100

qdatrain2 <- qda(as.factor(diagnosis)~., train2)
qdapred2 <- predict(qdatrain2, validation, prior=qdatrain2$prior)
qda_cm2 =CrossTable(x = validation2$diagnosis , y = qdapred2$class, prop.chisq = FALSE)
qda_accuracy2 <- ((qda_cm2$t[1]+qda_cm2$t[4])/sum(qda_cm2$t))*100

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

rndfor_train2 <- randomForest(as.factor(diagnosis)~., train2, ntree=500)
rndfor_pred2 <- predict(rndfor_train2, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm2 <- CrossTable(x = validation2$diagnosis , y = rndfor_pred2, prop.chisq = FALSE)
rndfor_accuracy2 <- ((rndfor_cm2$t[1]+rndfor_cm2$t[4])/sum(rndfor_cm2$t))*100

ada_train2 <- ada(as.factor(diagnosis)~., train2, 500)
ada_pred2 <- predict(ada_train2, validation, probability = FALSE, decision.values = TRUE)
ada_cm2 <- CrossTable(x = validation2$diagnosis , y = ada_pred2, prop.chisq = FALSE)
ada_accuracy2 <- ((ada_cm2$t[1]+ada_cm2$t[4])/sum(ada_cm2$t))*100

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

knn_pred3 <- knn(train,validation,train3$diagnosis, k=sqrt(nrow(bcd)), l=0, prob=FALSE, use.all=TRUE)
summary(knn_pred3)
knn_cm3 =CrossTable(x = validation3$diagnosis , y = knn_pred3, prop.chisq = FALSE)
knn_accuracy3 <- ((knn_cm3$t[1]+knn_cm3$t[4])/sum(knn_cm3$t))*100

NBtrain3 <-  naiveBayes(as.factor(diagnosis)~., train3, laplace = 0)
NBpred3 <- predict(NBtrain3, validation, probability = FALSE, decision.values = TRUE)
NB_cm3 =CrossTable(x = validation3$diagnosis , y = NBpred3, prop.chisq = FALSE)
NB_accuracy3 <- ((NB_cm3$t[1]+NB_cm3$t[4])/sum(NB_cm3$t))*100

ldatrain3 <- lda(as.factor(diagnosis)~., train3)
ldapred3 <- predict(ldatrain3, validation, prior=ldatrain3$prior, method=c("plug-in", "predictive", "debiased"))
lda_cm3 =CrossTable(x = validation3$diagnosis , y = ldapred3$class, prop.chisq = FALSE)
lda_accuracy3 <- ((lda_cm3$t[1]+lda_cm3$t[4])/sum(lda_cm3$t))*100

qdatrain3 <- qda(as.factor(diagnosis)~., train3)
qdapred3 <- predict(qdatrain3, validation, prior=qdatrain3$prior)
qda_cm3 =CrossTable(x = validation3$diagnosis , y = qdapred3$class, prop.chisq = FALSE)
qda_accuracy3 <- ((qda_cm3$t[1]+qda_cm3$t[4])/sum(qda_cm3$t))*100

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

rndfor_train3 <- randomForest(as.factor(diagnosis)~., train3, ntree=500)
rndfor_pred3 <- predict(rndfor_train3, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm3 <- CrossTable(x = validation3$diagnosis , y = rndfor_pred3, prop.chisq = FALSE)
rndfor_accuracy3 <- ((rndfor_cm3$t[1]+rndfor_cm3$t[4])/sum(rndfor_cm3$t))*100

ada_train3 <- ada(as.factor(diagnosis)~., train3, 500)
ada_pred3 <- predict(ada_train3, validation, probability = FALSE, decision.values = TRUE)
ada_cm3 <- CrossTable(x = validation3$diagnosis , y = ada_pred3, prop.chisq = FALSE)
ada_accuracy3 <- ((ada_cm3$t[1]+ada_cm3$t[4])/sum(ada_cm3$t))*100

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

knn_pred4 <- knn(train,validation,train4$diagnosis, k=sqrt(nrow(bcd)), l=0, prob=FALSE, use.all=TRUE)
summary(knn_pred4)
knn_cm4 =CrossTable(x = validation4$diagnosis , y = knn_pred4, prop.chisq = FALSE)
knn_accuracy4 <- ((knn_cm4$t[1]+knn_cm4$t[4])/sum(knn_cm4$t))*100

NBtrain4 <-  naiveBayes(as.factor(diagnosis)~., train4, laplace = 0)
NBpred4 <- predict(NBtrain4, validation, probability = FALSE, decision.values = TRUE)
NB_cm4 =CrossTable(x = validation4$diagnosis , y = NBpred4, prop.chisq = FALSE)
NB_accuracy4 <- ((NB_cm4$t[1]+NB_cm4$t[4])/sum(NB_cm4$t))*100

ldatrain4 <- lda(as.factor(diagnosis)~., train4)
ldapred4 <- predict(ldatrain4, validation, prior=ldatrain4$prior, method=c("plug-in", "predictive", "debiased"))
lda_cm4 =CrossTable(x = validation4$diagnosis , y = ldapred4$class, prop.chisq = FALSE)
lda_accuracy4 <- ((lda_cm4$t[1]+lda_cm4$t[4])/sum(lda_cm4$t))*100

qdatrain4 <- qda(as.factor(diagnosis)~., train4)
qdapred4 <- predict(qdatrain4, validation, prior=qdatrain4$prior)
qda_cm4 =CrossTable(x = validation4$diagnosis , y = qdapred4$class, prop.chisq = FALSE)
qda_accuracy4 <- ((qda_cm4$t[1]+qda_cm4$t[4])/sum(qda_cm4$t))*100

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

rndfor_train4 <- randomForest(as.factor(diagnosis)~., train4, ntree=500)
rndfor_pred4 <- predict(rndfor_train4, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm4 <- CrossTable(x = validation4$diagnosis , y = rndfor_pred4, prop.chisq = FALSE)
rndfor_accuracy4 <- ((rndfor_cm4$t[1]+rndfor_cm4$t[4])/sum(rndfor_cm4$t))*100

ada_train4 <- ada(as.factor(diagnosis)~., train4, 500)
ada_pred4 <- predict(ada_train4, validation, probability = FALSE, decision.values = TRUE)
ada_cm4 <- CrossTable(x = validation4$diagnosis , y = ada_pred4, prop.chisq = FALSE)
ada_accuracy4 <- ((ada_cm4$t[1]+ada_cm4$t[4])/sum(ada_cm4$t))*100

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

knn_pred5 <- knn(train,validation,train5$diagnosis, k=sqrt(nrow(bcd)), l=0, prob=FALSE, use.all=TRUE)
summary(knn_pred5)
knn_cm5 =CrossTable(x = validation5$diagnosis , y = knn_pred5, prop.chisq = FALSE)
knn_accuracy5 <- ((knn_cm5$t[1]+knn_cm5$t[4])/sum(knn_cm5$t))*100

NBtrain5 <-  naiveBayes(as.factor(diagnosis)~., train5, laplace = 0)
NBpred5 <- predict(NBtrain5, validation, probability = FALSE, decision.values = TRUE)
NB_cm5 =CrossTable(x = validation5$diagnosis , y = NBpred5, prop.chisq = FALSE)
NB_accuracy5 <- ((NB_cm5$t[1]+NB_cm5$t[4])/sum(NB_cm5$t))*100

ldatrain5 <- lda(as.factor(diagnosis)~., train5)
ldapred5 <- predict(ldatrain5, validation, prior=ldatrain5$prior, method=c("plug-in", "predictive", "debiased"))
lda_cm5 =CrossTable(x = validation5$diagnosis , y = ldapred5$class, prop.chisq = FALSE)
lda_accuracy5 <- ((lda_cm5$t[1]+lda_cm5$t[4])/sum(lda_cm5$t))*100

qdatrain5 <- qda(as.factor(diagnosis)~., train5)
qdapred5 <- predict(qdatrain5, validation, prior=qdatrain5$prior)
qda_cm5 =CrossTable(x = validation5$diagnosis , y = qdapred5$class, prop.chisq = FALSE)
qda_accuracy5 <- ((qda_cm5$t[1]+qda_cm5$t[4])/sum(qda_cm5$t))*100

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

rndfor_train5 <- randomForest(as.factor(diagnosis)~., train5, ntree=500)
rndfor_pred5 <- predict(rndfor_train5, validation, probability = FALSE, decision.values = TRUE)
rndfor_cm5 <- CrossTable(x = validation5$diagnosis , y = rndfor_pred5, prop.chisq = FALSE)
rndfor_accuracy5 <- ((rndfor_cm5$t[1]+rndfor_cm5$t[4])/sum(rndfor_cm5$t))*100

ada_train5 <- ada(as.factor(diagnosis)~., train5, 500)
ada_pred5 <- predict(ada_train5, validation, probability = FALSE, decision.values = TRUE)
ada_cm5 <- CrossTable(x = validation5$diagnosis , y = ada_pred5, prop.chisq = FALSE)
ada_accuracy5 <- ((ada_cm5$t[1]+ada_cm5$t[4])/sum(ada_cm5$t))*100

nn_train5 <- neuralnet(as.factor(diagnosis)~., train5, hidden=3, linear.output = FALSE)
nn_pred5 <- predict(nn_train5, validation, rep=1, all.units=FALSE)
nn_cm5 <- CrossTable(x = validation5$diagnosis , y = max.col(nn_pred5), prop.chisq = FALSE)
nn_accuracy5 <- ((nn_cm5$t[1]+nn_cm5$t[4])/sum(nn_cm5$t))*100
print(plot(nn_train5))

# Accuracies ----

knn_accuracies <- c(knn_accuracy1,knn_accuracy2,knn_accuracy3,knn_accuracy4,knn_accuracy5)
NB_accuracies <- c(NB_accuracy1,NB_accuracy2,NB_accuracy3,NB_accuracy4,NB_accuracy5)
lda_accuracies <- c(lda_accuracy1,lda_accuracy2,lda_accuracy3,lda_accuracy4,lda_accuracy5)
qda_accuracies <- c(qda_accuracy1,qda_accuracy2,qda_accuracy3,qda_accuracy4,qda_accuracy5)
svme_accuracies <- c(svme_accuracy1,svme_accuracy2,svme_accuracy3,svme_accuracy4,svme_accuracy5)
svmp_accuracies <- c(svmp_accuracy1,svmp_accuracy2,svmp_accuracy3,svmp_accuracy4,svmp_accuracy5)
svmr_accuracies <- c(svmr_accuracy1,svmr_accuracy2,svmr_accuracy3,svmr_accuracy4,svmr_accuracy5)
svmke_accuracies <- c(svmke_accuracy1,svmke_accuracy2,svmke_accuracy3,svmke_accuracy4,svmke_accuracy5)
svmkp_accuracies <- c(svmkp_accuracy1,svmkp_accuracy2,svmkp_accuracy3,svmkp_accuracy4,svmkp_accuracy5)
svmkr_accuracies <- c(svmkr_accuracy1,svmkr_accuracy2,svmkr_accuracy3,svmkr_accuracy4,svmkr_accuracy5)
rndfor_accuracies <- c(rndfor_accuracy1,rndfor_accuracy2,rndfor_accuracy3,rndfor_accuracy4,rndfor_accuracy5)
nn_accuracies <- c(nn_accuracy1,nn_accuracy2,nn_accuracy3,nn_accuracy4,nn_accuracy5)
ada_accuracies <- c(ada_accuracy1,ada_accuracy2,ada_accuracy3,ada_accuracy4,ada_accuracy5)

knn_accuracies
NB_accuracies
lda_accuracies
qda_accuracies
svme_accuracies
svmp_accuracies
svmr_accuracies
svmke_accuracies
svmkp_accuracies
svmkr_accuracies
rndfor_accuracies
ada_accuracies
nn_accuracies

# Validation scatterplots ----
scatter_validation1 <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 1
scatter_validation2 <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 2
scatter_validation3 <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 3
scatter_validation4 <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 4
scatter_validation5 <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + geom_point() # group 5

# Prediction scatterplots knn
knn_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(knn_pred1))) + geom_point() # group 1
knn_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(knn_pred2))) + geom_point() # group 2
knn_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(knn_pred3))) + geom_point() # group 3
knn_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(knn_pred4))) + geom_point() # group 4
knn_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(knn_pred5))) + geom_point() # group 5

# Prediction scatterplots naive bayes
NB_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(NBpred1))) + geom_point()
NB_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(NBpred2))) + geom_point()
NB_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(NBpred3))) + geom_point()
NB_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(NBpred4))) + geom_point()
NB_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(NBpred5))) + geom_point()

# Prediction scatterplots lda
lda_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(ldapred1$class))) + geom_point()
lda_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(ldapred2$class))) + geom_point()
lda_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(ldapred3$class))) + geom_point()
lda_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(ldapred4$class))) + geom_point()
lda_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(ldapred5$class))) + geom_point()

# Prediction scatterplots qda
qda_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(qdapred1$class))) + geom_point()
qda_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(qdapred2$class))) + geom_point()
qda_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(qdapred3$class))) + geom_point()
qda_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(qdapred4$class))) + geom_point()
qda_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(qdapred5$class))) + geom_point()

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

# Prediction scatterplots nn
nn_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred1)))) + geom_point()
nn_scatter_2_pred <- ggplot(validation2, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred2)))) + geom_point()
nn_scatter_3_pred <- ggplot(validation3, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred3)))) + geom_point()
nn_scatter_4_pred <- ggplot(validation4, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred4)))) + geom_point()
nn_scatter_5_pred <- ggplot(validation5, aes(area_worst, concave.points_worst, colour = factor(max.col(nn_pred5)))) + geom_point()

# export plots ----

pdf(file = 'ElenesAGN_TakeHomeFinalExam_GraphicOutputs.pdf')

multiplot(scatter_validation1,knn_scatter_1_pred)
multiplot(scatter_validation2,knn_scatter_2_pred)
multiplot(scatter_validation3,knn_scatter_3_pred)
multiplot(scatter_validation4,knn_scatter_4_pred)
multiplot(scatter_validation5,knn_scatter_5_pred)

multiplot(scatter_validation1,NB_scatter_1_pred)
multiplot(scatter_validation2,NB_scatter_2_pred)
multiplot(scatter_validation3,NB_scatter_3_pred)
multiplot(scatter_validation4,NB_scatter_4_pred)
multiplot(scatter_validation5,NB_scatter_5_pred)

multiplot(scatter_validation1,lda_scatter_1_pred)
multiplot(scatter_validation2,lda_scatter_2_pred)
multiplot(scatter_validation3,lda_scatter_3_pred)
multiplot(scatter_validation4,lda_scatter_4_pred)
multiplot(scatter_validation5,lda_scatter_5_pred)

multiplot(scatter_validation1,qda_scatter_1_pred)
multiplot(scatter_validation2,qda_scatter_2_pred)
multiplot(scatter_validation3,qda_scatter_3_pred)
multiplot(scatter_validation4,qda_scatter_4_pred)
multiplot(scatter_validation5,qda_scatter_5_pred)

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

multiplot(scatter_validation1,nn_scatter_1_pred)
multiplot(scatter_validation2,nn_scatter_2_pred)
multiplot(scatter_validation3,nn_scatter_3_pred)
multiplot(scatter_validation4,nn_scatter_4_pred)
multiplot(scatter_validation5,nn_scatter_5_pred)

dev.off()

png(filename="Pairs_correlograms.png",width=3840,height=2160)
multiplot(corrs)
dev.off()

# confusion matrices ----
knn_confusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(knn_pred1))
knn_confusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(knn_pred2))
knn_confusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(knn_pred3))
knn_confusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(knn_pred4))
knn_confusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(knn_pred5))

NBconfusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(NBpred1))
NBconfusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(NBpred2))
NBconfusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(NBpred3))
NBconfusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(NBpred4))
NBconfusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(NBpred5))

ldaconfusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(ldapred1$class))
ldaconfusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(ldapred2$class))
ldaconfusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(ldapred3$class))
ldaconfusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(ldapred4$class))
ldaconfusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(ldapred5$class))

qdaconfusionmatrix1 <- confusionMatrix(as.factor(validation1$diagnosis),as.factor(qdapred1$class))
qdaconfusionmatrix2 <- confusionMatrix(as.factor(validation2$diagnosis),as.factor(qdapred2$class))
qdaconfusionmatrix3 <- confusionMatrix(as.factor(validation3$diagnosis),as.factor(qdapred3$class))
qdaconfusionmatrix4 <- confusionMatrix(as.factor(validation4$diagnosis),as.factor(qdapred4$class))
qdaconfusionmatrix5 <- confusionMatrix(as.factor(validation5$diagnosis),as.factor(qdapred5$class))

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

nn_confusionmatrix1 <- confusionMatrix(as.factor(as.integer(as.factor(validation1$diagnosis))),as.factor(max.col(nn_pred1)))
nn_confusionmatrix2 <- confusionMatrix(as.factor(as.integer(as.factor(validation2$diagnosis))),as.factor(max.col(nn_pred2)))
nn_confusionmatrix3 <- confusionMatrix(as.factor(as.integer(as.factor(validation3$diagnosis))),as.factor(max.col(nn_pred3)))
nn_confusionmatrix4 <- confusionMatrix(as.factor(as.integer(as.factor(validation4$diagnosis))),as.factor(max.col(nn_pred4)))
nn_confusionmatrix5 <- confusionMatrix(as.factor(as.integer(as.factor(validation5$diagnosis))),as.factor(max.col(nn_pred5)))

# export text outputs ----

sink(file = 'ElenesAGN_TakeHomeFinalExam_TextOutputs.txt')

writeLines(" \nAbout the data")

writeLines(" \n\n The explanatory variables are the mean, standard deviation and the worst measurement of the following features: \n")
writeLines(" 1. Radius")
writeLines(" 2. Texture")
writeLines(" 3. Perimeter")
writeLines(" 4. Area")
writeLines(" 5. Smoothness")
writeLines(" 6. Compactness")
writeLines(" 7. Concavity")
writeLines(" 8. Concave points")
writeLines(" 9. Symmetry")
writeLines(" 10. Fractal dimension")
writeLines(" \n\n I can see some trends by visual inspection of the correlogram made with ggpairs. \n")
writeLines(" \n\n In the diagonal, I show the distribution of the data for malignant (blue) and benign (red) tumors; the more separated these red and blue distributions are, the better the data will help in predicting the malignant or benign property of the tumor. If these distributions overlap too much, data will not be as useful for classification. \n")
writeLines(" \n\n The overlap seems to be larger when comparing the standard deviation of the measurements, which means the variability of the data is similar in both cases. If we compare the overlap in mean values vs that of worst measurements, they seem to mirror each other, however, the overlap slightly decreases in worst measurements, indicating that these could be more useful than mean and standard deviation values for classification. \n")
writeLines(" \n\n Below the diagonal we see the specific correlations. Some strong correlations show that features are connected, for example, radius, perimeter and area are all measures of size, and thus they are strongly correlated. Texture and symmetry do not correlate with any other feature. Another group of features with a smaller and more dispersive correlation is that of smoothness, compactness, concavity and concave points. Fractal dimension does not correlate with any other feature except for smoothness, compactness and concavity, but it's a small correlation and can only be observed in the set of worst measurements. \n")
writeLines(" \n\n Given these observations, for the scatterplots showing output results I decided to use the worst measurements for one feature from each of the dominant groups: worst measurement of area, from the group of features related to size, and worst measurement of concave points, from the group of features related to shape. \n")

writeLines(" \n\n Confusion matrices for 5-fold cross validation of k nearest neighbors \n")
print(knn_confusionmatrix1)
print(knn_confusionmatrix2)
print(knn_confusionmatrix3)
print(knn_confusionmatrix4)
print(knn_confusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of Naive Bayes classifier \n")
print(NBconfusionmatrix1)
print(NBconfusionmatrix2)
print(NBconfusionmatrix3)
print(NBconfusionmatrix4)
print(NBconfusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of linear discriminant analysis \n")
print(ldaconfusionmatrix1)
print(ldaconfusionmatrix2)
print(ldaconfusionmatrix3)
print(ldaconfusionmatrix4)
print(ldaconfusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of quadratic discriminant analysis \n")
print(qdaconfusionmatrix1)
print(qdaconfusionmatrix2)
print(qdaconfusionmatrix3)
print(qdaconfusionmatrix4)
print(qdaconfusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of linear SVM \n")
print(svme_confusionmatrix1)
print(svme_confusionmatrix2)
print(svme_confusionmatrix3)
print(svme_confusionmatrix4)
print(svme_confusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of polynomial SVM \n")
print(svmp_confusionmatrix1)
print(svmp_confusionmatrix2)
print(svmp_confusionmatrix3)
print(svmp_confusionmatrix4)
print(svmp_confusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of radial SVM \n")
print(svmr_confusionmatrix1)
print(svmr_confusionmatrix2)
print(svmr_confusionmatrix3)
print(svmr_confusionmatrix4)
print(svmr_confusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of tuned linear SVM \n")
print(svmke_confusionmatrix1)
print(svmke_confusionmatrix2)
print(svmke_confusionmatrix3)
print(svmke_confusionmatrix4)
print(svmke_confusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of tuned qudratic SVM \n")
print(svmkp_confusionmatrix1)
print(svmkp_confusionmatrix2)
print(svmkp_confusionmatrix3)
print(svmkp_confusionmatrix4)
print(svmkp_confusionmatrix5)

writeLines(" \n\n Confusion matrices for 5-fold cross validation of tuned radial SVM \n")
print(svmkr_confusionmatrix1)
print(svmkr_confusionmatrix2)
print(svmkr_confusionmatrix3)
print(svmkr_confusionmatrix4)
print(svmkr_confusionmatrix5)

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

writeLines(" \n\n Confusion matrices for 5-fold cross validation of neural net model \n")
print(nn_confusionmatrix1)
print(nn_confusionmatrix2)
print(nn_confusionmatrix3)
print(nn_confusionmatrix4)
print(nn_confusionmatrix5)

writeLines(" \n\n Accuracies \n")
print(knn_accuracies)
print(NB_accuracies)
print(lda_accuracies)
print(qda_accuracies)
print(svme_accuracies)
print(svmp_accuracies)
print(svmr_accuracies)
print(svmke_accuracies)
print(svmkp_accuracies)
print(svmkr_accuracies)
print(rndfor_accuracies)
print(ada_accuracies)
print(nn_accuracies)

writeLines(" \n\n Mean accuracies \n")
print(mean(knn_accuracies))
print(mean(NB_accuracies))
print(mean(lda_accuracies))
print(mean(qda_accuracies))
print(mean(svme_accuracies))
print(mean(svmp_accuracies))
print(mean(svmr_accuracies))
print(mean(svmke_accuracies))
print(mean(svmkp_accuracies))
print(mean(svmkr_accuracies))
print(mean(rndfor_accuracies))
print(mean(ada_accuracies))
print(mean(nn_accuracies))


sink()

# Comparing methods overall ----

titles <- cbind('knn','Naive Bayes','linear discriminant analysis','quadratic discriminant analysis','linear SVM','quadratic SVM','radial SVM','tuned linear SVM','tuned quadratic SVM','tuned radial SVM','random forest','boosted gradient','neural network')

accuracies <- cbind(knn_accuracies,NB_accuracies,lda_accuracies,qda_accuracies,svme_accuracies,svmp_accuracies,svmr_accuracies,svmke_accuracies,svmkp_accuracies,svmkr_accuracies,rndfor_accuracies,ada_accuracies,nn_accuracies)

mean_accuracies <- colMeans(accuracies)
mean_accuracies_df <- as.data.frame(mean_accuracies,titles)
mean_accuracies_df$titles <-t(titles)
mean_accuracies_df_sorted <- mean_accuracies_df[order(-mean_accuracies),,drop=FALSE]

Overall_accuracy_bar <- ggplot(mean_accuracies_df_sorted, aes(x=reorder(titles,mean_accuracies),y=mean_accuracies)) +
  geom_bar(stat="identity",fill='gray90') +
  coord_flip() + 
  geom_text(aes(label=scales::percent(mean_accuracies/100), hjust=1.2), size=8) +
  geom_text(aes(y=0,label=reorder(titles,mean_accuracies)), hjust=0, size=10) +
  labs(title = "Overall accuracy", x = "", y = "") + 
  labs(x=NULL)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=20)
        )
Overall_accuracy_bar
png(filename="overall_accuracy.png",width=720,height=480)
multiplot(Overall_accuracy_bar)
dev.off()

# balance accuracy ----
svme_balancedacc1 <- svme_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
svme_balancedacc2 <- svme_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
svme_balancedacc3 <- svme_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
svme_balancedacc4 <- svme_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
svme_balancedacc5 <- svme_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

svmp_balancedacc1 <- svmp_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
svmp_balancedacc2 <- svmp_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
svmp_balancedacc3 <- svmp_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
svmp_balancedacc4 <- svmp_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
svmp_balancedacc5 <- svmp_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

svmr_balancedacc1 <- svmr_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
svmr_balancedacc2 <- svmr_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
svmr_balancedacc3 <- svmr_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
svmr_balancedacc4 <- svmr_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
svmr_balancedacc5 <- svmr_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

svmke_balancedacc1 <- svmke_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
svmke_balancedacc2 <- svmke_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
svmke_balancedacc3 <- svmke_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
svmke_balancedacc4 <- svmke_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
svmke_balancedacc5 <- svmke_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

svmkp_balancedacc1 <- svmkp_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
svmkp_balancedacc2 <- svmkp_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
svmkp_balancedacc3 <- svmkp_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
svmkp_balancedacc4 <- svmkp_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
svmkp_balancedacc5 <- svmkp_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

svmkr_balancedacc1 <- svmkr_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
svmkr_balancedacc2 <- svmkr_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
svmkr_balancedacc3 <- svmkr_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
svmkr_balancedacc4 <- svmkr_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
svmkr_balancedacc5 <- svmkr_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

rndfor_balancedacc1 <- rndfor_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
rndfor_balancedacc2 <- rndfor_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
rndfor_balancedacc3 <- rndfor_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
rndfor_balancedacc4 <- rndfor_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
rndfor_balancedacc5 <- rndfor_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

ada_balancedacc1 <- ada_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
ada_balancedacc2 <- ada_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
ada_balancedacc3 <- ada_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
ada_balancedacc4 <- ada_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
ada_balancedacc5 <- ada_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

nn_balancedacc1 <- nn_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
nn_balancedacc2 <- nn_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
nn_balancedacc3 <- nn_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
nn_balancedacc4 <- nn_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
nn_balancedacc5 <- nn_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

knn_balancedacc1 <- knn_confusionmatrix1[["byClass"]][["Balanced Accuracy"]]
knn_balancedacc2 <- knn_confusionmatrix2[["byClass"]][["Balanced Accuracy"]]
knn_balancedacc3 <- knn_confusionmatrix3[["byClass"]][["Balanced Accuracy"]]
knn_balancedacc4 <- knn_confusionmatrix4[["byClass"]][["Balanced Accuracy"]]
knn_balancedacc5 <- knn_confusionmatrix5[["byClass"]][["Balanced Accuracy"]]

NB_balancedacc1 <- NBconfusionmatrix1[["byClass"]][["Balanced Accuracy"]]
NB_balancedacc2 <- NBconfusionmatrix2[["byClass"]][["Balanced Accuracy"]]
NB_balancedacc3 <- NBconfusionmatrix3[["byClass"]][["Balanced Accuracy"]]
NB_balancedacc4 <- NBconfusionmatrix4[["byClass"]][["Balanced Accuracy"]]
NB_balancedacc5 <- NBconfusionmatrix5[["byClass"]][["Balanced Accuracy"]]

lda_balancedacc1 <- ldaconfusionmatrix1[["byClass"]][["Balanced Accuracy"]]
lda_balancedacc2 <- ldaconfusionmatrix2[["byClass"]][["Balanced Accuracy"]]
lda_balancedacc3 <- ldaconfusionmatrix3[["byClass"]][["Balanced Accuracy"]]
lda_balancedacc4 <- ldaconfusionmatrix4[["byClass"]][["Balanced Accuracy"]]
lda_balancedacc5 <- ldaconfusionmatrix5[["byClass"]][["Balanced Accuracy"]]

qda_balancedacc1 <- qdaconfusionmatrix1[["byClass"]][["Balanced Accuracy"]]
qda_balancedacc2 <- qdaconfusionmatrix2[["byClass"]][["Balanced Accuracy"]]
qda_balancedacc3 <- qdaconfusionmatrix3[["byClass"]][["Balanced Accuracy"]]
qda_balancedacc4 <- qdaconfusionmatrix4[["byClass"]][["Balanced Accuracy"]]
qda_balancedacc5 <- qdaconfusionmatrix5[["byClass"]][["Balanced Accuracy"]]

svme_balancedacc <- c(svme_balancedacc1,svme_balancedacc2,svme_balancedacc3,svme_balancedacc4,svme_balancedacc5)
svmp_balancedacc <- c(svmp_balancedacc1,svmp_balancedacc2,svmp_balancedacc3,svmp_balancedacc4,svmp_balancedacc5)
svmr_balancedacc <- c(svmr_balancedacc1,svmr_balancedacc2,svmr_balancedacc3,svmr_balancedacc4,svmr_balancedacc5)
svmke_balancedacc <- c(svmke_balancedacc1,svmke_balancedacc2,svmke_balancedacc3,svmke_balancedacc4,svmke_balancedacc5)
svmkp_balancedacc <- c(svmkp_balancedacc1,svmkp_balancedacc2,svmkp_balancedacc3,svmkp_balancedacc4,svmkp_balancedacc5)
svmkr_balancedacc <- c(svmkr_balancedacc1,svmkr_balancedacc2,svmkr_balancedacc3,svmkr_balancedacc4,svmkr_balancedacc5)
rndfor_balancedacc <- c(rndfor_balancedacc1,rndfor_balancedacc2,rndfor_balancedacc3,rndfor_balancedacc4,rndfor_balancedacc5)
ada_balancedacc <- c(ada_balancedacc1,ada_balancedacc2,ada_balancedacc3,ada_balancedacc4,ada_balancedacc5)
nn_balancedacc <- c(nn_balancedacc1,nn_balancedacc2,nn_balancedacc3,nn_balancedacc4,nn_balancedacc5)
knn_balancedacc <- c(knn_balancedacc1,knn_balancedacc2,knn_balancedacc3,knn_balancedacc4,knn_balancedacc5)
NB_balancedacc <- c(NB_balancedacc1,NB_balancedacc2,NB_balancedacc3,NB_balancedacc4,NB_balancedacc5)
lda_balancedacc <- c(lda_balancedacc1,lda_balancedacc2,lda_balancedacc3,lda_balancedacc4,lda_balancedacc5)
qda_balancedacc <- c(qda_balancedacc1,qda_balancedacc2,qda_balancedacc3,qda_balancedacc4,qda_balancedacc5)

svme_balancedacc
svmp_balancedacc
svmr_balancedacc
svmke_balancedacc
svmkp_balancedacc
svmkr_balancedacc
rndfor_balancedacc
ada_balancedacc
nn_balancedacc
knn_balancedacc
NB_balancedacc
lda_balancedacc
qda_balancedacc

svme_balancedacc_avg <- mean(svme_balancedacc)
svmp_balancedacc_avg <- mean(svmp_balancedacc)
svmr_balancedacc_avg <- mean(svmr_balancedacc)
svmke_balancedacc_avg <- mean(svmke_balancedacc)
svmkp_balancedacc_avg <- mean(svmkp_balancedacc)
svmkr_balancedacc_avg <- mean(svmkr_balancedacc)
rndfor_balancedacc_avg <-mean(rndfor_balancedacc)
ada_balancedacc_avg <- mean(ada_balancedacc)
nn_balancedacc_avg <- mean(nn_balancedacc)
knn_balancedacc_avg <- mean(knn_balancedacc)
NB_balancedacc_avg <- mean(NB_balancedacc)
lda_balancedacc_avg <- mean(lda_balancedacc)
qda_balancedacc_avg <- mean(qda_balancedacc)

svme_balancedacc_avg 
svmp_balancedacc_avg
svmr_balancedacc_avg
svmke_balancedacc_avg
svmkp_balancedacc_avg
svmkr_balancedacc_avg
rndfor_balancedacc_avg 
ada_balancedacc_avg 
nn_balancedacc_avg 
knn_balancedacc_avg
NB_balancedacc_avg
lda_balancedacc_avg
qda_balancedacc_avg


# Comparing sensitivity


# sensitivity ----
svme_sensitivity1 <- svme_confusionmatrix1[["byClass"]][["Sensitivity"]]
svme_sensitivity2 <- svme_confusionmatrix2[["byClass"]][["Sensitivity"]]
svme_sensitivity3 <- svme_confusionmatrix3[["byClass"]][["Sensitivity"]]
svme_sensitivity4 <- svme_confusionmatrix4[["byClass"]][["Sensitivity"]]
svme_sensitivity5 <- svme_confusionmatrix5[["byClass"]][["Sensitivity"]]

svmp_sensitivity1 <- svmp_confusionmatrix1[["byClass"]][["Sensitivity"]]
svmp_sensitivity2 <- svmp_confusionmatrix2[["byClass"]][["Sensitivity"]]
svmp_sensitivity3 <- svmp_confusionmatrix3[["byClass"]][["Sensitivity"]]
svmp_sensitivity4 <- svmp_confusionmatrix4[["byClass"]][["Sensitivity"]]
svmp_sensitivity5 <- svmp_confusionmatrix5[["byClass"]][["Sensitivity"]]

svmr_sensitivity1 <- svmr_confusionmatrix1[["byClass"]][["Sensitivity"]]
svmr_sensitivity2 <- svmr_confusionmatrix2[["byClass"]][["Sensitivity"]]
svmr_sensitivity3 <- svmr_confusionmatrix3[["byClass"]][["Sensitivity"]]
svmr_sensitivity4 <- svmr_confusionmatrix4[["byClass"]][["Sensitivity"]]
svmr_sensitivity5 <- svmr_confusionmatrix5[["byClass"]][["Sensitivity"]]

svmke_sensitivity1 <- svmke_confusionmatrix1[["byClass"]][["Sensitivity"]]
svmke_sensitivity2 <- svmke_confusionmatrix2[["byClass"]][["Sensitivity"]]
svmke_sensitivity3 <- svmke_confusionmatrix3[["byClass"]][["Sensitivity"]]
svmke_sensitivity4 <- svmke_confusionmatrix4[["byClass"]][["Sensitivity"]]
svmke_sensitivity5 <- svmke_confusionmatrix5[["byClass"]][["Sensitivity"]]

svmkp_sensitivity1 <- svmkp_confusionmatrix1[["byClass"]][["Sensitivity"]]
svmkp_sensitivity2 <- svmkp_confusionmatrix2[["byClass"]][["Sensitivity"]]
svmkp_sensitivity3 <- svmkp_confusionmatrix3[["byClass"]][["Sensitivity"]]
svmkp_sensitivity4 <- svmkp_confusionmatrix4[["byClass"]][["Sensitivity"]]
svmkp_sensitivity5 <- svmkp_confusionmatrix5[["byClass"]][["Sensitivity"]]

svmkr_sensitivity1 <- svmkr_confusionmatrix1[["byClass"]][["Sensitivity"]]
svmkr_sensitivity2 <- svmkr_confusionmatrix2[["byClass"]][["Sensitivity"]]
svmkr_sensitivity3 <- svmkr_confusionmatrix3[["byClass"]][["Sensitivity"]]
svmkr_sensitivity4 <- svmkr_confusionmatrix4[["byClass"]][["Sensitivity"]]
svmkr_sensitivity5 <- svmkr_confusionmatrix5[["byClass"]][["Sensitivity"]]

rndfor_sensitivity1 <- rndfor_confusionmatrix1[["byClass"]][["Sensitivity"]]
rndfor_sensitivity2 <- rndfor_confusionmatrix2[["byClass"]][["Sensitivity"]]
rndfor_sensitivity3 <- rndfor_confusionmatrix3[["byClass"]][["Sensitivity"]]
rndfor_sensitivity4 <- rndfor_confusionmatrix4[["byClass"]][["Sensitivity"]]
rndfor_sensitivity5 <- rndfor_confusionmatrix5[["byClass"]][["Sensitivity"]]

ada_sensitivity1 <- ada_confusionmatrix1[["byClass"]][["Sensitivity"]]
ada_sensitivity2 <- ada_confusionmatrix2[["byClass"]][["Sensitivity"]]
ada_sensitivity3 <- ada_confusionmatrix3[["byClass"]][["Sensitivity"]]
ada_sensitivity4 <- ada_confusionmatrix4[["byClass"]][["Sensitivity"]]
ada_sensitivity5 <- ada_confusionmatrix5[["byClass"]][["Sensitivity"]]

nn_sensitivity1 <- nn_confusionmatrix1[["byClass"]][["Sensitivity"]]
nn_sensitivity2 <- nn_confusionmatrix2[["byClass"]][["Sensitivity"]]
nn_sensitivity3 <- nn_confusionmatrix3[["byClass"]][["Sensitivity"]]
nn_sensitivity4 <- nn_confusionmatrix4[["byClass"]][["Sensitivity"]]
nn_sensitivity5 <- nn_confusionmatrix5[["byClass"]][["Sensitivity"]]

knn_sensitivity1 <- knn_confusionmatrix1[["byClass"]][["Sensitivity"]]
knn_sensitivity2 <- knn_confusionmatrix2[["byClass"]][["Sensitivity"]]
knn_sensitivity3 <- knn_confusionmatrix3[["byClass"]][["Sensitivity"]]
knn_sensitivity4 <- knn_confusionmatrix4[["byClass"]][["Sensitivity"]]
knn_sensitivity5 <- knn_confusionmatrix5[["byClass"]][["Sensitivity"]]

NB_sensitivity1 <- NBconfusionmatrix1[["byClass"]][["Sensitivity"]]
NB_sensitivity2 <- NBconfusionmatrix2[["byClass"]][["Sensitivity"]]
NB_sensitivity3 <- NBconfusionmatrix3[["byClass"]][["Sensitivity"]]
NB_sensitivity4 <- NBconfusionmatrix4[["byClass"]][["Sensitivity"]]
NB_sensitivity5 <- NBconfusionmatrix5[["byClass"]][["Sensitivity"]]

lda_sensitivity1 <- ldaconfusionmatrix1[["byClass"]][["Sensitivity"]]
lda_sensitivity2 <- ldaconfusionmatrix2[["byClass"]][["Sensitivity"]]
lda_sensitivity3 <- ldaconfusionmatrix3[["byClass"]][["Sensitivity"]]
lda_sensitivity4 <- ldaconfusionmatrix4[["byClass"]][["Sensitivity"]]
lda_sensitivity5 <- ldaconfusionmatrix5[["byClass"]][["Sensitivity"]]

qda_sensitivity1 <- qdaconfusionmatrix1[["byClass"]][["Sensitivity"]]
qda_sensitivity2 <- qdaconfusionmatrix2[["byClass"]][["Sensitivity"]]
qda_sensitivity3 <- qdaconfusionmatrix3[["byClass"]][["Sensitivity"]]
qda_sensitivity4 <- qdaconfusionmatrix4[["byClass"]][["Sensitivity"]]
qda_sensitivity5 <- qdaconfusionmatrix5[["byClass"]][["Sensitivity"]]

svme_sensitivity <- c(svme_sensitivity1,svme_sensitivity2,svme_sensitivity3,svme_sensitivity4,svme_sensitivity5)
svmp_sensitivity <- c(svmp_sensitivity1,svmp_sensitivity2,svmp_sensitivity3,svmp_sensitivity4,svmp_sensitivity5)
svmr_sensitivity <- c(svmr_sensitivity1,svmr_sensitivity2,svmr_sensitivity3,svmr_sensitivity4,svmr_sensitivity5)
svmke_sensitivity <- c(svmke_sensitivity1,svmke_sensitivity2,svmke_sensitivity3,svmke_sensitivity4,svmke_sensitivity5)
svmkp_sensitivity <- c(svmkp_sensitivity1,svmkp_sensitivity2,svmkp_sensitivity3,svmkp_sensitivity4,svmkp_sensitivity5)
svmkr_sensitivity <- c(svmkr_sensitivity1,svmkr_sensitivity2,svmkr_sensitivity3,svmkr_sensitivity4,svmkr_sensitivity5)
rndfor_sensitivity <- c(rndfor_sensitivity1,rndfor_sensitivity2,rndfor_sensitivity3,rndfor_sensitivity4,rndfor_sensitivity5)
ada_sensitivity <- c(ada_sensitivity1,ada_sensitivity2,ada_sensitivity3,ada_sensitivity4,ada_sensitivity5)
nn_sensitivity <- c(nn_sensitivity1,nn_sensitivity2,nn_sensitivity3,nn_sensitivity4,nn_sensitivity5)
knn_sensitivity <- c(knn_sensitivity1,knn_sensitivity2,knn_sensitivity3,knn_sensitivity4,knn_sensitivity5)
NB_sensitivity <- c(NB_sensitivity1,NB_sensitivity2,NB_sensitivity3,NB_sensitivity4,NB_sensitivity5)
lda_sensitivity <- c(lda_sensitivity1,lda_sensitivity2,lda_sensitivity3,lda_sensitivity4,lda_sensitivity5)
qda_sensitivity <- c(qda_sensitivity1,qda_sensitivity2,qda_sensitivity3,qda_sensitivity4,qda_sensitivity5)

svme_sensitivity
svmp_sensitivity
svmr_sensitivity
svmke_sensitivity
svmkp_sensitivity
svmkr_sensitivity
rndfor_sensitivity
ada_sensitivity
nn_sensitivity
knn_sensitivity
NB_sensitivity
lda_sensitivity
qda_sensitivity

svme_sensitivity_avg <- mean(svme_sensitivity)
svmp_sensitivity_avg <- mean(svmp_sensitivity)
svmr_sensitivity_avg <- mean(svmr_sensitivity)
svmke_sensitivity_avg <- mean(svmke_sensitivity)
svmkp_sensitivity_avg <- mean(svmkp_sensitivity)
svmkr_sensitivity_avg <- mean(svmkr_sensitivity)
rndfor_sensitivity_avg <-mean(rndfor_sensitivity)
ada_sensitivity_avg <- mean(ada_sensitivity)
nn_sensitivity_avg <- mean(nn_sensitivity)
knn_sensitivity_avg <- mean(knn_sensitivity)
NB_sensitivity_avg <- mean(NB_sensitivity)
lda_sensitivity_avg <- mean(lda_sensitivity)
qda_sensitivity_avg <- mean(qda_sensitivity)

svme_sensitivity_avg 
svmp_sensitivity_avg
svmr_sensitivity_avg
svmke_sensitivity_avg
svmkp_sensitivity_avg
svmkr_sensitivity_avg
rndfor_sensitivity_avg 
ada_sensitivity_avg 
nn_sensitivity_avg 
knn_sensitivity_avg
NB_sensitivity_avg
lda_sensitivity_avg
qda_sensitivity_avg

# specificity ----
svme_specificity1 <- svme_confusionmatrix1[["byClass"]][["Specificity"]]
svme_specificity2 <- svme_confusionmatrix2[["byClass"]][["Specificity"]]
svme_specificity3 <- svme_confusionmatrix3[["byClass"]][["Specificity"]]
svme_specificity4 <- svme_confusionmatrix4[["byClass"]][["Specificity"]]
svme_specificity5 <- svme_confusionmatrix5[["byClass"]][["Specificity"]]

svmp_specificity1 <- svmp_confusionmatrix1[["byClass"]][["Specificity"]]
svmp_specificity2 <- svmp_confusionmatrix2[["byClass"]][["Specificity"]]
svmp_specificity3 <- svmp_confusionmatrix3[["byClass"]][["Specificity"]]
svmp_specificity4 <- svmp_confusionmatrix4[["byClass"]][["Specificity"]]
svmp_specificity5 <- svmp_confusionmatrix5[["byClass"]][["Specificity"]]

svmr_specificity1 <- svmr_confusionmatrix1[["byClass"]][["Specificity"]]
svmr_specificity2 <- svmr_confusionmatrix2[["byClass"]][["Specificity"]]
svmr_specificity3 <- svmr_confusionmatrix3[["byClass"]][["Specificity"]]
svmr_specificity4 <- svmr_confusionmatrix4[["byClass"]][["Specificity"]]
svmr_specificity5 <- svmr_confusionmatrix5[["byClass"]][["Specificity"]]

svmke_specificity1 <- svmke_confusionmatrix1[["byClass"]][["Specificity"]]
svmke_specificity2 <- svmke_confusionmatrix2[["byClass"]][["Specificity"]]
svmke_specificity3 <- svmke_confusionmatrix3[["byClass"]][["Specificity"]]
svmke_specificity4 <- svmke_confusionmatrix4[["byClass"]][["Specificity"]]
svmke_specificity5 <- svmke_confusionmatrix5[["byClass"]][["Specificity"]]

svmkp_specificity1 <- svmkp_confusionmatrix1[["byClass"]][["Specificity"]]
svmkp_specificity2 <- svmkp_confusionmatrix2[["byClass"]][["Specificity"]]
svmkp_specificity3 <- svmkp_confusionmatrix3[["byClass"]][["Specificity"]]
svmkp_specificity4 <- svmkp_confusionmatrix4[["byClass"]][["Specificity"]]
svmkp_specificity5 <- svmkp_confusionmatrix5[["byClass"]][["Specificity"]]

svmkr_specificity1 <- svmkr_confusionmatrix1[["byClass"]][["Specificity"]]
svmkr_specificity2 <- svmkr_confusionmatrix2[["byClass"]][["Specificity"]]
svmkr_specificity3 <- svmkr_confusionmatrix3[["byClass"]][["Specificity"]]
svmkr_specificity4 <- svmkr_confusionmatrix4[["byClass"]][["Specificity"]]
svmkr_specificity5 <- svmkr_confusionmatrix5[["byClass"]][["Specificity"]]

rndfor_specificity1 <- rndfor_confusionmatrix1[["byClass"]][["Specificity"]]
rndfor_specificity2 <- rndfor_confusionmatrix2[["byClass"]][["Specificity"]]
rndfor_specificity3 <- rndfor_confusionmatrix3[["byClass"]][["Specificity"]]
rndfor_specificity4 <- rndfor_confusionmatrix4[["byClass"]][["Specificity"]]
rndfor_specificity5 <- rndfor_confusionmatrix5[["byClass"]][["Specificity"]]

ada_specificity1 <- ada_confusionmatrix1[["byClass"]][["Specificity"]]
ada_specificity2 <- ada_confusionmatrix2[["byClass"]][["Specificity"]]
ada_specificity3 <- ada_confusionmatrix3[["byClass"]][["Specificity"]]
ada_specificity4 <- ada_confusionmatrix4[["byClass"]][["Specificity"]]
ada_specificity5 <- ada_confusionmatrix5[["byClass"]][["Specificity"]]

nn_specificity1 <- nn_confusionmatrix1[["byClass"]][["Specificity"]]
nn_specificity2 <- nn_confusionmatrix2[["byClass"]][["Specificity"]]
nn_specificity3 <- nn_confusionmatrix3[["byClass"]][["Specificity"]]
nn_specificity4 <- nn_confusionmatrix4[["byClass"]][["Specificity"]]
nn_specificity5 <- nn_confusionmatrix5[["byClass"]][["Specificity"]]

knn_specificity1 <- knn_confusionmatrix1[["byClass"]][["Specificity"]]
knn_specificity2 <- knn_confusionmatrix2[["byClass"]][["Specificity"]]
knn_specificity3 <- knn_confusionmatrix3[["byClass"]][["Specificity"]]
knn_specificity4 <- knn_confusionmatrix4[["byClass"]][["Specificity"]]
knn_specificity5 <- knn_confusionmatrix5[["byClass"]][["Specificity"]]

NB_specificity1 <- NBconfusionmatrix1[["byClass"]][["Specificity"]]
NB_specificity2 <- NBconfusionmatrix2[["byClass"]][["Specificity"]]
NB_specificity3 <- NBconfusionmatrix3[["byClass"]][["Specificity"]]
NB_specificity4 <- NBconfusionmatrix4[["byClass"]][["Specificity"]]
NB_specificity5 <- NBconfusionmatrix5[["byClass"]][["Specificity"]]

lda_specificity1 <- ldaconfusionmatrix1[["byClass"]][["Specificity"]]
lda_specificity2 <- ldaconfusionmatrix2[["byClass"]][["Specificity"]]
lda_specificity3 <- ldaconfusionmatrix3[["byClass"]][["Specificity"]]
lda_specificity4 <- ldaconfusionmatrix4[["byClass"]][["Specificity"]]
lda_specificity5 <- ldaconfusionmatrix5[["byClass"]][["Specificity"]]

qda_specificity1 <- qdaconfusionmatrix1[["byClass"]][["Specificity"]]
qda_specificity2 <- qdaconfusionmatrix2[["byClass"]][["Specificity"]]
qda_specificity3 <- qdaconfusionmatrix3[["byClass"]][["Specificity"]]
qda_specificity4 <- qdaconfusionmatrix4[["byClass"]][["Specificity"]]
qda_specificity5 <- qdaconfusionmatrix5[["byClass"]][["Specificity"]]

svme_specificity <- c(svme_specificity1,svme_specificity2,svme_specificity3,svme_specificity4,svme_specificity5)
svmp_specificity <- c(svmp_specificity1,svmp_specificity2,svmp_specificity3,svmp_specificity4,svmp_specificity5)
svmr_specificity <- c(svmr_specificity1,svmr_specificity2,svmr_specificity3,svmr_specificity4,svmr_specificity5)
svmke_specificity <- c(svmke_specificity1,svmke_specificity2,svmke_specificity3,svmke_specificity4,svmke_specificity5)
svmkp_specificity <- c(svmkp_specificity1,svmkp_specificity2,svmkp_specificity3,svmkp_specificity4,svmkp_specificity5)
svmkr_specificity <- c(svmkr_specificity1,svmkr_specificity2,svmkr_specificity3,svmkr_specificity4,svmkr_specificity5)
rndfor_specificity <- c(rndfor_specificity1,rndfor_specificity2,rndfor_specificity3,rndfor_specificity4,rndfor_specificity5)
ada_specificity <- c(ada_specificity1,ada_specificity2,ada_specificity3,ada_specificity4,ada_specificity5)
nn_specificity <- c(nn_specificity1,nn_specificity2,nn_specificity3,nn_specificity4,nn_specificity5)
knn_specificity <- c(knn_specificity1,knn_specificity2,knn_specificity3,knn_specificity4,knn_specificity5)
NB_specificity <- c(NB_specificity1,NB_specificity2,NB_specificity3,NB_specificity4,NB_specificity5)
lda_specificity <- c(lda_specificity1,lda_specificity2,lda_specificity3,lda_specificity4,lda_specificity5)
qda_specificity <- c(qda_specificity1,qda_specificity2,qda_specificity3,qda_specificity4,qda_specificity5)

svme_specificity
svmp_specificity
svmr_specificity
svmke_specificity
svmkp_specificity
svmkr_specificity
rndfor_specificity
ada_specificity
nn_specificity
knn_specificity
NB_specificity
lda_specificity
qda_specificity

svme_specificity_avg <- mean(svme_specificity)
svmp_specificity_avg <- mean(svmp_specificity)
svmr_specificity_avg <- mean(svmr_specificity)
svmke_specificity_avg <- mean(svmke_specificity)
svmkp_specificity_avg <- mean(svmkp_specificity)
svmkr_specificity_avg <- mean(svmkr_specificity)
rndfor_specificity_avg <-mean(rndfor_specificity)
ada_specificity_avg <- mean(ada_specificity)
nn_specificity_avg <- mean(nn_specificity)
knn_specificity_avg <- mean(knn_specificity)
NB_specificity_avg <- mean(NB_specificity)
lda_specificity_avg <- mean(lda_specificity)
qda_specificity_avg <- mean(qda_specificity)

svme_specificity_avg 
svmp_specificity_avg
svmr_specificity_avg
svmke_specificity_avg
svmkp_specificity_avg
svmkr_specificity_avg
rndfor_specificity_avg 
ada_specificity_avg 
nn_specificity_avg 
knn_specificity_avg
NB_specificity_avg
lda_specificity_avg
qda_specificity_avg

# Comparing specificity ----
# titles <- cbind('knn','Naive Bayes','linear discriminant analysis','quadratic discriminant analysis','linear SVM','quadratic SVM','radial SVM','tuned linear SVM','tuned quadratic SVM','tuned radial SVM','random forest','boosted gradient','neural network')

specificities <- cbind(knn_specificity,NB_specificity,lda_specificity,qda_specificity,svme_specificity,svmp_specificity,svmr_specificity,svmke_specificity,svmkp_specificity,svmkr_specificity,rndfor_specificity,ada_specificity,nn_specificity)

mean_specificities <- colMeans(specificities)

mean_specificities_df <- as.data.frame(mean_specificities,titles)
mean_specificities_df$titles <- t(titles)
mean_specificities_df_sorted <- mean_specificities_df[order(-mean_specificities),,drop=FALSE]

mean_specificities_plot <- ggplot(mean_specificities_df_sorted, aes(x=reorder(titles,mean_specificities), y=mean_specificities)) +
  geom_bar(stat="identity", fill='gray90') +
  coord_flip() + 
  geom_text(aes(label=scales::percent(round(mean_specificities,2)), hjust=1.2), size=8) +
  geom_text(aes(y=0,label=reorder(titles,mean_specificities)), hjust=0, size=10) +
  labs(title = "Mean of malignant misdiagnosed as benign", x = "", y = "") + 
  labs(x=NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=20)
  )
mean_specificities_plot

png(filename="mean_specificity.png",width=720,height=480)
multiplot(mean_specificities_plot)
dev.off()

# Scatterplots QDA & svme ----

scatter_validation1 <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(diagnosis))) + 
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.direction = "horizontal"
        )
scatter_validation1

svme_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(svme_pred1))) + 
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.direction = "horizontal"
  )
svme_scatter_1_pred

qda_scatter_1_pred <- ggplot(validation1, aes(area_worst, concave.points_worst, colour = factor(qdapred1$class))) + 
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.direction = "horizontal"
  )
qda_scatter_1_pred

png(filename="scatterplotsquadsvm&linearSVM.png",width=720,height=480)
multiplot(scatter_validation1,scatter_validation1,qda_scatter_1_pred,svme_scatter_1_pred,cols = 2)
dev.off()