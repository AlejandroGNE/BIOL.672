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

# Cross-validation ----
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

# Accuracies ----

knn_accuracies <- c(knn_accuracy1,knn_accuracy2,knn_accuracy3,knn_accuracy4,knn_accuracy5)
NB_accuracies <- c(NB_accuracy1,NB_accuracy2,NB_accuracy3,NB_accuracy4,NB_accuracy5)
lda_accuracies <- c(lda_accuracy1,lda_accuracy2,lda_accuracy3,lda_accuracy4,lda_accuracy5)
qda_accuracies <- c(qda_accuracy1,qda_accuracy2,qda_accuracy3,qda_accuracy4,qda_accuracy5)

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

# export plots ----

pdf(file = 'ElenesAGN_A2_GraphicOutputs.pdf')

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

dev.off()

png(filename="Pairs_correlograms.png",width=3840,height=2160)
multiplot(corrs)
dev.off()

# export confusion matrices ----
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

sink(file = 'ElenesAGN_A2_TextOutputs.txt')
print(knn_confusionmatrix1)
print(knn_confusionmatrix2)
print(knn_confusionmatrix3)
print(knn_confusionmatrix4)
print(knn_confusionmatrix5)
print(NBconfusionmatrix1)
print(NBconfusionmatrix2)
print(NBconfusionmatrix3)
print(NBconfusionmatrix4)
print(NBconfusionmatrix5)
print(ldaconfusionmatrix1)
print(ldaconfusionmatrix2)
print(ldaconfusionmatrix3)
print(ldaconfusionmatrix4)
print(ldaconfusionmatrix5)
print(qdaconfusionmatrix1)
print(qdaconfusionmatrix2)
print(qdaconfusionmatrix3)
print(qdaconfusionmatrix4)
print(qdaconfusionmatrix5)
sink()

