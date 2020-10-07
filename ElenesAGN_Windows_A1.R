# Alejandro G.N. Elenes
# Problem 1 ----

# Code written in Windows 10 Home

# R packages needed to run this script:
library(ggplot2)
library(psych)
library(readxl)
library(openxlsx)
library(reshape2)
library(dplyr)
library(tidyverse)
library(rstatix)
library(Rmisc)
library(grid)
library(MASS)
library(mixtools)

# Data files necessary to run this script:

# 3_1_Generator_Y2018.xlsx # included with script, can also be accessed through https://www.eia.gov/electricity/data/eia860/archive/xls/eia8602018.zip

# NEM_CurrentlyInterconnectedDataset_2020-04-30.csv # can be accessed through here: https://www.californiadgstats.ca.gov/download/interconnection_rule21_projects/NEM_CurrentlyInterconnectedDataset_2020-04-30.zip

# Problem 2 ----

# Create a dataset of 5000 random numbers generated in R using any model distribution you like. Write an R script that reports the sample mean and sample standard deviation and then plots a histogram from this data and that includes both a density line laid over the histogram bars and a normal curve fitted to the data and overlaid as well.

random_data <- rnorm(5000) # create variable with 5000 random numbers using normal distribution
describe(random_data) # use describe function from psych package to report mean and st dev
hist(random_data, prob = TRUE) # plot a histogram of random data
lines(density(random_data,na.rm=T),col="red",lwd=2) # over histogram, add density line
curve(dnorm(x, mean=mean(random_data), sd=sd(random_data)), add=TRUE) # over histogram, plot normal distribution based on random data

random <- data.frame(random_data)
hist_random <- ggplot(random, aes(x=random_data)) + geom_histogram(aes(y=..density..),colour='black',fill='gray') + geom_density(color='red') + stat_function(fun=dnorm, args=list(mean = mean(random_data), sd = sd(random_data))) + labs(title='Histogram of random data') + theme_bw() # using ggplot

# Problem 3 ----

# To automate the printing of stat test output in R you will want to explore the 'sink' and 'unsink' functions in R. Use your script to create a file called 'desc.txt' and print your mean and SD from #2 to it. Add lines to your script to rename and save Rplots.pdf as histo.pdf.

sink(file = 'desc.txt') # trigger printing to text file
stats_random <- describe(random_data) 
print(stats_random) # print stats of random data
sink() # stop printing to text file
pdf(file = 'histo.pdf') # trigger printing to pdf file
multiplot(hist_random)
dev.off() # stop printing to pdf

# Problem 4 ----

# 1 - Find an appropriate data set for a one-way ANOVA. Shape the data so you can use the read.table function from an input data file (as demo during class). 
# 2 - Extend your R script so that that it conducts a oneway ANOVA (i.e. oneway.test(values, category)).
# 3 - Add to this script, an error bar chart output of your categories with colored bars of your choice.
# 4 - Lastly add pairwise t tests for each category pair and use both Bonferroni and Benjamini-Hochberg multiple test correction using the p.adjust function in R. 
# 5 - As in #3, export output to appropriately named text and pdf files.
# 6 - Be sure to include verbal interpretation printed to screen and/or output files.

  # 4.1 - about the dataset

# I will look for trends in power plants retirements in the US. e.g:how old are power plants at age of retirement? how much do plants physical life differs from their originally planned lifetime? Is there a trend, and is this trend different among different technologies?

# I will use the dataset from survey form EIA-860, which collects generator-level specific information about existing and planned generators and associated environmental equipment at electric power plants with 1 MW or greater of combined nameplate capacity. Raw data comes in a zip folder with multiple excel files. I use the sheet 'Retired and Canceled' from the file '3_1_Generator_Y2018.xlsx'

  # 4.2 - code
  # 4.2.1 - load data, create variables

Retirement_dataset <- read.xlsx("3_1_Generator_Y2018.xlsx", sheet = "Retired and Canceled", startRow = 2, na.strings = " ")
Retirement_dataset <- Retirement_dataset[complete.cases(Retirement_dataset$Retirement.Year),] # drop missing values (cancelled power plants)
Online <- Retirement_dataset$Operating.Year # Year plant went online / was born
Retired <- Retirement_dataset$Retirement.Year # Year plant went offline / was retired / died
Source <- Retirement_dataset$Energy.Source.1 # Primary fuel
Capacity <- Retirement_dataset$`Nameplate.Capacity.(MW)` # Plant capacity in MW
Age <- Retired - Online # Plant phyisical life / age
Fuel_unique_categories <- unique(Source)

  # 4.2.2 - Group into higher level categories based on fuel type for easier tractability

Fuel_categories <- matrix(0,length(Retired),1)
Fuel_categories[which(Source == "ANT" )] <- "COAL"
Fuel_categories[which(Source == "BIT" )] <- "COAL"
Fuel_categories[which(Source == "LIG" )] <- "COAL"
Fuel_categories[which(Source == "SGC" )] <- "COAL"
Fuel_categories[which(Source == "SUB" )] <- "COAL"
Fuel_categories[which(Source == "WC"  )] <- "COAL"
Fuel_categories[which(Source == "RC"  )] <- "COAL"
Fuel_categories[which(Source == "DFO" )] <- "OIL"
Fuel_categories[which(Source == "JF"  )] <- "OIL"
Fuel_categories[which(Source == "KER" )] <- "OIL"
Fuel_categories[which(Source == "PC"  )] <- "OIL"
Fuel_categories[which(Source == "PG"  )] <- "OIL"
Fuel_categories[which(Source == "RFO" )] <- "OIL"
Fuel_categories[which(Source == "SGP" )] <- "OIL"
Fuel_categories[which(Source == "WO"  )] <- "OIL"
Fuel_categories[which(Source == "BFG" )] <- "GAS"
Fuel_categories[which(Source == "NG"  )] <- "GAS"
Fuel_categories[which(Source == "OG"  )] <- "GAS"
Fuel_categories[which(Source == "AB"  )] <- "BIO"
Fuel_categories[which(Source == "MSW" )] <- "BIO"
Fuel_categories[which(Source == "OBS" )] <- "BIO"
Fuel_categories[which(Source == "WDS" )] <- "BIO"
Fuel_categories[which(Source == "OBL" )] <- "BIO"
Fuel_categories[which(Source == "SLW" )] <- "BIO"
Fuel_categories[which(Source == "BLQ" )] <- "BIO"
Fuel_categories[which(Source == "WDL" )] <- "BIO"
Fuel_categories[which(Source == "LFG" )] <- "BIO"
Fuel_categories[which(Source == "OBG" )] <- "BIO"
Fuel_categories[which(Source == "SUN" )] <- "SOLAR"
Fuel_categories[which(Source == "WND" )] <- "WIND"
Fuel_categories[which(Source == "GEO" )] <- "GEOTHERMAL"
Fuel_categories[which(Source == "WAT" )] <- "HYDRO"
Fuel_categories[which(Source == "NUC" )] <- "NUCLEAR"
Fuel_categories[which(Source == "PUR" )] <- "OTHER"
Fuel_categories[which(Source == "WH"  )] <- "OTHER"
Fuel_categories[which(Source == "TDF" )] <- "OTHER"
Fuel_categories[which(Source == "MWH" )] <- "OTHER"
Fuel_categories[which(Source == "OTH" )] <- "OTHER"

# Went from 31 to 10 categories
  # 4.2.3 - One way ANOVA

age.anova <- oneway.test(Age~Fuel_categories)
age.anova

  # 4.2.4 - Error bar chart

Retirement_dataset$Fuel.Categories <- Fuel_categories
Retirement_dataset$Age <- Age
stats <- summarySE(Retirement_dataset, measurevar = "Age", groupvars = "Fuel_categories")
stats
ErrorBar_AgeVSFuelCat <- ggplot(stats, aes(x=Fuel_categories, y=Age)) + 
     geom_bar(stat="identity", fill="black", width=0.8) +
     geom_errorbar(aes(ymin=Age-se, ymax=Age+se), width=0.2, colour='red') +
    labs(title = "Error bar chart", x = "\n Technology", y = "Physical life [years]\n") +
  theme(axis.text.x=element_text(angle=15))+ theme_bw()
ErrorBar_AgeVSFuelCat # print plot

  # 4.2.5 - Pairwise t-tests

# Bonferroni
age.pairwise.bonf <- pairwise.t.test(Age, Fuel_categories, p.adj = "bonf")
age.pairwise.bonf

# Benjamini-Hochberg
age.pairwise.benjhoch <- pairwise.t.test(Age, Fuel_categories, p.adj = "BH")
age.pairwise.benjhoch

BonfPvalues <- age.pairwise.bonf$p.value # create table to export for exogenous interpretaion
BHPvalues <- age.pairwise.benjhoch$p.value # idem

# Problem 5 ----

# 1 - Use a Kruskal Wallis test applied to your ANOVA input data to examine it without assumptions of normality.
# 2 - Choose two or more of your categories and test the correlation between them using both Pearson and Spearman rank methods.
# 3 - Make scatterplots showing these relations. 
# 4 - Run a sample KS test to test your assumption of normality.
# 5 - Did this assumption hold? Did outcomes of parametric and nonparametric tests appear consistent? Why or why not? 
# 6 - Be sure to include verbal interpretation printed to screen and/or output files.

  # 5.1. - Kruskal Wallis test

age.kruskal <- kruskal.test(Age~Fuel_categories) # p-value < 0.05 = there are significant differences among some of the groups

  # 5.2. - Pearson & Spearman correlation

age.pearson <- cor.test(x=Capacity,y=Age,method="pearson")
age.spearman <- cor.test(x=Capacity,y=Age,method="spearman")

  # 5.3. - Pearson & Spearman correlation scatterplots

Pearson_scatter <- ggplot(Retirement_dataset, aes(x=Capacity,y=Age)) + geom_point() + geom_smooth(method=lm, se=FALSE) + labs (title="Correlation scatterplot",x= "\nCapacity [MW]",y="Physical life [years]")+ theme_bw()
Pearson_scatter # small correlation in aggregated data, many low-capacity observations; uncomment following lines to see colored categories

# Scatter_colorfuels <- ggplot(Retirement_dataset, aes(x=Capacity,y=Age,color=Fuel_categories)) + geom_point() + geom_smooth(method=lm) + labs (x= "\nCapacity [MW]",y="Physical life [years]")
# Scatter_colorfuels # correlation changes with category

  # 5.4. - Kolmogorov-Smirnov (KS) normality test

age.norm <- ks.test(Age,"pnorm") # p-value < 0.05 = data does not follow indicated distribution
size.norm <- ks.test(Capacity,"pnorm")

# Problem 6 ----

# 1 - Run a simple linear regresssion (lm function in R) on the same comparison from your correlations in #5.
# 2 - How does result compare?
# 3 - When is regression appropriate instead of correlation?
# 4 - Be sure to include verbal interpretation printed to screen and/or outputfiles. Be sure to plot your results. 
# 5 - IMPORTANT: from this point forward and throughout the rest of the course, always use the ggplot2 package to create your plots.

YageXcapacity <- lm(Age~Capacity, data=Retirement_dataset)
OLStest <- summary(YageXcapacity)
ggplotRegression <- function(fit) { 
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
  
  }
OLSplot <- ggplotRegression(YageXcapacity)

# Problem 7 ----

# 1 - Find an appropriate data set from your field (or use your own thesis data) for a series of multivariate statistical tests, similar to the Iris data set demo in class. It should linearly include 4 to 6 measured quantitative variables that are intercorrelated and 1 categorical variable.
# 2 - IMPORTANT: shape the data so you can use read.table function to input data similarly to how I demonstrate this in class.
# 3 - NOTE: you may create a second R script for this latter part of the Unit Assignment that is importing and analyzing your own data.

# download database from the following link and move it into same folder as this script: https://www.californiadgstats.ca.gov/download/interconnection_rule21_projects/NEM_CurrentlyInterconnectedDataset_2020-04-30.zip 
PVdata <- read.csv("NEM_CurrentlyInterconnectedDataset_2020-04-30.csv",header=TRUE)
summary(PVdata) # check numerical variables, NA's

# Balance dataset, drop zero values
DATAPV <- PVdata[complete.cases(PVdata$Total.System.Cost),]
DATAPV$Total.System.Cost[DATAPV$Total.System.Cost==0] <- NA
DATAPV$System.Size.AC[DATAPV$System.Size.AC<0] <- NA
DATAPV$System.Size.AC[DATAPV$System.Size.AC==0] <- NA
DATAPV <- DATAPV[complete.cases(DATAPV$Total.System.Cost),]
DATAPV <- DATAPV[complete.cases(DATAPV$System.Size.AC),]
DATAPV <- DATAPV[complete.cases(DATAPV$System.Size.DC),]
DATAPV <- DATAPV[complete.cases(DATAPV$Electric.Vehicle.Count),]
DATAPV <- DATAPV[complete.cases(DATAPV$Itc.Cost.Basis),]

# quantitative variables
cost <- DATAPV$Total.System.Cost
systemsizedc <- DATAPV$System.Size.DC
systemsizeac <- DATAPV$System.Size.AC
taxcredit <- DATAPV$Itc.Cost.Basis
EVcount <- DATAPV$Electric.Vehicle.Count

# qualitative variable
self <- DATAPV$Self.Installer
dataset <- data.frame(cost,systemsizeac,systemsizedc,taxcredit,EVcount,self)

# Problem 8 ----

# 1 - Use the 'manova' function in R to add a MANOVA to this script using the command 'cbind' to combine the measurements (the dependent variables) into a single test of significance across categories.
# 2 - Be sure to include verbal interpretation printed to screen and/or output files.

manovaY <- cbind(cost,systemsizeac,systemsizedc,taxcredit,EVcount)
manovatest <- manova(manovaY ~ self)
manova <- summary(manovatest)

# Problem 9 ----

# 1 - Utilize the lm function in R to conduct a multiple regression that tells how one variable is predicted by the remaining variables.
# 2 - In your interpretation state: which one is the best predictor?
# 3 - Conduct the same test again but within one of your categories instead of across all of them.
# 4 - Be sure to include verbal interpretation printed to screen and/or output files.

mlr <- lm(cost ~ systemsizeac + taxcredit + EVcount)
stats_mlr <- summary(mlr)

mlr2 <- lm(cost ~ systemsizeac + taxcredit + EVcount, data=subset(dataset,self=='Yes'))
stats_mlr2 <- summary(mlr2)

# Problem 10 ----

# 1 - Create a composite variable (e.g. like 'size' in Iris dataset as = pl*pw + sl*sw) that makes sense to you using your own data set. Rations combining variables are also indicative of many features
# 2 - Use the aov function in R to set up an ANCOVA that compares the relation of one variable predicted by another while controlling for your composite feature.
# 3 - Be sure to include verbal interpretation printed to screen and/or output files.

unitcost <- (cost/systemsizeac)
ancova <- aov(cost ~ systemsizeac*unitcost)
stats_ancova <- summary(ancova)

# Use the following R functions princomp() and factanal()

# Problem 11 ----

# 1 - Conduct a principal components analysis of the individual measurements and output that covers the following questions.
# 2 - NOTE: do not include categorical variables in your dataframe.
# 3 - How successful was the reduction of the data? (explain using the scree plot as a reference)
# 4 - How would you interpret the loadings of the original measurements on each of the PC's?
# 5 - Which one best describes overall size variation (i.e. PC with all positive loadings)?
# 6 - Do any of the PC's appear strongly driven by one particular measurement (see loadings on PC)?

dataset2 <- data.frame(cost,systemsizeac,taxcredit,EVcount)

pca <- prcomp(dataset2, scale=TRUE)
stats_pca <- summary(pca)

pcadata <-
  data.frame(
    PC = paste0("PC", 1:4),
    variance_explained = round((pca$sdev) ^ 2 / sum((pca$sdev) ^ 2) * 100))

pc <- 
  data.frame(
    pc1 = pca$x[, 1],
    pc2 = pca$x[, 2],
    pc3 = pca$x[, 3],
    pc4 = pca$x[, 4]
  )

pca_screeplot <-   ggplot(pcadata, aes(x=PC,y=variance_explained,group=1))+
  geom_bar(stat="identity", fill='black')+
  labs(title ="Scree plot of PCA", x = "\n Principal Component", y = "Percent of variation explained\n")+ theme_bw()

pca_rotations <- (pca[["rotation"]])
pca_rotations

pca_1vs2 <- ggplot(pc, aes(x = pc1, y = pc2)) + geom_point() +
  labs(title = "Principal components", x = "\n PC1", y = "PC2\n") + theme_bw()
pca_2vs3 <- ggplot(pc, aes(x = pc2, y = pc3)) + geom_point() +
  labs(title = "Principal components", x = "\n PC2", y = "PC3\n") + theme_bw()
pca_3vs4 <- ggplot(pc, aes(x = pc3, y = pc4)) + geom_point() +
  labs(title = "Principal components", x = "\n PC3", y = "PC4\n") + theme_bw()

# Problem 12 ----

# 1 - Conduct a Factor Analysis (using same or new data if applicable).
# 2 - Can you determine any common latent underlying 'factors' in your data?
# 3 - How many of such factors are significant?
# 4 - What important characteristics can you find that tend to cluster together?
# 5 - If two traits are far apart along the axis of a significant factor, what does this indicate?
# 6 - What if they are close together?
# 7 - Was this factor analysis particularly successful in identifying a 1 or 2 large underlying latent variables that define something about your data?

factanalysis <- factanal(dataset2,1,rotation="varimax")

# Problem 13 ----

# 1 - Create a scatterplot of two of your most interesting variables (or can use PC 1 against PC 2 from #11 if you prefer) and color the points of the plot according to your main categorical variable.
# 2 - Visually determine how many potential 'clusters' of data points you see in your plot. This is 'k'.
# 3 - Use the kmeans function to run k means clustering (i.e. unsupervised machine learning feature extraction) and create a plot where the points are colored by cluster membership.

color_scatter <- ggplot(dataset, aes(x=systemsizeac,y=cost, shape=self, color=self)) + geom_point() + labs (title="Scatterplot to inspect for k-means",x= "\nSystem size [kW]",y="System cost [$]")+ theme_bw()
color_scatter

kmeans <-kmeans(dataset2,centers=2)
dataset2$cluster <- as.character(kmeans$cluster)

kmeans_scatter <- ggplot(dataset2, aes(x=systemsizeac,y=cost, shape=cluster, color=cluster)) + geom_point() + labs (title="K-means clustering",x= "\nSystem size [kW]",y="System cost [$]")+ theme_bw()
kmeans_scatter

# Problem 14 ----

# Large complex data sets are becoming more common in this 'age of information'. Simple models of probability density, which often assume that underlying variability is due to a single process (e.g. sampling error) often can fail to describe fully the underlying 'latent' or unobservable aspects of this kind of data. Therefore we need to employ more sophisticated algorithms to fit models to complex data. Multimodal density distributions are often indicative of underlying latent variability. One common solution to fitting density functions to multimodal density is to use a Gaussian Mixture Model (GMM) fitted by the Expectation-Maximization (EM) algorithm as we discussed in class. 
# 1 - In this last part of the assignment you will fit three simple probability density functions (normal, lognormal and another of your own choosing) to the distributions of one of the independent variables (or PC1 representing the most reduced form of your data if your prefer). 
# 2 - You will also fit a GMM as a fourth model of probability density.
# 3 - All models will be compared and evaluated using the Bayesian Information Criterion (BIC) as a method of multimodel inference. 
# 4 - Your goal is to determine which model of the four (normal, lognormal, yourchoice and GMM) best fits the distribution (Iris flower size will be demo in class)

# Problem 15 ----

# Write an R script that:
# 1 - Imports or loads the data set,
# 2 - Fits the model to your distribution,
# 3 - Calculates likelihoods, 
# 4 - Compares BIC, and 
# 5 - Plots histogram(s) of size with each model fitted to the data.
# 6 - Did the model testing indicate the presence of latency (via GMM having best fit)?
# 7 - Be sure to include verbal interpretation printed to screen and/or output files.

ind_variable <- unitcost # define independent variable to be tested

normal <- fitdistr(ind_variable, densfun='normal')
lognormal <- fitdistr(ind_variable, densfun='log-normal')
exponential <- fitdistr(ind_variable, densfun='exponential')
GMM <- normalmixEM(ind_variable)
GMM_loglikelihood <- GMM$loglik

BIC_GMM <- -2*GMM_loglikelihood+4*log(150)
BIC_fits <- BIC(normal, lognormal, exponential)

hist_normal <- ggplot(dataset, aes(x=ind_variable)) + geom_histogram(aes(y=..density..)) + geom_density() + stat_function(fun=dnorm, color="red", args=list(mean = normal$estimate[1], sd = normal$estimate[2])) + labs(title='Normal fit',x='Unit cost [$/kW]', y='Density')+ theme_bw()

hist_lognormal <- ggplot(dataset, aes(x=ind_variable)) + geom_histogram(aes(y=..density..)) + geom_density() + stat_function(fun=dlnorm, color="red", args=list(meanlog = lognormal$estimate[1], sdlog = lognormal$estimate[2])) + labs(title='Log-normal fit',x='Unit cost [$/kW]', y='Density')+ theme_bw()

hist_exponential <- ggplot(dataset, aes(x=ind_variable)) + geom_histogram(aes(y=..density..)) + geom_density() + stat_function(fun=dexp, color="red", args=list(rate = exponential$estimate[1]))+ labs(title='Exponential fit',x='Unit cost [$/kW]', y='Density')+ theme_bw()

hist_GMM <- ggplot(dataset, aes(x=ind_variable)) + geom_histogram(aes(y=2*(..density..))) + geom_density(aes(y=2*(..density..))) + stat_function(fun=dnorm, color="red", args=list(mean = GMM$mu[1], sd = GMM$sigma[1])) + stat_function(fun=dnorm, color="red", args=list(mean = GMM$mu[2], sd = GMM$sigma[2]))+ labs(title='GMM fit',x='Unit cost [$/kW]', y='2*Density')+ theme_bw()

# RUN OUTPUTS SECTION AGAIN (AFTER SOURCING) ----
    # Text outputs ----
sink(file = 'ElenesAGN_A1_TextOutputs&Interpretations.txt')
# Problem 3 ----
writeLines(" \nProblem 3")
writeLines(" \n\n Mean and standard deviation \n")
print(stats_random)
# Problem 4 ----
writeLines(" \n\nProblem 4")
writeLines(" \nOne-way ANOVA test")
print(age.anova)
writeLines("The one-way analysis of variance (ANOVA) compares the means between three or more independent (unrelated) groups. Here, the groups are defined by the technology categories (biomass plants make up group 1, coal plants group 2, etc...) and the variable whose mean we are comparing is the power plant age or lifespan.   \n \nThe null hypothesis assumes the mean power plant lifespan between groups is statistically the same i.e. that any differences between them are random, while the alternative hypothesis states that the differences between at least 2 of the groups compared are not random i.e. power plants lifespan is consistently different at least between 2 of the groups.   \n\nSince the p-value is less than 0.05, we must reject the null hypothesis. This indicates that power plants lifespan is consistently different at least between two of the technologies analyzed. \n")
writeLines(" \n\nBonferroni pairwise t-test \n")
writeLines("In the previous one-way ANOVA test we rejected the null hypothesis, implying that there are significant differences between at least two technologies. But this does not tell us which technologies are significantly different; for this purpose, we can use the Bonferroni method.   \n\nBonferroni's method helps to answer 'which technologies are significantly different from each other?' by performing a pairwise comparison of the means i.e. by comparing all the means against each other e.g. given 3 technologies: coal, biomass & gas, compare all pairs of them: coal vs biomass, coal vs gas, and biomass vs gas. For our case, the number of pairs is given by k = (a)*(a-1)*(1/2) where a = number of groups. Since a = 10, then here we will compare k = (10)*(10-1)*(1/2) = (5)*(9) = 45 pairs, as shown below.")
print(age.pairwise.bonf)
writeLines("\nMany of the comparisons yielded significant differences. Nuclear plants have the least significant differences in average lifespan with respect to the rest of technologies and it only differs significantly from that of hydroelectric plants, which on the other hand are significantly different from all the other groups, followed by coal & gas with 8, oil & other with 5, and biomass, wind, solar & geothermal with 4.")
writeLines(" \n\n\nBenjamini-Hochberg pairwise t-test \n")
writeLines("The multiple test correction in the previous Bonferroni pairwise t-test sets alpha equal to 0.05/k where k is the number of tests run, 45. This correction is more conservative and quickly reduces statistical power with larger sets of tests. Another way to perform correction when doing multiple testing is to determine our own rate of false discovery (FDR) with the Benjamini-Hochberg (BH) correction. BH correction controls the expected proportion of false positives by modifying p-values based on the rank of the data.")
print(age.pairwise.benjhoch)
writeLines("\nWe obtain 32 significant pairwise differences with BH correction, 5 more than those found with the Bonferroni method. Nuclear plants are now different from solar, coal, and wind ones, and solar plants from biomass & geothermal. 13 of the 45 pairwise comparison remain, which may indicate similar lifespans; most of these belong to comparisons between renewable/clean technologies.")
# Problem 5 ----
writeLines(" \n\n\nProblem 5")
writeLines(" \n\nDid this assumption (of normality) hold? \n")
print(age.norm)
print(size.norm)
writeLines(" \nThe KS tests for both variables indicate that neither the power plants physical life nor their capacities follow a normal distribution. \n")
writeLines(" \nDid outcomes of parametric and nonparametric tests appear consistent? Why or why not? \n")
writeLines(" \nKruskal Wallis test")
print(age.kruskal)
writeLines(" \nOne-way ANOVA test")
print(age.anova)
writeLines("\nYes, they are consistent. Both the parametric one-way ANOVA test and the non-parametric Kruskal Wallis test indicate that power plant lifespan is significantly different between some of the technology categories compared. \n")
# Problem 6 ----
writeLines(" \n\nProblem 6")
writeLines(" \n\nHow does result (from regression) compare (to correlation)? \n")
writeLines(" \nCorrelation tests")
print(age.spearman)
print(age.pearson)
writeLines(" \nLinear regression")
print(OLStest)
writeLines(" \nThe rejection of the null hypothesis in the correlation tests indicates that age and capacity are positively correlated. The regression shows a significant but small parameter for the hypothesis that capacity causes age.  \n")
writeLines(" \n\nWhen is regression appropriate instead of correlation? \n")
writeLines(" \nRegression is appropriate when testing a relation of causation. Correlation only tests if there is some relationship or connection, without assuming causation in any direction.  \n")
# Problem 8 ----
writeLines(" \n\nProblem 8")
writeLines(" \n\nMANOVA test \n")
print(manova)
writeLines(" \n\nThe results from MANOVA indicate the differences between the groups are significant. \n")
# Problem 9 ----
writeLines(" \n\nProblem 9")
writeLines(" \n\nWhich one is the best predictor (from the multiple linear regression)? \n")
writeLines(" \nMultiple linear regression test")
print(stats_mlr)
writeLines(" \n\nThe best predictor for system cost is EVcount, the number of electric vehicles in the facility, with the largest coefficient, followed by system size. Tax credit causes a decrease in system cost. All the explanatory variables had significant parameters.  \n")
writeLines(" \nMultiple linear regression test for a subgroup of self-installers")
print(stats_mlr2)
writeLines(" \n\nFor this subset the best predictor is system size. EV count is not signficant anymore and the parameter for tax credit is smaller.  \n")
# Problem 10 ----
writeLines(" \n\nProblem 10")
writeLines(" \nANCOVA")
print(stats_ancova)
writeLines(" \n\nSystem size remains a significant predictor of system cost even when controlling for the unit cost.  \n")
# Problem 11 ----
writeLines(" \n\nProblem 11")
writeLines(" \nHow successful was the reduction of the data?\n")
print(stats_pca)
writeLines(" \n\nBy looking at the proportion of variance explained by each component above (also shown in the Scree Plot, see pdf document with graphical outputs) around 40% of the variation is explained by PC 1, while PC2 and PC3 explain about 25% each, and PC4 only 8%. With only 3 components we can explain 92% of the variation in the data.  \n")
writeLines(" \nHow would you interpret the loadings of the original measurements on each of the PC's?\n")
print(pca_rotations)
writeLines(" \n\nPV system cost and size influence PC1 the most. PC1 could relate to the overall size of the system, as a bigger system means greater costs and a larger PV system size. Tax credit is weakly negatively correlated and EVcount has little influecence.
            \nPC2 is strongly negatively correlated with tax credit, and moderately correlated with EV count; cost and PV system size have a negative sign but low influence. This means that PC2 decreases with larger tax credits, and increases with more electric vehicles in the facilities. 
            
            \nPC3 also has strong parameters for tax credit and EV count, but they share the sign: PC3 decreases as tax credits and electric vehicles increase.
           
           \nIn PC4 cost and size have strong loadings but opposite signs. The principal components 2 through 4 may relate to characteristics of specific subgroups of population. \n")
writeLines(" \nWhich one best describes overall size variation?")
writeLines(" \n\nNone of the principal components presented all-positive loadings, but PC1 probably reflects better the overall size variation.  \n")
writeLines(" \nDo any of the PC's appear strongly driven by one particular measurement (see loadings on PC)?")
writeLines(" \n\nYes, PC2 is strongly driven by taxcredit, and PC3 by EVcount.  \n")
# Problem 12 ----
writeLines(" \n\nProblem 12")

writeLines(" \nFactor Analysis")
print(factanalysis)

writeLines(" \n\nI was able to carry out the factor analysis with one factor only, and the null hypothesis that this single factor is sufficient must be rejected. The one factor only explains 35% of the variation. We obtained very high values for uniqueness in EVcount and tax credit, and relatively low for sistem size and cost. The loadings were high for system size and cost, and negligible for EVcount and tax credit.  \n")

writeLines(" \nCan you determine any common latent underlying 'factors' in your data?")
writeLines(" \n\nThe result is inconclusive: only 1 factor was obtained that does not explain much and was not significant.  \n")

writeLines(" \nHow many of such factors are significant?")
writeLines(" \n\nNone of them.  \n")

writeLines(" \nWhat important characteristics can you find that tend to cluster together?")
writeLines(" \n\nSystem size dominated the single factor, followed by cost, but the analysis did not show any particular trend.  \n")

writeLines(" \nIf two traits are far apart along the axis of a significant factor, what does this indicate?")
writeLines(" \n\nThat they are inversely correlated. \n")

writeLines(" \nWhat if they are close together?")
writeLines(" \n\nThat could mean they are positively correlated. \n")

writeLines(" \nWas this factor analysis particularly successful in identifying a 1 or 2 large underlying latent variables that define something about your data?")
writeLines(" \n\nNo, because results were inconclusive. \n")
# Problem 13 ----
writeLines(" \n\nProblem 13")
writeLines(" \nVisually determine: how many potential 'clusters' of data points you see in your plot?")
writeLines(" \n\nThere aren't obvious clusters. Because of the coloring, I could say there are two clusters. \n")
# Problem 15 ----
writeLines(" \n\nProblem 15")
writeLines(" \nDid the model testing indicate the presence of latency (via GMM having best fit)?")
writeLines(" \nBIC for GMM")
print(BIC_GMM)
writeLines(" \nBIC for the rest of the models tested")
print(BIC_fits)
writeLines(" \n\nYes, the model with the best fit was the Gaussian Mixture Model (GMM) that showed the least value for the Bayesian Information Criterion (BIC). \n")
# end of text outputs ----
sink()

    # Graphical outputs ----
pdf(file = 'ElenesAGN_A1_GraphicOutputs.pdf')
# Problem 3 ----
multiplot(hist_random)
# Problem 4 ----
multiplot(ErrorBar_AgeVSFuelCat)
# Problem 5 ----
multiplot(Pearson_scatter)
# Problem 6 ----
multiplot(OLSplot)
# Problem 11 ----
multiplot(pca_screeplot,pca_1vs2, pca_2vs3,pca_3vs4,cols=2)
# Problem 13 ----
multiplot(color_scatter,kmeans_scatter)
# Problem 15 ----
multiplot(hist_normal, hist_lognormal,hist_exponential, hist_GMM, cols=2)
# end of graphical outputs ----
dev.off()