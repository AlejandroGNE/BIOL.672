# Alejandro G.N. Elenes
# Problem 1 ----
# Code written in macOS Catalina version 10.15.6
# R packages needed to run this script:
library(ggplot2)
library(psych)
library(readxl)
library(openxlsx)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(Rmisc)
# Data files necessary to run this script:
# 3_1_Generator_Y2018.xlsx
# Problem 2 ----
# Create a dataset of 5000 random numbers generated in R 
# using any model distribution you like. Write an R script 
# that reports the sample mean and sample standard deviation 
# and then plots a histogram from this data and that includes 
# both a density line laid over the histogram bars and a 
# normal curve fitted to the data and overlaid as well.

random_data <- rnorm(5000) # create variable with 5000 random numbers using normal distribution
describe(random_data) # use describe function from psych package to report mean and st dev
hist(random_data, prob = TRUE) # plot a histogram of random data
lines(density(random_data,na.rm=T),col="red",lwd=2) # over histogram, add density line
curve(dnorm(x, mean=mean(random_data), sd=sd(random_data)), add=TRUE) # over histogram, plot normal distribution based on random data

# Problem 3 ----
# To automate the printing of stat test output in R you will 
# want to explore the ‘sink’ and ‘unsink’ functions in R. 
# Use your script to create a file called ‘desc.txt’ and 
# print your mean and SD from #2 to it. Add lines to your 
# script to rename and save Rplots.pdf as histo.pdf.

sink(file = 'desc.txt') # trigger printing to text file
describe(random_data) # print stats of random data
sink() # stop printing to text file
pdf(file = 'histo.pdf', width = 4, height = 4) # trigger printing to pdf file
hist(random_data, prob = TRUE) # plot histogram again, this time prints into pdf
lines(density(random_data,na.rm=T),col="red",lwd=2) 
curve(dnorm(x, mean=mean(random_data), sd=sd(random_data)), add=TRUE)
dev.off() # stop printing to pdf

# Problem 4 ----

# Find an appropriate data set for a one-way ANOVA. Shape 
# the data so you can use the read.table function from an 
# input data file (as demo during class). 

# Extend your R script so that that it conducts a oneway 
# ANOVA (i.e. oneway.test(values, category)). 

# Add to this script, an error bar chart output of your 
# categories with colored bars of your choice. 

# Lastly add pairwise t tests for each category pair and 
# use both Bonferroni and Benjamini-Hochberg multiple test 
# correction using the p.adjust function in R. 

# As in #3, export output to appropriately named text and pdf 
# files.

# Be sure to include verbal interpretation printed to screen 
# and/or output files.

  # 4.1 - about the dataset ----

# I will look for trends in power plants retirements in the US. e.g:
# how old are power plants at age of retirement?
# how much do plants physical life differs from their originally planned lifetime? 
# Is there a trend, and is this trend different among different technologies?

# I will use the dataset from survey form EIA-860, which
# collects generator-level specific information about 
# existing and planned generators and associated
# environmental equipment at electric power plants with
# 1 MW or greater of combined nameplate capacity.
# Raw data comes in a zip folder with multiple excel 
# files. I use the sheet 'Retired and Canceled' from 
# the file '3_1_Generator_Y2018.xlsx'

  # 4.2 - code ----
  # 4.2.1 - load data, create variables ----
Retirement_dataset <- read.xlsx("3_1_Generator_Y2018.xlsx", sheet = "Retired and Canceled", startRow = 2, na.strings = " ")
Retirement_dataset <- Retirement_dataset[complete.cases(Retirement_dataset$Retirement.Year),] # drop missing values (cancelled power plants)
Online <- Retirement_dataset$Operating.Year # Year plant went online / was born
Retired <- Retirement_dataset$Retirement.Year # Year plant went offline / was retired / died
Source <- Retirement_dataset$Energy.Source.1 # Primary fuel
Capacity <- Retirement_dataset$`Nameplate.Capacity.(MW)` # Plant capacity in MW
Age <- Retired - Online # Plant phyisical life / age
Fuel_unique_categories <- unique(Source)
  # 4.2.2 - Group into higher level categories based on fuel type for easier tractability----
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
  # 4.2.3 - One way ANOVA ----
age.anova <- oneway.test(Age~Fuel_categories)
age.anova
  # 4.2.4 - Error bar chart ----

# Try using geom_errorbar
std.err.plot.agevsfuelcategories <- ggplot(Retirement_dataset, aes(x = Fuel_categories, y = Age)) +
  geom_errorbar(aes(
      ymin=Age-sd(Age)/sqrt(length(Age)),
      ymax=Age+sd(Age)/sqrt(length(Age))), 
    width=.2, 
    position=position_dodge(.9))
std.err.plot.agevsfuelcategories # print plot
# didn't work, geom_errorbar creates subgroups within my categories 
# and plots error bars for each subgroup

# Try using Rmisc
Retirement_dataset$Fuel.Categories <- Fuel_categories
Retirement_dataset$Age <- Age
stats <- summarySE(Retirement_dataset, measurevar = "Age", groupvars = "Fuel_categories")
stats
ErrorBar_AgeVSFuelCat <- ggplot(stats, aes(x=Fuel_categories, y=Age)) + 
     geom_bar(stat="identity", fill="black", width=0.8) +
     geom_errorbar(aes(ymin=Age-se, ymax=Age+se), width=0.2, colour='red') +
    labs(x = "\n Technology", y = "Physical life [years]\n") +
  theme(axis.text.x=element_text(angle=15))
ErrorBar_AgeVSFuelCat # print plot

  # 4.2.5 - Pairwise t-tests ----
# Bonferroni
age.pairwise.bonf <- pairwise.t.test(Age, Fuel_categories, p.adj = "bonf")
age.pairwise.bonf
# Benjamini-Hochberg
age.pairwise.benjhoch <- pairwise.t.test(Age, Fuel_categories, p.adj = "BH")
age.pairwise.benjhoch
# Some descriptive statistics
box.plot.agevsfuelcategories <- ggplot(Retirement_dataset, mapping = aes(x = Fuel_categories, y = Age)) + geom_boxplot() + labs(x = 'Fuel categories', y = 'Phyiscal life [years]') + theme(axis.text.x=element_text(angle=15), panel.background = element_rect(fill = 'grey50'))
box.plot.agevsfuelcategories
  # 4.2.6 - Export & interpretation ----
sink(file = 'A1Problem4_Text.txt') # used new text file to avoid rewriting of previous exercises
age.anova
writeLines("The one-way analysis of variance (ANOVA) compares the means between three or more independent (unrelated) groups. Here, the groups are defined by the technology categories (biomass plants make up group 1, coal plants group 2, etc...) and the variable whose mean we are comparing is the power plant age or lifespan.   \n \nThe null hypothesis assumes the mean power plant lifespan between groups is statistically the same i.e. that any differences between them are random, while the alternative hypothesis states that the differences between at least 2 of the groups compared are not random i.e. power plants lifespan is consistently different at least between 2 of the groups.   \n\nSince the p-value is less than 0.05, we must reject the null hypothesis, which indicates that power plants lifespan is consistently different at least between two of the technologies analyzed. \n")
writeLines("In the previous one-way ANOVA test we rejected the null hypothesis, implying that there are significant differences between at least two technologies. But this does not tell us which technologies are significantly different; for this purpose, we can use the Bonferroni method.   \n\nBonferroni's method helps to answer 'which technologies are significantly different from each other?' by performing a pairwise comparison of the means i.e. by comparing all the means against each other e.g. given 3 technologies: coal, biomass & gas, compare all pairs of them: coal vs biomass, coal vs gas, and biomass vs gas. For our case, the number of pairs is given by k = (a)*(a-1)*(1/2) where a = number of groups. Since a = 10, then here we will compare k = (10)*(10-1)*(1/2) = (5)*(9) = 45 pairs, as shown below.   \n\n")
age.pairwise.bonf
writeLines("Many of the comparisons yielded significant differences. Nuclear plants have the least significant differences in average lifespan with respect to the rest of technologies and it only differs significantly from that of hydroelectric plants, which on the other hand are significantly different from all the other groups, followed by coal & gas with 8, oil & other with 5, and biomass, wind, solar & geothermal with 4.")
age.pairwise.benjhoch
writeLines("")
sink()
pdf(file = 'A1Problem4_Plots.pdf') # trigger printing to pdf file
ErrorBar_AgeVSFuelCat
print("")
dev.off() # stop printing to pdf

