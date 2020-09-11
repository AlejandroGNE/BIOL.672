#### Alejandro G.N. Elenes ####
# Code written in macOS Catalina version 10.15.6
# R packages needed to run this script:
install.packages("psych")
install.packages("readxl")
install.packages("openxlsx")
install.packages("reshape2")
install.packages("deplyr")
library(psych)
library(readxl)
library(openxlsx)
library(reshape2)
library(dplyr)
# Data files necessary to run this script:

#### 2 ####
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

#### 3 ####
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

#### 4 ####
# Find an appropriate data set for a one-way ANOVA. Shape 
# the data so you can use the read.table function from an 
# input data file (as demo during class). Extend your R 
# script so that that it conducts a oneway ANOVA (i.e. 
# oneway.test(values, category)). Add to this script, an 
# error bar chart output of your categories with colored 
# bars of your choice. Lastly add pairwise t tests for each
# category pair and use both Bonferroni and 
# Benjamini-Hochberg multiple test correction using the 
# p.adjust function in R. As in #3, export output to 
# appropriately named text and pdf files.Be sure to include 
# verbal interpretation printed to screen and/or output files.

#### 4 - about the dataset ####

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

#### 4 - code ####

Retirement_dataset <- read.xlsx("3_1_Generator_Y2018.xlsx", sheet = "Retired and Canceled", startRow = 2, na.strings = " ")
Retirement_dataset <- Retirement_dataset[complete.cases(Retirement_dataset$Retirement.Year),] # drop missing values (cancelled power plants)
Online <- Retirement_dataset$Operating.Year
Retired <- Retirement_dataset$Retirement.Year
Source <- Retirement_dataset$Energy.Source.1
Capacity <- Retirement_dataset$`Nameplate.Capacity.(MW)`
Fuel_unique_categories <- unique(Source)
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

