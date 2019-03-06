# Title     : Repeated Measures ANOVA
# Objective : Model selection for MMN-Dataset
# Created by: MEkramnia
# Created on: 11/2/19

# Install packages if not already installed
if(!require(psych)){install.packages("psych")}
if(!require(nlme)){install.packages("nlme")}
if(!require(car)){install.packages("car")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(rcompanion)){install.packages("rcompanion")}

install.packages("afex")
install.packages("emmeans")

# Imports
library(psych)
library(nlme)
library(car)
library(rcompanion)
library(multcompView)
library(emmeans) # emmeans must now be loaded explicitly for follow-up tests.
library(rcompanion)
library(ggplot2)
library(tidyr)
library(afex) # needed for ANOVA functions.aov_ez
library(dplyr)
library(broom)
library(multcompView)

MMNdf <- read.csv(file = "/Users/minaekramnia/Github/badbaby/badbaby/data/Ds-mmn-2mos_meg_df.tsv",
                  header=TRUE, sep="\t")

#Model 2: Using aov_ez() function to fit a mixed ANOVA with treatment as between-subjects factor and phase as within-subject factors.
#It seems that aov_ez() codes the IVs with effects coding by default, rather than dummy coding as in base R's lm aov; 
#The message Contrasts set to contr.sum for the following variables: factor_a, factor_b is a warning of this. 
#When we set check.contrasts to FALSE, aov_ez uses dummy coding and we get concordant Fs and Ps 
#using afex package: 
#To use this data in an ANOVA one needs to aggregate the data such that only one observation per participant 
#and cell of the design (i.e., combination of all factors) remains. 
#afex does this automatically for us. 

fit_nice <- aov_ez("subjId", "latencies", MMNdf, between = c("sesGroup", "birthWeight", "nSibs", "maternalEdu", "paternalEdu", "age"), within = c("hemisphere", "oddballCond"), return = "nice" )
fit_nice

####Behavioural Data
aov_cdi <- aov_ez("subject", "m3l", MMNcdi, between = c("sesGroup", "age", "maternalEdu", "paternalEdu", "nSibs"),  
                  return = "nice" )
aov_cdi
######
#Results summary:
#This function is equivalent to model 1: aov_car and produce exactly the same results in F and Pvalue: Not significant. 
