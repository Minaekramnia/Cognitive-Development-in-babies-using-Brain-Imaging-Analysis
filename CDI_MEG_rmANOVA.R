# Title     : Repeated Measures ANOVA
# Objective : Model selection for MMN-Dataset
# Created by: Mekramnia
# Created on: 11/5/18

###########################
# Repeated Measures ANOVA #
###########################
# When an experimental design takes measurements on the same experimental unit # over time, the
# analysis of the data must take into account the probability that measurements for a given
# experimental unit will be correlated in some way. One approach to deal with non-independent
# observations is to include an autocorrelation structure in the model using
# the nmle package, e.g. treating Subject_ID as a random variable to take into account
# native differences among subjects
# https://rcompanion.org/handbook/I_09.html

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
library(emmeans)
library(rcompanion)
library(ggplot2)
library(tidyr)
library(afex)
library(dplyr)
library(broom)
library(multcompView)

#Reading the data 
#MEG
MMNdf <- read.csv(file = "/Users/minaekramnia/Github/badbaby/badbaby/data/Ds-mmn-2mos_meg_df.tsv",
header=TRUE, sep="\t")

#Behavioral
MMNcdi <- read.csv(file = "/Users/minaekramnia/Github/badbaby/badbaby/data/Ds-mmn-2mos_cdi_df.tsv",
                  header=TRUE, sep="\t")
MMNdf.mean <- aggregate(MMNdf$auc,
by=list(MMNdf$subjId, MMNdf$sesGroup,
MMNdf$oddballCond, MMNdf$hemisphere),
FUN='mean')

# Order factors by the order in data frame, otherwise, R will alphabetize them
#MMNdf$Ses = factor(MMNdf$SES,levels = unique(MMNdf$SES))

# Check the data frame
headTail(MMNdf)
str(MMNdf)
summary(MMNdf)

##############################
# Mina WIP Analysis of Variance
##############################
# ANOVA Mixed:
# IV between: sesGroup, birthWeight, nSibs, maternalEdu, paternalEdu
# IV within (Repeated Measures): oddballCond, hemisphere
# DV MEG:         auc, latencies
# DV Behavioral:  vocab, m3l

#plot summary of the data: 
knitr::kable(describe(MMNdf), digits = 2, caption = 'MEG Data Summary')
knitr::kable(describe(MMNcdi), digits = 2, caption = 'Behavioural Data Summary')
Latencies <- MMNdf$latencies
hist(Latencies , n=20)
Number_of_Siblings<-MMNdf$nSibs
hist(Number_of_Siblings, n=10)

Paternal_Education<- MMNdf$paternalEdu
Maternal_Education<- MMNdf$maternalEdu
hist(Paternal_Education)
hist(Maternal_Education)

#Model1: In this model, we use function aov_car() from the car package
#The formula passed to aov_car basically needs to be the same as for standard aov with a few differences:
#It must have an Error term specifying the column containing the participant (or unit of observation) identifier (e.g., minimally +Error(id)). 
#This is necessary to allow the automatic aggregation even in designs without repeated-measures factor.
#Since 'auc' is very small, we get the log. 'ges' provides generalized eta-squared, ‘the recommended effect size statistics for repeated measure designs’

#The two responses per cell of the design and participants are aggregated for the analysis as indicated by the warning message. 
# Sphericity correction method: GG 

aov.bww <- aov(latencies ~ oddballCond * hemisphere * sesGroup 
                   + Error(subjId/(oddballCond*hemisphere)), 
                   data=MMNdf ) #default sphericity corrections disabled
summary(aov.bww)
#MMNdf$log_auc <- log(MMNdf$auc) 
aov.bww_auc <- aov(auc ~ oddballCond * hemisphere * sesGroup  
                   + Error(subjId/(oddballCond*hemisphere)), 
                  data=MMNdf,anova_table = list(correction = 'none'), 
                       fun_aggregate = mean) #default sphericity corrections disabled
summary(aov.bww_auc)

#Alternatively, using the nlme package for a lme style ANOVA:
lme_latencies<-lme(latencies ~  oddballCond * hemisphere * sesGroup , random = ~1|subjId/oddballCond, data=MMNdf)
summary(glht(lme_latencies, linfct = mcp(sesGroup = "Tukey")))

#Pairwise comparisons
#To interpret the interaction we use pairwise comparison. 
#we start by getting the marginal means, and use emmeans 
#To test the significance of the mean differences, we will use the pairs function. 
#Correct the significance values, using holm adjustment, because it is uniformly more powerful than Bonferroni.

within_fitted_interaction <- emmeans(aov.bww, ~sesGroup|hemisphere|oddballCond)
within_fitted_interaction2 <- emmeans(aov.bww_auc, ~sesGroup|hemisphere|oddballCond)
ref <-pairs(within_fitted_interaction2, adjust = "holm") #pairwise comparison with no correction

#########
#Results: This model shows no significant effect for the IVs on auc and latencies. 

########Plots : 
Sum_lat = groupwiseMean(latencies ~ sesGroup + hemisphere,
                        data=MMNdf, conf=0.95, digits=3, traditional=FALSE, percentile=TRUE)
Sum_lat

Sum_auc = groupwiseMean(auc ~ sesGroup + hemisphere + oddballCond,
                        data=MMNdf, conf=0.95, digits=3, traditional=FALSE, percentile=TRUE)
Sum_auc
pd = position_dodge(.2)
ggplot(Sum_auc, aes(x=hemisphere, y=within_fitted_interaction, color=sesGroup)) +
  geom_errorbar(aes(ymin=Percentile.lower,
                    ymax=Percentile.upper),
                width=.2, size=0.7, position=pd) +
  theme_bw() +
  theme(axis.title=element_text(face="bold")) +
  ylab("Peak Magnitude")


#############
# Behavioural Data Analysis: 
# IV betwen = sesGroup
# IV within = maternalEdu, paternalEdu, nSibs, birthWeight
# outcome variable = m3l, vocab
aov.bww_cdi <- aov(m3l ~ maternalEdu * paternalEdu * sesGroup *nSibs * birthWeight , 
                   data=MMNcdi, anova_table = list(correction = 'none'), 
                   return = "nice" , fun_aggregate = mean) #default sphericity corrections disabled
summary(aov.bww_cdi)

posthoc <- emmeans(aov.bww_cdi, ~maternalEdu)

aov.bww_vocab <-aov(vocab ~ maternalEdu * paternalEdu * sesGroup * nSibs , data=MMNcdi)
summary(aov.bww_cdi)

fiber.lm <- lm(m3l ~  maternalEdu * birthWeight, data=MMNcdi)
posthoc_m3l <- emtrends(fiber.lm, pairwise ~ maternalEdu, var = "birthWeight")

#ggplot(Sum_m3l, aes(x=age, y=posthoc_m3l, color=sesGroup)) +
 # geom_errorbar(aes(ymin=Percentile.lower,
  #                  ymax=Percentile.upper),
  #              width=.2, size=0.7, position=pd) +
  #theme_bw() +
  #theme(axis.title=element_text(face="bold")) +
  #ylab("Mean Lenght of Utterance")

###M3l-age/ses
aov.bww_m3l <- aov(m3l ~  sesGroup * birthWeight * nSibs * maternalEdu * paternalEdu * age,
                   data=MMNcdi,anova_table = list(correction = 'none'), 
                   fun_aggregate = mean)
summary(aov.bww_m3l)

#within_fitted_interaction3 <- TukeyHSD(aov.bww_m3l, "nSibs")

Sum_m3l = groupwiseMean(m3l ~ sesGroup + age,
                        data=MMNcdi, conf=0.95, digits=3, traditional=FALSE, percentile=TRUE)
Sum_m3l
ggplot(Sum_m3l, aes(x=age, y=within_fitted_interaction3, color=sesGroup)) +
  geom_errorbar(aes(ymin=Percentile.lower,
                    ymax=Percentile.upper),
                width=.2, size=0.7, position=pd) +
  theme_bw() +
  theme(axis.title=element_text(face="bold")) +
  ylab("Mean Lenght of Utterance")

################################
#OLS for MEG and Behavioral data
################################
ols_meg <-glm(latencies ~ oddballCond + hemisphere + sesGroup + oddballCond*hemisphere + hemisphere*sesGroup + 
           sesGroup*oddballCond + oddballCond*hemisphere*sesGroup, data = MMNdf)
ols_meg
summary(ols_meg) # sesGrouplow is signifacantly correlated with auc

ols_meg_lat <-glm(latencies ~ oddballCond * hemisphere * sesGroup * age * birthWeight , data = MMNdf)
summary(ols_meg_lat)

ols_meg_auc <-glm(auc ~ oddballCond * hemisphere * sesGroup * age * birthWeight , data = MMNdf)
summary(ols_meg_auc)

ols_cdi <-glm(vocab ~ age + sesGroup + maternalEdu + paternalEdu + nSibs + birthWeight + age:sesGroup, data = MMNcdi)
summary(ols_cdi)

ols_cdi_m3l <-glm(m3l ~ age + sesGroup + maternalEdu + paternalEdu + nSibs + birthWeight + age:sesGroup, data = MMNcdi)
summary(ols_cdi_m3l)

###Result: birthweight has significant effect on the number of vocabs. Also the paternal education C and the 
#interaction of age and sesGroup low

#M3l is significantly correlated with age, ses group low, paternal educ f and g and k, number of siblings,
#birthweight, and the interaction of age and ses group low

#Vocab is significantly corelated with paterna education c and birthweight and the interaction of age and ses group

#latencies is significantly correlated with ses group low, interaction of age and ses group low, 
#interaction of ses and birthweight and the interaction of age, bw and ses

#auc is not significantly affected by any correleted variable

######################
# Mixed effect model #
######################
# Fit a linear mixed-effects model in the formulation described in Laird and
# Ware (1982) but allowing for nested random effects. The within-group errors
# are allowed to be correlated and/or have unequal variances. Model is fit by
# maximizing the restricted log-likelihood (REML).
# Covariance structure corAR1 with form ~1|subjvar is used to indicate a
# temporal autocorrelation structure of order one, where time covariate is
# order of observations, and Subject_ID indicates the variable for
# experimental units.  Autocorrelation is modeled within levels of the
# subjvar, and not between them.
# NB Autocorrelation structures can be chosen by either of two methods. The
# first is to choose a structure based on theoretical expectations of how the
# data should be correlated. The second is to try various autocorrelation
# structures and compare the resulting models with a criterion like AIC,
# AICc, or BIC to choose the structure that best models the data.

fmla <- as.formula("latencies ~ ses_label + condition_label + hem_label +
                    ses_label * condition_label + ses_label * hem_label +
                    condition_label * hem_label + ses_label *
                    condition_label * hem_label")

# Mixed effect model
model.mixed = lme(fmla, random=~ 1 | Subject_ID,
correlation=corAR1(form=~ 1 | Subject_ID), data=MMNdf, method="REML")

summary(model.mixed)
Anova(model.mixed)

######################
# Fixed effect model #
######################
model.fixed = gls(fmla, correlation=corAR1(form=~1|Subject_ID),
data=MMNdf, method='REML')

summary(model.fixed)
Anova(model.fixed)

##############################
# Model comparison/selection #
##############################
# The random effects in the model can be tested by comparing the model to a
# model fitted with just the fixed effects and excluding the random effects.
# Because there are not random effects in this second model, the gls function
# in the nlme package is used to fit this model.

anova(model.mixed, model.fixed) # compares the reduction in the residual sum of squares
# NB the residual sum of squares (RSS), also known as the sum of squared
# residuals (SSR) or the sum of squared errors of prediction (SSE), is the sum
# of the squares of residuals (deviations predicted from actual empirical
# values of data). It is a measure of the discrepancy between the data and an
# estimation model. A small RSS indicates a tight fit of the model to the data.
# It is used as an optimality criterion in parameter selection and model
# selection.

###########################################################
# Pseudo R^2 Measures of Fit for fixed model of latencies #
###########################################################
# The nagelkerke function can be used to calculate a p-value and pseudo
# R-squared value for the model. One approach is to define the null model as
# one with no fixed effects except for an intercept, indicated with a 1 on
# the right side of the ~.  And to also include the random effects, in this
# case 1|Subject_ID.
# A signficant pseudo R-squared indicates this model better fits the outcome
# data than the Null model. While pseudo R-squareds cannot be interpreted
# independently or compared across datasets, they are valid and useful in
# evaluating multiple models predicting the same outcome on the same dataset.
# In other words, a pseudo R-squared statistic without context has little
# meaning. A pseudo R-squared only has meaning when compared to another
# pseudo R-squared of the same type, on the same data, predicting the same
# outcome.  In this situation, the higher pseudo R-squared indicates which
# model better predicts the outcome.
# https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/

s

#####################
# Post-Hoc analysis #
#####################
marginal = emmeans(model.fixed, ~ ses_label : hem_label)
cld(marginal, alpha=0.05, Letters=letters, ### Use lower-case letters for .group
    adjust="tukey")     ###  Tukey-adjusted comparisons

Sum = groupwiseMean(latencies ~ ses_label + hem_label,
    data=MMNdf, conf=0.95, digits=3, traditional=FALSE, percentile=TRUE)
Sum

pd = position_dodge(.2)
ggplot(Sum_lat, aes(x=hem_label, y=emmeans, color=ses_label)) +
    geom_errorbar(aes(ymin=Percentile.lower,
        ymax=Percentile.upper),
        width=.2, size=0.7, position=pd) +
    geom_point(shape=15, size=4, position=pd) +
    theme_bw() +
    theme(axis.title=element_text(face="bold")) +
    ylab("Peak Latency")
