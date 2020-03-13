# Script for LMMs and other auxillary stats for the MW_HDDM study

# Start by making quantile plots for state overall. 

# Load required packages

library(lme4) # freq LMMs
library(brms) # bayesian LMMs
library(bayesplot) # plotting MCMC posteriors
library(sjstats) # std_beta for standardised betas
library(lmerTest) # extends lme4 to give p values
library(ggplot2) # plotting
library(plotrix)
library(scales)
library(dplyr)

d1 <- read.csv('./Data/Gotrials_rtbin5.csv') 

subject_means <- group_by(d1,subj_idx,State,rt_bin5) %>% summarize(accuracy = mean(Corr), rt=mean(rt), Vig=mean(Vig),na.rm=T)
subject_means

summary <- subject_means %>% group_by(State,rt_bin5) %>% 
  summarise_each(funs(mean,sd,std.error),na.rm=T)
summary

# RT
p<- ggplot(summary, aes(x=as.factor(rt_bin5), y=rt_mean, group=State, color=State)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=rt_mean-rt_std.error, ymax=rt_mean+rt_std.error), width=.2,
                position=position_dodge(0.05))
print(p)

# Accuracy
p<- ggplot(summary, aes(x=rt_bin5, y=accuracy_mean, group=State, color=State)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=accuracy_mean-accuracy_std.error, ymax=accuracy_mean+accuracy_std.error), width=.2,
                position=position_dodge(0.05))
print(p)
