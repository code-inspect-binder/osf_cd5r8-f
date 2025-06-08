# Data pre-processing ----

## Clear console and environment
rm(list=ls())

## set working directory
getwd()
#setwd("~/Desktop/Hui&Wu-R/raw_data")
#setwd( "G:/My Drive/HuiWu_Reliability/Data_Packet/raw_data")
#setwd("Volumes/GoogleDrive/My Drive/HuiWu_Reliability/Data_Packet/raw_data")

## load data
d = read.csv("./bind_SPR_English_L2_version1-2.csv", fileEncoding="UTF-8-BOM")

#load packages
library(tidyr)
library(tidyverse)
library(psych)
library(lme4)
library(lmerTest)
library(car)
library(dplyr)
library(rstatix)
library(Rmisc)

## because of conflicts between different packages, we are address them by setting our preferences
library(conflicted)
conflict_prefer("group_by", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")


#Unify the labels for different columns
names(d)[names(d) == "name"] <- "Subject"
names(d)[names(d) == "RT"] <- "Reaction.Time"

##Follow Fang&Wu(2022) data treatment -- excluding participants with comprehension question accuracy < 0.8
# Exp1 should be excluded as CQs are critically manipulated not for attention checking
d = d %>% filter(!str_detect(Condition,"Monkey-Exp1")) 
## delete rows we are not interested in 
d <- droplevels(subset(d, Condition!="consent" & Condition!="background" & Condition!="intro" & Condition!="practice"& Condition!="debrief"))

## filter comprehension question
CQ1 <- d %>% filter(Controller=='QuestionAlt')
names(CQ1)[names(CQ1) == "Reaction.Time"] <- "accuracy"
CQ1$accuracy <- as.numeric(as.character(CQ1$accuracy))

## Compute mean accuracy for each participant for the experimental items 
CQ_acc = summarySE(CQ1, measurevar="accuracy", groupvars=c("Subject"))

## Exclude participants with accuracy less than 0.8
d_lowAcc =  CQ_acc %>% filter(accuracy < 0.8)
d = d[!d$Subject %in% d_lowAcc$Subject,]

rm(list=ls()[! ls() %in% c("d")])

## select and filter useful columns and rows
d = d %>% 
  filter(str_detect(Condition,"Exp3")) %>% 
  filter(Controller != 'QuestionAlt' ) %>% 
  select(7:10, 12,13)

## Changing variables to factors
d$Subject = as.factor(d$Subject)
d$Condition = as.factor(d$Condition)
d$region = as.factor(d$region)

#Create a new column with Simpler Condition Names "Or & EitherOr"
d[substr(d$Condition, 13,21)=="either-no","Cond"]<-"Or"
d[substr(d$Condition,20,21)!="no","Cond"]<-"EitherOr"

#Keep only the critical region (numbered as 3) & Unify the names for labeling
d = d %>% filter(region=='3') %>% select(2,4:7)
names(d)[names(d) == "Cond"] = "Condition"
names(d)[names(d) == "Item"] = "ItemNo"

# fix ItemNo such that each item is a duplet
d <- d %>%
        arrange(Sentence) %>%
        mutate(ItemNo = rep(1:20, each=122))

# only removing those impossible trials (RT < 0)
d_untrimmed<- d %>%
            filter (Reaction.Time >0) %>%
            select (1,2,4,5)

# removing the outliers in a basic way -- RT is 2.5 SD away
# keep only useful columns for analysis
d_trimmed <-  d %>% 
  group_by(Subject) %>% 
  filter(abs(Reaction.Time - mean(Reaction.Time)) < (sd(Reaction.Time) * 2.5)) %>%
  select (1,2,4,5)

rm(d)

## write a function to apply the correction according to Krus and Helmstadter 1993
KH_Correct <- function(r) {
  corrected <- r*-1 / (0.5 * (1-r))
  return(corrected)
}

# Approach 1: Raw RT difference ----

#similar procedure as in the data analysis for Hui et al.
## WITHOUT TRIMMING ----

d <- d_untrimmed

## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
d$ItemNo = parse_number(as.character(d$ItemNo))
de = d %>% filter(ItemNo %% 2 == 0)
do = d %>% filter(ItemNo %% 2 == 1)

#### AVERAGE BY ITEM ----
mean.e <- de %>%
      group_by(ItemNo, Condition) %>%
      summarise(meanRT = mean(Reaction.Time))%>%
      pivot_wider(names_from = Condition, values_from = meanRT) %>%
      mutate(rt.diff = Or - EitherOr)

mean.o <- do %>%
  group_by(ItemNo, Condition) %>%
  summarise(meanRT = mean(Reaction.Time))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

## Table 2 in the main text
### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE)   
KH_Correct(-0.17) #.29

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T) 
KH_Correct(-0.1886) # .32

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT----
mean.e <- de %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(Reaction.Time))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

mean.o <- do %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(Reaction.Time))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

## Table 2 in the main text
### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE) #.16  

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T) #.05

## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

## WITH TRIMMING ----
d <- d_trimmed

## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
d$ItemNo = parse_number(as.character(d$ItemNo))
de = d %>% filter(ItemNo %% 2 == 0)
do = d %>% filter(ItemNo %% 2 == 1)

#### AVERAGE BY ITEM ----
mean.e <- de %>%
  group_by(ItemNo, Condition) %>%
  summarise(meanRT = mean(Reaction.Time))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

mean.o <- do %>%
  group_by(ItemNo, Condition) %>%
  summarise(meanRT = mean(Reaction.Time))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

## Table 2 in the main text
### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE)   #.10

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T) #.02

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT----
mean.e <- de %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(Reaction.Time))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

mean.o <- do %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(Reaction.Time))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

## Table 2 in the main text
### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE) 
KH_Correct(-0.09) #.17

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T)
KH_Correct(-0.0163) #.03

## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

#Approach 2: Z-Transformed RT difference ----

## WITHOUT TRIMMING ----
d<- d_untrimmed

# split the data into two halves
d$ItemNo = parse_number(as.character(d$ItemNo))
de = d %>%  filter(ItemNo %% 2 == 0)
do = d %>%  filter(ItemNo %% 2 == 1)


#Z-score transform RTs, compute zRT difference
## AVERAGE by item ----
mean.e <- de %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(ItemNo, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

mean.o <- do  %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(ItemNo, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE) #.36

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T) #.36


## AVERAGE by subject ----
mean.e <- de %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)
  
mean.o <- do  %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE) 
KH_Correct(-0.05) #.10

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T)
KH_Correct(-0.03) #.06

## WITH TRIMMING ----
## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

d <- d_trimmed

# split the data into two halves
d$ItemNo = parse_number(as.character(d$ItemNo))
de = d %>%  filter(ItemNo %% 2 == 0)
do = d %>%  filter(ItemNo %% 2 == 1)


#Z-score transform RTs, compute zRT difference
## AVERAGE by item ----
mean.e <- de %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(ItemNo, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

mean.o <- do  %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(ItemNo, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE) #.15

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T) #.15


## AVERAGE by subject ----
mean.e <- de %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

mean.o <- do  %>%
  group_by(Subject) %>%
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) %>% 
  select (1:5) %>%
  group_by(Subject, Condition) %>%
  summarise(meanRT = mean(zRT))%>%
  pivot_wider(names_from = Condition, values_from = meanRT) %>%
  mutate(rt.diff = Or - EitherOr)

### Pearson
print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE) 
KH_Correct(-0.08) #.15

#### Robust
WRS2::pbcor(mean.o$rt.diff, mean.e$rt.diff, ci = T)
KH_Correct(-0.092) #.17


# Approach 3: LMM----
## WITH TRIMMING ----
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

d<-d_trimmed

# Split the data in to two halves 
d$ItemNo = parse_number(as.character(d$ItemNo))
de <- d %>% filter(ItemNo %% 2 == 0)
do <- d %>% filter(ItemNo %% 2 == 1)


## fit models for even and odd numbered items (Bayesian & Generalized LMM)
library(optimx)
library(nloptr)
library(blme)

m_de = blmer(1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de, 
             control=lmerControl(optimizer = "nloptwrap",
                                 optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))

m_do = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = do, 
             control=lmerControl(optimizer = "nloptwrap",
                                 optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))

###get the slopes for all the individuals
#even-numbered
person_rand_e<-data.frame(ranef(m_de)$Subject)
person_rand_e$Subject<- row.names(person_rand_e)

#odd-numbered
person_rand_o<-data.frame(ranef(m_do)$Subject)
person_rand_o$Subject<- row.names(person_rand_o)
#combine them and create the same data frame
person <- inner_join(person_rand_e, person_rand_o, by = "Subject", suffix = c("_even", "_odd"))

#check the correlation between even- and odd-numbered sets
### Pearson
print(corr.test(person$ConditionOr_even, person$ConditionOr_odd), short = FALSE) 
KH_Correct(-.09) #.16

### Robust
WRS2::pbcor(person$ConditionOr_even, person$ConditionOr_odd, ci = T) 
KH_Correct(-.2629) #.42

## WITHOUT TRIMMING ----
## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

d<-d_untrimmed

# Split the data in to two halves 
d$ItemNo = parse_number(as.character(d$ItemNo))
de <- d %>% filter(ItemNo %% 2 == 0)
do <- d %>% filter(ItemNo %% 2 == 1)


## fit models for even and odd numbered items (Bayesian & Generalized LMM)
library(optimx)
library(blme)

m_de = blmer(1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de, 
             control=lmerControl(optimizer = "nloptwrap",
                                 optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))

m_do = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = do, 
             control=lmerControl(optimizer = "nloptwrap",
                                 optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))

###get the slopes for all the individuals
#even-numbered
person_rand_e<-data.frame(ranef(m_de)$Subject)
person_rand_e$Subject<- row.names(person_rand_e)

#odd-numbered
person_rand_o<-data.frame(ranef(m_do)$Subject)
person_rand_o$Subject<- row.names(person_rand_o)
#combine them and create the same data frame
person <- inner_join(person_rand_e, person_rand_o, by = "Subject", suffix = c("_even", "_odd"))

#check the correlation between even- and odd-numbered sets
### Pearson
print(corr.test(person$ConditionOr_even, person$ConditionOr_odd), short = FALSE) 
KH_Correct(-.11) #.20

### Robust
WRS2::pbcor(person$ConditionOr_even, person$ConditionOr_odd, ci = T) 
KH_Correct(-0.24) #.39

## With model criticism ----
## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed","de", "do", "KH_Correct", "m_de", "m_do")])

## model criticism (remove data points with abs(standardized residuals) over 2.5)
de_trim<- de[abs(scale(resid(m_de)))<2.5,]
do_trim<- do[abs(scale(resid(m_do)))<2.5,]

m_de_mc = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de_trim, 
                control=lmerControl(optimizer = "nloptwrap",
                                    optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))

m_do_mc = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = do_trim, 
                control=lmerControl(optimizer = "nloptwrap",
                                    optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))

###get the slopes for all the individuals
#even-numbered
person_rand_e<-data.frame(ranef(m_de_mc)$Subject)
person_rand_e$Subject<- row.names(person_rand_e)
#odd-numbered
person_rand_o<-data.frame(ranef(m_do_mc)$Subject)
person_rand_o$Subject<- row.names(person_rand_o)

#combine them and create the same data frame
person <- inner_join(person_rand_e, person_rand_o, by = "Subject", suffix = c("_even", "_odd"))

## Correlations
### Pearson
print(corr.test(person$ConditionOr_even, person$ConditionOr_odd), short = FALSE) #.02

### Robust
WRS2::pbcor(person$ConditionOr_even, person$ConditionOr_odd, ci = T)
KH_Correct(-.0844) #.16


### fit lmm to whole data
m <- blmer((-1/Reaction.Time)*1000 ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = d_trimmed, 
                     control=lmerControl(optimizer = "nloptwrap",
                                         optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))
summary(m)
