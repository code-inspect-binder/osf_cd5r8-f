# Data pre-processing ----

## Clear console and environment
rm(list=ls())

## set working directory
getwd()
#setwd("~/Desktop/Hui&Wu-R/raw_data")
#setwd( "G:/My Drive/HuiWu_Reliability/Data_Packet/raw_data")

## load packages
library(tidyr)
library(tidyverse)
library(psych)
library(lme4)
library(lmerTest)
library(car)
library(rstatix)
## because of conflicts between different packages, we are address them by setting our preferences
library(conflicted)
conflict_prefer("rename", "dplyr")

## Pre-process data according to Buffinton's script
## load data (both ASRT and participant removal data)

DS <- read.csv('./raw_data/MasterData_OutliersRemoved.csv')
head(DS)

## use this variable in case data need to be subset for participants in the final dataset.
FinalIDs <- DS$ID

rm(DS)

d <- read.csv("./raw_data/ASRT_MasterData.csv")
head(d)

## screen participants (only keep those that are in the DS dataset)
## filter out practice trials, code trial types

d <- subset(d, (Subject %in% FinalIDs) & (TrialType == 'P' | TrialType == 'R'),
                     select = c(Subject, TrialType, firstRT))

## Assign itemNo. for mixed-effects model
d$ItemNo <- rep(1:800, each = 2)


## create trimmed and untrimmed data for our own analysis
d_trimmed <- d %>%
            group_by(Subject) %>%
            filter(firstRT > 100) %>%
            mutate(UpperOutlier = mean(firstRT) + 3*sd(firstRT)) %>%
            filter(firstRT < UpperOutlier) %>%
            select (1:4)%>%
           rename (Condition = TrialType, Reaction.Time = firstRT)

d_untrimmed<-d %>%
              group_by(Subject) %>%
              filter(firstRT > 100) %>%
              rename (Condition = TrialType, Reaction.Time = firstRT)

## write a function to apply the correction according to Krus and Helmstadter 1993
KH_Correct <- function(r) {
  corrected <- r*-1 / (0.5 * (1-r))
    return(corrected)
}

## clearn the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])


# all similar to procedures in the data analysis of the other two data sets
# Approach 1: Raw RT difference ----
## WITHOUT TRIMMING ----

d <- d_untrimmed

## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
d$ItemNo = parse_number(as.character(d$ItemNo))
de = d %>% filter(ItemNo %% 2 == 0)
do = d %>% filter(ItemNo %% 2 == 1)

#compute RT difference, remove NAs
de <- de |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time) |>
  mutate (rt.diff = R - P) |>
  filter(!is.na(rt.diff))

do <- do |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time) |>
  mutate (rt.diff = R - P) |>
  filter(!is.na(rt.diff))

### AVERAGE BY ITEM ----
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(rt.diff))

## Table 2 in the main text
### Pearson
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = FALSE)   #.01

#### Robust
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) # .05


rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT----
mean.e = de %>% group_by(Subject) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(rt.diff))


## Table 2 in the main text
### Pearson
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = FALSE) # .06 

#### Robust
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) #.44
KH_Correct(-0.2818)

## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

## WITH TRIMMING ----
d <- d_trimmed

## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
d$ItemNo = parse_number(as.character(d$ItemNo))
de <- d %>% filter(ItemNo %% 2 == 0)
do = d %>% filter(ItemNo %% 2 == 1)

# turn it wide, compute raw difference, and remove NAs
de <- de |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time) |>
  mutate (rt.diff = R - P) |>
  filter(!is.na(rt.diff))

# turn it wide, compute raw difference, and remove NAs
do <- do |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time) |>
  mutate (rt.diff = R - P) |>
  filter(!is.na(rt.diff))

### AVERAGE BY ITEM ----
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(rt.diff))

### Pearson 
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = FALSE)  #.02

### Robust
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) 
KH_Correct(-0.0066) #.01

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT ----
mean.e = de %>% group_by(Subject) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(rt.diff))

### Pearson 
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = F)
KH_Correct(-0.31) # .47

### Robust
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) 
KH_Correct(-0.3426) # .51


#Approach 2: Z-Transformed RT difference ----

## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

## WITHOUT TRIMMING ----
 d<- d_untrimmed

# split the data into two halves
d$ItemNo = parse_number(as.character(d$ItemNo))
de = d %>%  filter(ItemNo %% 2 == 0)
do = d %>%  filter(ItemNo %% 2 == 1)


#Z-score transform RTs, compute zRT difference, remove NAs for the 2 groups
de <- de |>
  group_by(Subject) |>
  mutate (zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,4,5) |>
  pivot_wider(names_from = Condition, values_from = zRT) |>
  mutate (zrt.diff = R-P) |>
  filter(!is.na(zrt.diff))

do <- do |> 
  group_by(Subject) |>
  mutate (zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,4,5) |>
  pivot_wider(names_from = Condition, values_from = zRT) |>
  mutate (zrt.diff = R-P) |>
  filter(!is.na(zrt.diff))

### AVERAGE BY ITEM ----
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))

### Pearson
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) #.09

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T) # .10

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT ----
mean.e = de %>% group_by(Subject) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(zrt.diff))

### Pearson 
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) 
KH_Correct(-0.32) #.48

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T) 
KH_Correct(-0.2983) #.46


## WITH TRIMMING ----
## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

d <- d_trimmed

## SPLITHALF AGAIN INTO EVEN- & ODD-NUMBERED ITEMS
d$ItemNo = parse_number(as.character(d$ItemNo))
de = d %>% filter(ItemNo %% 2 == 0)
do = d %>% filter(ItemNo %% 2 == 1)

#Z-score transform RTs, compute zRT difference, remove NAs for the 2 groups
de <- de |>
  group_by(Subject) |>
  mutate (zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,4,5) |>
  pivot_wider(names_from = Condition, values_from = zRT) |>
  mutate (zrt.diff = R-P) |>
  filter(!is.na(zrt.diff))

# turn it wide, compute raw difference, and remove NAs
do <- do |>
  group_by(Subject) |>
  mutate (zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,4,5) |>
  pivot_wider(names_from = Condition, values_from = zRT) |>
  mutate (zrt.diff = R-P) |>
  filter(!is.na(zrt.diff))

### AVERAGE BY ITEM ----
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))

### Pearson
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) #.07

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T) #.06

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT ----
mean.e = de %>% group_by(Subject) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(zrt.diff))

### Pearson
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) 
KH_Correct(-0.34) #.51

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T) # .01
KH_Correct(-0.3075) #


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

m_de = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de, 
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
print(corr.test(person$ConditionR_even, person$ConditionR_odd), short = FALSE)
KH_Correct(-.28) #.44

### Robust
WRS2::pbcor(person$ConditionR_even, person$ConditionR_odd, ci = T) #.30

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

m_de = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de, 
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
print(corr.test(person$ConditionR_even, person$ConditionR_odd), short = FALSE) 
KH_Correct(-.31) #.47

### Robust
WRS2::pbcor(person$ConditionR_even, person$ConditionR_odd, ci = T) #.30

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
print(corr.test(person$ConditionR_even, person$ConditionR_odd), short = FALSE)
KH_Correct(-.33) #.50

### Robust
WRS2::pbcor(person$ConditionR_even, person$ConditionR_odd, ci = T)
KH_Correct(-.3249) #.49


# fit lmm to the whole dataset
m <- blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = d_trimmed, 
      control=lmerControl(optimizer = "nloptwrap",
                          optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 2e5)))
summary(m)
