# Data pre-processing ----

## First, clear our console and environment
rm(list=ls())

## Then, set working directory
getwd() #check whether the directory is what you need. 
# If not, find the right path and change it using setwd()
#setwd("~/Desktop/Hui&Wu-R/raw_data/raw_data")
#setwd( "G:/My Drive/HuiWu_Reliability/Data_Packet/raw_data")

## load the needed dataset
d = read.csv("data_exp_35094-v18_task-uoc2.csv", fileEncoding="UTF-8-BOM")

## load packages
library(tidyr)
library(tidyverse) # Notice there are conflicts between these two packages
# when using one of the conflicted functions, remember to add specific name in the front
# e.g., instead of writing "filter()", write specifically "dplyr::filter()" or "stats::filter()"
library(psych)
library(lme4)
library(lmerTest)
library(car)
library(rstatix)

## Select and filter useful columns and rows from the dataset
## The essential columns that are needed: Participant ID, RT, Accuracy, TargetType, Condition
## Data collected during non-reaction periods need to be removed
d = d %>% select(12, 37, 38, 43, 44, 54:57) %>% 
          filter(Zone.Type == "response_keyboard_single", Primeword != "")
## Change the participantID to a shorter & uniform name "Subject"
names(d)[names(d) == "Participant.Public.ID"] <- "Subject"

## inspect false alarm rate (the rate of responding yes incorrectly)
## creating a table for different subjects with large false alarm rates
d_fa = d %>% filter(TargetType == "non") %>% group_by(Subject) %>% mutate(fa = mean (Incorrect)) %>% 
  slice (1) %>% arrange (desc(fa)) %>% filter (fa > 0.50) %>% select (1, 10)
## remove participants with an FA > 0.5 (n = 16)
d = d[!d$Subject %in% d_fa$Subject,]
## remove dataset in the environment that is not needed later
rm(d_fa)
## keep only the columns that are necessary
d = d %>% select(1,3:9)
## remove rows that have NAs
d = d %>% na.omit()

## Changing variables to factors for accurate further data analysis
d$Correct = as.factor(d$Correct)
d$Incorrect = as.factor(d$Incorrect)
d$Condition = as.factor(d$Condition)

## subset data to only include critical trials 
d = d %>% filter (TargetType == "word")

### filter correct responses by first creating an accuracy rate table for subjects and items
d.acc = d%>% select (1,3,7,8)%>%
  group_by(Subject, ItemNo) %>% summarise(acc=sum(Correct==1))

### add the accuracy rates for different subjects and items to the bigger dataset
d = left_join(d.acc, d, by = c("Subject", "ItemNo"))

## remove dataset in the environment that is not needed later
rm(d.acc)

### create trimmed and untrimmed data sets
## For untrimmed set, as long as the RT is plausible enough, it is kept.
## We are not taking out any outliers. 
## RT <= 300ms is taken out because it is technically not possible to respond that fast.  
d_untrimmed = d %>% filter (acc == 2, Reaction.Time > 300)

## For trimmed set, we remove any RT that is too long (>=2500ms)
d_trimmed = d %>% filter (acc == 2, Reaction.Time > 300, Reaction.Time < 2500)

## To get the descriptives (means & standard deviations of the two conditions) of these two datasets
d_untrimmed %>% group_by(Condition) %>% get_summary_stats(Reaction.Time, type = "mean_sd")
d_trimmed %>% group_by(Condition) %>% get_summary_stats(Reaction.Time, type = "mean_sd")
## both trimmed and untrimmed datasets show a faster RT to the related condition than the unrelated condition

## remove dataset in the environment that is not needed later
rm(d)

## write a function to apply the correction according to Krus and Helmstadter 1993
KH_Correct <- function(r) {
  corrected <- r*-1 / (0.5 * (1-r))
    return(corrected)
}

# Approach 1: Raw RT difference ----
## WITHOUT TRIMMING ----

# make the name simpler
d <- d_untrimmed

## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
# we are using the splithalf method to check the reliability
# spliting the items into 2 halves and compare the two
d$ItemNo = parse_number(d$ItemNo) # get the pure numbers out
de = d %>% filter(ItemNo %% 2 == 0) # even-numbered items
do = d %>% filter(ItemNo %% 2 == 1) # odd-numbered items

# inverse transform the data set so we can easily compute RT difference
# then we remove NAs
de <- de |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time, id_cols = c(1:2)) |>
  mutate(rt.diff = unrelated - related) |>
  filter(!is.na(rt.diff))

do <- do |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time, id_cols = c(1:2)) |>
  mutate(rt.diff = unrelated - related) |>
  filter(!is.na(rt.diff))

### AVERAGE BY ITEM ----
# get the mean values of RT difference by items
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(rt.diff))

## Table 2 in the main text

### Pearson
# reliability calculation - basic correlation test (-0.06)
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = FALSE)   
# correction of the basic test
KH_Correct(-0.06) #.11

#### Robust (get the robust correlation)
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) 
KH_Correct(-0.0213) #.04

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT----
# get the mean values of RT differences by participants
mean.e = de %>% group_by(Subject) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(rt.diff))

#row numbers do not match, so we remove the 2 extra rows in the odd-numbered results
anti_join(mean.o, mean.e, by = c("Subject" = "Subject"))
mean.o = mean.o[-c(47, 120),]

## Table 2 in the main text
### Pearson
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = FALSE) # .07 

#### Robust
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) # .08

rm(d, de, do)

#### Plotting ----

## create a dataset that contains means of both odd-numbered & even-numbered items' RTs
raw_untrimmed<-inner_join(x= mean.o, y=mean.e, by="Subject")
names(raw_untrimmed) <- c("Subject", "Odd", "Even")


## scatter plot
ggplot(data=raw_untrimmed, aes(x = Even, y =  Odd)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  scale_x_continuous(name =  "Raw RT Differences \n (Even-Numbered Sub-Dataset)", limits = c(-250,250)) +
  scale_y_continuous(name =  "Raw RT Differences \n (Odd-Numbered Sub-Dataset)", limits = c(-250,250))

## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

## WITH TRIMMING ----
# same procedures as above
d <- d_trimmed

## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
d$ItemNo = parse_number(d$ItemNo)
de <- d %>% filter(ItemNo %% 2 == 0)
do = d %>% filter(ItemNo %% 2 == 1)

# turn it wide, compute raw difference, and remove NAs
de <- de |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time, id_cols = c(1:2)) |>
  mutate(rt.diff = unrelated - related) |>
  filter(!is.na(rt.diff))

# turn it wide, compute raw difference, and remove NAs
do <- do |> 
  pivot_wider(names_from = Condition, values_from = Reaction.Time, id_cols = c(1:2)) |>
  mutate(rt.diff = unrelated - related) |>
  filter(!is.na(rt.diff))

### AVERAGE BY ITEM ----
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(rt.diff))

### Pearson 
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = FALSE)  #.36

### Robust
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) # .23

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT ----
mean.e = de %>% group_by(Subject) %>% summarise(mean(rt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(rt.diff))

#row numbers do not match, so we remove 3 extra rows in the odd-numbered results
anti_join(mean.o, mean.e, by = c("Subject" = "Subject"))
mean.o = mean.o[-c(47, 71, 120),]

### Pearson 
print(corr.test(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`), short = F)
KH_Correct(-0.06) # .11

### Robust
WRS2::pbcor(mean.o$`mean(rt.diff)`, mean.e$`mean(rt.diff)`, ci = T) 
KH_Correct(-0.0927) # .17


#Approach 2: Z-Transformed RT difference ----

## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

## WITHOUT TRIMMING ----
# similar procedure as above other than the z-transformation
 d<- d_untrimmed

# split the data into two halves
d$ItemNo = parse_number(d$ItemNo)
de = d %>%  filter(ItemNo %% 2 == 0)
do = d %>%  filter(ItemNo %% 2 == 1)


#Z-score transform RTs, compute zRT difference, remove NAs for the 2 groups
de <- de |>
  group_by(Subject) |>
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,9,10) |>
  pivot_wider(names_from = Condition, values_from = zRT, id_cols = c(1:2)) |>
  mutate(zrt.diff = unrelated - related) |>
  filter(!is.na(zrt.diff))

do <- do |> 
  group_by(Subject) |>
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,9,10) |>
  pivot_wider(names_from = Condition, values_from = zRT, id_cols = c(1:2)) |>
  mutate(zrt.diff = unrelated - related) |>
  filter(!is.na(zrt.diff))

### AVERAGE BY ITEM ----
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))

#remove a row because of NA
mean.o = mean.o[-c(9),]
mean.e = mean.e[-c(9),]

### Pearson
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) #.15

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T) # .04

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT ----
mean.e = de %>% group_by(Subject) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(zrt.diff))

#again, remove 3 extra rows in the odd-numbered results
anti_join(mean.o, mean.e, by = c("Subject" = "Subject"))
mean.o = mean.o[-c(47, 120),]

### Pearson 
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) #.15

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T) # .15

## WITH TRIMMING ----
## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

d <- d_trimmed

## SPLITHALF AGAIN INTO EVEN- & ODD-NUMBERED ITEMS
d$ItemNo = parse_number(d$ItemNo)
de = d %>% filter(ItemNo %% 2 == 0)
do = d %>% filter(ItemNo %% 2 == 1)

#Z-score transform RTs, compute zRT difference, remove NAs for the 2 groups
de <- de |>
  group_by(Subject) |>
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,9,10) |>
  pivot_wider(names_from = Condition, values_from = zRT, id_cols = c(1:2)) |>
  mutate(zrt.diff = unrelated - related) |>
  filter(!is.na(zrt.diff))

# turn it wide, compute raw difference, and remove NAs
do <- do |>
  group_by(Subject) |>
  mutate(zRT = scale(Reaction.Time, center = T, scale = T), 
          meanRT = mean(Reaction.Time)) |> 
  select (1,2,9,10) |>
  pivot_wider(names_from = Condition, values_from = zRT, id_cols = c(1:2)) |>
  mutate(zrt.diff = unrelated - related) |>
  filter(!is.na(zrt.diff))

### AVERAGE BY ITEM ----
mean.e = de %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(ItemNo) %>% summarise(mean(zrt.diff))
#remove a row because of NA
mean.o = mean.o[-c(9),]
mean.e = mean.e[-c(9),]

### Pearson
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) #.02

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T)
KH_Correct(-0.0598) # .11

rm(mean.e, mean.o)

### AVERAGE BY PARTICIPANT ----
mean.e = de %>% group_by(Subject) %>% summarise(mean(zrt.diff))
mean.o = do %>% group_by(Subject) %>% summarise(mean(zrt.diff))

#again, remove 3 extra rows in the odd-numbered results
mean.o = mean.o[-c(47, 71, 120),]

### Pearson
print(corr.test(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`), short = FALSE) #.03

### Robust
WRS2::pbcor(mean.o$`mean(zrt.diff)`, mean.e$`mean(zrt.diff)`, ci = T) # .01

# Approach 3: Linear Mixed-Effects Modeling (LMM)----


## WITH TRIMMING ----
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

d<-d_trimmed

# Split the data in to two halves 
d$ItemNo = parse_number(d$ItemNo)
de <- d %>% filter(ItemNo %% 2 == 0)
do <- d %>% filter(ItemNo %% 2 == 1)

## fit models for even and odd numbered items (Bayesian & Generalized LMM)
library(optimx)
library(blme)

m_de = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de, 
             control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

m_do = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = do, 
             control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

###get the slopes for all the individuals
#even-numbered
person_rand_e<-data.frame(ranef(m_de)$Subject)
person_rand_e$Subject<- row.names(person_rand_e)

#odd-numbered
person_rand_o<-data.frame(ranef(m_do)$Subject)
person_rand_o$Subject<- row.names(person_rand_o)

#combine them into one data frame
person <- inner_join(person_rand_e, person_rand_o, by = "Subject", suffix = c("_even", "_odd"))

#check the correlation between even- and odd-numbered sets
### Pearson
print(corr.test(person$Conditionunrelated_even, person$Conditionunrelated_odd), short = FALSE) #.85

### Robust
WRS2::pbcor(person$Conditionunrelated_even, person$Conditionunrelated_odd, ci = T) # .81

## WITHOUT TRIMMING ----
## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed", "KH_Correct")])

d<-d_untrimmed

# Split the data in to two halves 
d$ItemNo = parse_number(d$ItemNo)
de <- d %>% filter(ItemNo %% 2 == 0)
do <- d %>% filter(ItemNo %% 2 == 1)

## fit models for even and odd numbered items (Bayesian & Generalized LMM)
library(optimx)
library(blme)

m_de = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de, 
             control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

m_do = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = do, 
             control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

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
print(corr.test(person$Conditionunrelated_even, person$Conditionunrelated_odd), short = FALSE) #.81

### Robust
WRS2::pbcor(person$Conditionunrelated_even, person$Conditionunrelated_odd, ci = T) # .79

#### Plotting ----
library (ggplot2)

## scatter plot with random slopes
ggplot(data=person, aes(x = Conditionunrelated_even, y = Conditionunrelated_odd)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  scale_x_continuous(name =  "Participant Random Slopes \n (Even-Numbered Sub-Dataset)") +
  scale_y_continuous(name =  "Participant Random Slopes \n (Odd-Numbered Sub-Dataset)")

## With model criticism ----
## clean the environment
rm(list=ls()[! ls() %in% c("d_trimmed", "d_untrimmed","de", "do", "KH_Correct", "m_de", "m_do")])

## model criticism (remove data points with abs(standardized residuals) over 2.5)
de_trim<- de[abs(scale(resid(m_de)))<2.5,]
do_trim<- do[abs(scale(resid(m_do)))<2.5,]

m_de_mc = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = de_trim, 
             control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

m_do_mc = blmer(-1/Reaction.Time ~ Condition + (1+Condition|Subject) + (1+Condition|ItemNo), data = do_trim, 
             control=lmerControl(optimizer="optimx",optCtrl=list(method="nlminb")))

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
print(corr.test(person$Conditionunrelated_even, person$Conditionunrelated_odd), short = FALSE) #.88

### Robust
WRS2::pbcor(person$Conditionunrelated_even, person$Conditionunrelated_odd, ci = T) # .87

## put boxplots together
library(ggpubr)

ggarrange(z_o, p_lmm_mc_o, z_e, p_lmm_mc_e)

 