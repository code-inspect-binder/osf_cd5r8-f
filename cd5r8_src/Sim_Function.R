### We wrote this script based on the useful tutorial by Lisa DeBruine
### https://debruine.github.io/tutorials/sim-lmer.html#rslope-function

rm(list=ls())

#Setup
library(tidyverse) # for data wrangling, pipes, and good dataviz
library(lmerTest)  # for mixed effect models
library(GGally)    # makes it easy to plot relationships between variables
library(faux)      # for simulating correlated variables
library(psych)
library(optimx)
library(blme)

library(conflicted)
conflict_prefer("recode", "dplyr")
conflict_prefer("filter", "dplyr")

options("scipen"=10, "digits"=4) # control scientific notation in lmer output
## Setting a seed is common in data simulation
## With the same number, you can make sure you get the same results every time you run the script
set.seed(852) 

# simulate data ----

# we first create a practice function for simulation of data of subjects, stimuli & trials

# this function is from Line 28 to 105

sim <- function (

# specify parameters
sub_n            = 100, # number of subjects in this simulation
sub_sd           = .04, # SD (standard deviation) for the subjects' random intercept
stim_n           = 800,  # number of stimuli in this simulation
stim_sd          = .02,  # SD for the stimuli's random intercept
grand_i          = .30, # overall mean DV (dependent variables)
stim_version_eff = -.0005,  # mean difference between versions: related - unrelated
error_sd         = .1,  # residual (error) SD

sub_version_sd = .006,  # subject slopes SD
sub_i_version_cor = -0.02, # subject intercept-slope correlation

stim_version_sd = .001,  # SD for the stimuli's random slope for stim_version
stim_i_cor = -0.3 # correlations between intercept and slopes

) {
## subjects 
sub = faux::rnorm_multi(
  n = sub_n, 
  vars = 2, # number of variables
  r = sub_i_version_cor, # correlations among variables
  mu = 0, # means of random intercepts and slopes are always 0
  sd = c(sub_sd, sub_version_sd),
  varnames = c("sub_i", "sub_version_slope") # naming the variables (subject intercepts & slopes)
) %>%
  mutate(
    sub_id = 1:sub_n) # creating the ID No. from 1 to 100

## stimuli 
# specify correlations for rnorm_multi (one of several methods)
stim_cors = stim_i_cor
stim = rnorm_multi(
  n = stim_n, 
  vars = 2, 
  r = stim_cors, 
  mu = 0, # means of random intercepts and slopes are always 0
  sd = c(stim_sd, stim_version_sd),
  varnames = c("stim_i", "stim_version_slope")
) %>%
  mutate(
    stim_id = 1:stim_n
  )


## trials 
trials = crossing(
  sub_id = sub$sub_id, # get subject IDs from the sub data table
  stim_id = stim$stim_id, # get stimulus IDs from the stim data table
  stim_version = c("related", "unrelated") # all subjects see both related and unrelated versions of all stimuli
) %>%
  left_join(sub, by = "sub_id") %>% # includes the intercept, slope, and condition for each subject
  left_join(stim, by = "stim_id")   # includes the intercept and slopes for each stimulus

# compute the effects of trials on the overall results and the error terms 
dat = trials %>%
  mutate(
    # effect-code subject condition and stimulus version
    stim_version.e = recode(stim_version, "related" = 1, "unrelated" = 0),
    # calculate trial-specific effects by adding overall effects and slopes
    version_eff = stim_version_eff + stim_version_slope + sub_version_slope,
    # calculate error term (normally distributed residual with SD set above)
    err = rnorm(nrow(.), 0, error_sd),
    # calculate DV from intercepts, effects, and error
    dv = grand_i + sub_i + stim_i + err +
      (stim_version.e * version_eff) 
  )

## lmm analysis 
mod = blmer(dv ~ stim_version + (1+stim_version|sub_id) + (1+stim_version|stim_id), data = dat, 
            control=lmerControl(optimizer = "nloptwrap",
                              optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 5e6)))

mod.sum = summary(mod)

return(mod.sum)
}


# ADD RAW SCORE RELIABILITY ----
# first few components stay the same to create 
# then split into 2 halves and calculate the raw differences to check the reliability
# this function runs from Line 113 to Line 209

sim_raw <- function (
    # specify parameters
  sub_n            = 120, # number of subjects in this simulation
  sub_sd           = 0.29, # SD for the subjects' random intercept
  stim_n           = 40,  # number of stimuli in this simulation
  stim_sd          = 0.08,  # SD for the stimuli's random intercept
  grand_i          = -1.5, # overall mean DV
  stim_version_eff = 0.1,  # mean difference between versions: related - unrelated
  error_sd         = 0.2,  # residual (error) SD
  
  sub_version_sd = 0.03,  # subject slopes SD
  sub_i_version_cor = -0.97, # subject intercept-slope correlation
  
  stim_version_sd = 0.03,  # SD for the stimuli's random slope for stim_version
  stim_i_cor = -0.18 # correlations between intercept and slopes
) {
  set.seed(852)   
  ## subjects 
    sub = faux::rnorm_multi(
      n = sub_n, 
      vars = 2, 
      r = sub_i_version_cor,
      mu = 0, # means of random intercepts and slopes are always 0
      sd = c(sub_sd, sub_version_sd),
      varnames = c("sub_i", "sub_version_slope")
    ) %>%
      mutate(
        sub_id = 1:sub_n)
    
    ## stimuli 
    # specify correlations for rnorm_multi (one of several methods)
    stim_cors = stim_i_cor
    stim = rnorm_multi(
      n = stim_n, 
      vars = 2, 
      r = stim_cors, 
      mu = 0, # means of random intercepts and slopes are always 0
      sd = c(stim_sd, stim_version_sd),
      varnames = c("stim_i", "stim_version_slope")
    ) %>%
      mutate(
        stim_id = 1:stim_n
      )
    
    
    ## trials 
    trials = crossing(
      sub_id = sub$sub_id, # get subject IDs from the sub data table
      stim_id = stim$stim_id, # get stimulus IDs from the stim data table
      stim_version = c("related", "unrelated") # all subjects see both related and unrelated versions of all stimuli
    ) %>%
      left_join(sub, by = "sub_id") %>% # includes the intercept, slope, and condition for each subject
      left_join(stim, by = "stim_id")   # includes the intercept and slopes for each stimulus
    
    # compute dat = trials %>%
    dat = trials %>%
      mutate(
        # effect-code subject condition and stimulus version
        stim_version.e = recode(stim_version, "related" = 1, "unrelated" = 0),
        # calculate trial-specific effects by adding overall effects and slopes
        version_eff = stim_version_eff + stim_version_slope + sub_version_slope,
        # calculate error term (normally distributed residual with SD set above)
        err = rnorm(nrow(.), 0, error_sd),
        # calculate DV from intercepts, effects, and error
        dv = grand_i + sub_i + stim_i + err +
          (stim_version.e * version_eff) 
      )
## split half reliability 
## create another data set with a different but short name 
## for easier coding & keeping the original set the same
d<-dat
    
## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
d$stim_id = parse_number(as.character(d$stim_id))
de <- d %>% filter(stim_id %% 2 == 0)
do = d %>% filter(stim_id %% 2 == 1)
    
#### AVERAGE BY PARTICIPANT 
mean.e <- de%>%
      group_by(sub_id, stim_version) %>%
      summarise(meanRT = mean(dv))%>%
      pivot_wider(names_from = stim_version, values_from = meanRT) %>%
      mutate(rt.diff = related - unrelated)
    
mean.o <- do%>%
      group_by(sub_id, stim_version) %>%
      summarise(meanRT = mean(dv))%>%
      pivot_wider(names_from = stim_version, values_from = meanRT) %>%
      mutate(rt.diff = related - unrelated)
    
### Pearson


out<-print(corr.test(mean.o$rt.diff, mean.e$rt.diff), short = FALSE)
    
return(out)
}



# ADD MODEL_BASED RELIABILITY ----
# this function runs from Line 213 to Line 315

sim_model_based <- function (
   
  # specify parameters
  sub_n            = 120, # number of subjects in this simulation
  sub_sd           = 0.29, # SD for the subjects' random intercept
  stim_n           = 40,  # number of stimuli in this simulation
  stim_sd          = 0.08,  # SD for the stimuli's random intercept
  grand_i          = -1.5, # overall mean DV
  stim_version_eff = 0.1,  # mean difference between versions: related - unrelated
  error_sd         = 0.2,  # residual (error) SD
  
  sub_version_sd = 0.03,  # subject slopes SD
  sub_i_version_cor = -0.97, # subject intercept-slope correlation
  
  stim_version_sd = 0.03,  # SD for the stimuli's random slope for stim_version
  stim_i_cor = -0.18 # correlations between intercept and slopes
) {
  set.seed(860) 
   ## subjects 
  sub = faux::rnorm_multi(
    n = sub_n, 
    vars = 2, 
    r = sub_i_version_cor,
    mu = 0, # means of random intercepts and slopes are always 0
    sd = c(sub_sd, sub_version_sd),
    varnames = c("sub_i", "sub_version_slope")
  ) %>%
    mutate(
      sub_id = 1:sub_n)
  
  ## stimuli 
  # specify correlations for rnorm_multi (one of several methods)
  stim_cors = stim_i_cor
  stim = rnorm_multi(
    n = stim_n, 
    vars = 2, 
    r = stim_cors, 
    mu = 0, # means of random intercepts and slopes are always 0
    sd = c(stim_sd, stim_version_sd),
    varnames = c("stim_i", "stim_version_slope")
  ) %>%
    mutate(
      stim_id = 1:stim_n
    )
  
  
  ## trials 
  trials = crossing(
    sub_id = sub$sub_id, # get subject IDs from the sub data table
    stim_id = stim$stim_id, # get stimulus IDs from the stim data table
    stim_version = c("related", "unrelated") # all subjects see both related and unrelated versions of all stimuli
  ) %>%
    left_join(sub, by = "sub_id") %>% # includes the intercept, slope, and condition for each subject
    left_join(stim, by = "stim_id")   # includes the intercept and slopes for each stimulus
  
  # compute dat = trials %>%
  dat = trials %>%
    mutate(
      # effect-code subject condition and stimulus version
      stim_version.e = recode(stim_version, "related" = 1, "unrelated" = 0),
      # calculate trial-specific effects by adding overall effects and slopes
      version_eff = stim_version_eff + stim_version_slope + sub_version_slope,
      # calculate error term (normally distributed residual with SD set above)
      err = rnorm(nrow(.), 0, error_sd),
      # calculate DV from intercepts, effects, and error
      dv = grand_i + sub_i + stim_i + err +
        (stim_version.e * version_eff) 
    )
  ## split half reliability  
  d<-dat
  
  ## SPLITHALF INTO EVEN- & ODD-NUMBERED ITEMS
  d$stim_id = parse_number(as.character(d$stim_id))
  de = d %>% filter(stim_id %% 2 == 0)
  do = d %>% filter(stim_id %% 2 == 1)
  
  ## build models
  m_de = blmer(dv ~ stim_version + (1+stim_version|sub_id) + (1+stim_version|stim_id), data = de, 
               control=lmerControl(optimizer = "nloptwrap",
                                   optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 5e8)))
  
  m_do = blmer(dv ~ stim_version + (1+stim_version|sub_id) + (1+stim_version|stim_id), data = do, 
               control=lmerControl(optimizer = "nloptwrap",
                                   optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", maxit = 5e8)))
  
  ###get the slopes for all the individuals
  #even-numbered
  person_rand_e<-data.frame(ranef(m_de)$sub_id)
  person_rand_e$sub_id<- row.names(person_rand_e)
  
  #odd-numbered
  person_rand_o<-data.frame(ranef(m_do)$sub_id)
  person_rand_o$sub_id<- row.names(person_rand_o)
  #combine them and create the same data frame
  person <- inner_join(person_rand_e, person_rand_o, by = "sub_id", suffix = c("_even", "_odd"))
  
  #check the correlation between even- and odd-numbered sets
  ### Pearson
  out<-print(corr.test(person$stim_versionunrelated_even, person$stim_versionunrelated_odd), short = FALSE) 
  
  
  return(out)
}


# baseline (based on Hui et al) ----
# run the functions to check the results (correlation matrix & confidence intervals)
# the baseline number of stimuli is 40
sim_raw() #.23 [.05, .39]
sim_model_based() #.95 [.93, .97]
# we can see a clear advantage for model-based simulation


# but maybe the advantage is based on the number of stimuli
# what if we change the number of stimuli from 40 to 20
# K changing ----
sim_raw(stim_n = 20)  #.03 [-.15, .21]
sim_model_based(stim_n = 20) #.79 [.71, .85]
# both number sets have decreased
# but the results of the model-based analysis is still quite impressive and still shows a strong reliability


# let's try another number: 30
sim_raw(stim_n = 30)  #.12 [-.06, .29]
sim_model_based(stim_n = 30) #.90 [.86, .93]
# both number sets have increased a little but not much
# similarly, the model-based analysis still returns a strong reliability


# what about bigger numbers of stimuli? 
# How would a significant increase of stimuli number affect the reliability from the 2 methods?

sim_raw(stim_n = 80)  #.32 [.15, .47]
sim_model_based(stim_n = 80) #.95 [.93, .97]

sim_raw(stim_n = 120)  #.35 [.18, .50]
sim_model_based(stim_n = 120) #.98 [.97, .99]

sim_raw(stim_n = 160)  #.50 [.35, .62]
sim_model_based(stim_n = 160) #.98 [.98, .99] # convergence issues

sim_raw(stim_n = 400)  #.69 [.59, .78]
sim_model_based(stim_n = 400) #.98 [.97, .99] # convergence issues

sim_raw(stim_n = 600)  #.79 [.72, .85]
sim_model_based(stim_n = 600) #.99 [.99, .99] # convergence issues

## as shown in the above results
## when the number of items are much bigger, the basic method of raw differences is also quite powerful
## even though model based method returns better results, raw difference comparison also appears to be reliable enough


# changing levels of noise
### the other possible factor at play is the level of noise, or error term
# the baseline level of noise is 0.2
# let's try smaller and bigger ones

sim_raw(error_sd = 0.05)  # .81 [.73, .86]
sim_model_based(error_sd = 0.05) # .99 [.98, .99] # convergence issues

sim_raw(error_sd = 0.1)  #.52 [.38, .64]
sim_model_based(error_sd = 0.1) # .98 [.98, .99]

# smaller error terms return better reliability levels
# but they affect the raw difference results much more than the model-based method


# Now we try the other direction

sim_raw(error_sd = 0.4)  #.07 [-.11, .25]
sim_model_based(error_sd = 0.4) # .82[.76, .87]

sim_raw(error_sd = 0.6)  #.03 [-.15, .21]
sim_model_based(error_sd = 0.6) #.63 [.50, .72] 

sim_raw(error_sd = 0.8)  #.02 [-.16, .19]
sim_model_based(error_sd = 0.8) #.39 [.23, .53] 

sim_raw(error_sd = 1)  #.01 [-.17, .19]
sim_model_based(error_sd = 1) #.18 [.00. .35]

# Once the error term reaches over 0.2, the reliability plummets for the raw difference method
# method-based approach can handle bigger error terms much better, maintaining a reasonable level of reliability till after the error term reaches over 0.8

