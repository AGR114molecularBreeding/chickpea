# Dependencies
library(dplyr)
library(stringr)
library(tidyr)

# load function for analysis
source("functions2.R")

# load dataset
load("dat/rip1.rda")


# Phenotipic Data 
# -----------------------------------------------------------------------------

# Case study 1 : mean of average of repetitions ------------------------------------------
rip %>% 
  filter(RIL %in% c(52, 67)) %>% 
  group_by(RIL, repetition) %>% 
  summarize("dias" = mean(dias)) %>% 
  ungroup() %>% 
  group_by(RIL) %>% 
  summarize(across("dias", list(mean = mean, sd = sd)))

rip %>% 
  group_by(RIL, repetition) %>% 
  summarize("dias" = mean(dias)) %>% 
  ungroup() %>% 
  group_by(RIL) %>% 
  summarize(across("dias", list(mean = mean, sd = sd)))
  arrange(desc(dias_mean))

# Case study 2 : mean of 1st day to flower per repition  -----------------------
rip %>% 
  filter(RIL %in% c(52, 67)) %>% 
  group_by(RIL, repetition) %>% 
  summarise(dias = min(dias)) %>% 
  ungroup() %>% 
  group_by(RIL) %>% 
  summarize(across("dias", list(mean = mean, sd = sd))) %>% 
  arrange(desc(dias_mean))

rip %>% 
  group_by(RIL, repetition) %>% 
  summarise(dias = min(dias)) %>% 
  ungroup() %>% 
  group_by(RIL) %>% 
  summarize(across("dias", list(mean = mean, sd = sd))) %>% 
  arrange(desc(dias_mean))


