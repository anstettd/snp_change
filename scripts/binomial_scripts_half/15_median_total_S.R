##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Author Daniel Anstett
## 
## 
## Last Modified Aug 14, 2023
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#HistPop unique SNP ascorss all env
#Updated for binomial data
obs_env_unique <- read_csv("data/binomial_data_half/slope_obs_all_unique.csv") %>% 
                  filter(SE<5.5) %>% mutate(abs_slope = abs(Slope))

#Get slope median
median_pop <- obs_env_unique %>% group_by(Site) %>% 
                                  summarise(median = median(Slope, na.rm = TRUE)) 

#Get sum of the absolute value of all slopes
abs_pop <- obs_env_unique %>% group_by(Site) %>% 
                                  summarise(abs_slope = sum(abs_slope, na.rm = TRUE))

#Get sum of all positive slopes
pos_pop <- obs_env_unique %>% filter(Slope>0) %>% group_by(Site) %>% 
                                  summarise(pos_slope = sum(Slope, na.rm = TRUE))

slope.summary <- cbind(median_pop,abs_pop[,2],pos_pop[,2])
  
write_csv(slope.summary, "data/binomial_data_half/slope_summary.csv")

#write_csv(median_pop, "data/binomial_data_half/median_pop.csv")
