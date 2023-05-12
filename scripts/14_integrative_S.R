##################################################################################
## Get cumulative slope measure
## 
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
###################################################################################
#Import libraries
library(tidyverse)
library(car) 

###################################################################################
#Import files
snp_list_unique <- read.csv("data/snp_list_unique.csv")

#snp_list <- read.csv("Genomics_scripts/Data/snp_list.csv")
#snp_list_mat <- snp_list  %>% filter(env=="MAT")
#snp_list_map <- snp_list  %>% filter(env=="MAP")
#snp_list_cmd <- snp_list  %>% filter(env=="CMD")

#offset_pop <- read_csv("Genomics_scripts/Data/offset_pop_lambda.csv")

###################################################################################
#Produce integrative slope

slope.cumul <- data.frame() 
cumul_unique <- snp_list_unique %>% group_by(Site) %>% summarise(cumul_slope = sum(Slope))

#Populate blank dataframe 
#slope.cumul[1,1] <- cumul_unique$cumul_slope[1]
#slope.cumul[2,1] <- 0
#slope.cumul[3,1] <- cumul_unique$cumul_slope[2]
#slope.cumul[4,1] <- cumul_unique$cumul_slope[3]
#slope.cumul[5,1] <- cumul_unique$cumul_slope[4]
#slope.cumul[6,1] <- 0
#slope.cumul[7,1] <- 0
#slope.cumul[8,1] <- 0
#slope.cumul[9,1] <- cumul_unique$cumul_slope[5]
#slope.cumul[10,1] <- 0
#slope.cumul[11,1] <- cumul_unique$cumul_slope[6]
#slope.cumul[12,1] <- 0

#colnames(slope.cumul) <- "cumul_slope"
#slope.cumul <- slope.cumul %>% mutate(Paper_ID = row_number())

###################################################################################
## Integrate with demographic data

offset_pop <- read_csv("data/offset_pop_9var.csv") %>% filter(Paper_ID<13)
demo_cumul <- left_join(offset_pop,cumul_unique,by=c("Paper_ID"="Site"))

##Integrate with timeseries
timeseries <- read_csv("data/offset_pop_timeseries.csv")
time_cumul <-left_join(timeseries,cumul_unique,by=c("Paper_ID"="Site"))
time_cumul$Region <- c("South",
                       "South",
                       "Center",
                       "Center",
                       "Center",
                       "Center",
                       "Center",
                       "North",
                       "North",
                       "North",
                       "North",
                       "South")


#Export
#write_csv(demo_cumul,"data/demo_cumul.csv")
write_csv(time_cumul,"data/time_cumul.csv")





