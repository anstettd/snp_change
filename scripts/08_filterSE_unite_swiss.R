##################################################################################
## Remove high SE values for all snp glm in each 100 k region
## 
## Author Daniel Anstett
## 
## 
## Last Modified May 3, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Import and filter for high SE

swiss_glm_1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_1.csv")

write_csv(swiss_glm_1,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_1.csv")






