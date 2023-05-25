##################################################################################
## Calc frequency for all neutral snps
## For env 1 (MAT)
## Author Daniel Anstett
## 
## 
## Last Modified April 26, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions

#Melt glm
slope_melt <- function(df) {
  freq_1.melted <- reshape2::melt(df, id.vars = c("Site", "Year"))
  colnames(freq_1.melted)[3] <- "snp_ID"
  colnames(freq_1.melted)[4] <- "freq"
  
  freq_1.melted_mod <- na.omit(freq_1.melted)
  #  freq_1.melted_mod <- freq_1.melted
  
  freq_1.slope <- group_by(freq_1.melted_mod, Site, snp_ID) %>%
    arrange(Site, Year) %>%
    summarize(Slope = glm(as.numeric(freq) ~ as.numeric(Year), family = binomial)$coefficients[2]*2,
              SE = ifelse(is.na(glm(as.numeric(freq) ~ as.numeric(Year), family = binomial)$coefficients[2]))==F,
              summary(glm(as.numeric(freq) ~ as.numeric(Year), family = binomial))$coefficients[2,2],NA)
  
  return(freq_1.slope)
}
###################################################################################
#Import
test_data <- read_csv("data/swiss_test.csv")

#Fun glm melt function
slope_melt(test_data)




