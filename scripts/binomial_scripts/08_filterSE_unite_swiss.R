##########################################################################################################
## Remove high SE values for all snp glm in each 100 k region
## 
## Author Daniel Anstett
## 
## 
## Last Modified May 3, 2023
##########################################################################################################
#Import libraries
library(tidyverse)

##########################################################################################################
#Import and filter for high SE

swiss_glm_1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_1.csv")
swiss_glm_2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_2.csv")
swiss_glm_3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_3.csv")
swiss_glm_4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_4.csv")
swiss_glm_5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_5.csv")
swiss_glm_6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_6.csv")
swiss_glm_7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_7.csv")
swiss_glm_8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_8.csv")
swiss_glm_9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_9.csv")
swiss_glm_10 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_10.csv")
swiss_glm_11 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_11.csv")
swiss_glm_12 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_12.csv")
swiss_glm_13 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_13.csv")
swiss_glm_14 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_14.csv")
swiss_glm_15 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_15.csv")
swiss_glm_16 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_16.csv")
swiss_glm_17 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_17.csv")
swiss_glm_18 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_18.csv")
swiss_glm_19 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_19.csv")
swiss_glm_20 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_20.csv")

#Bind glms
swiss_glm <- rbind(swiss_glm_1,
                   swiss_glm_2,
                   swiss_glm_3,
                   swiss_glm_4,
                   swiss_glm_5,
                   swiss_glm_6,
                   swiss_glm_7,
                   swiss_glm_8,
                   swiss_glm_9,
                   swiss_glm_10,
                   swiss_glm_11,
                   swiss_glm_12,
                   swiss_glm_13,
                   swiss_glm_14,
                   swiss_glm_15,
                   swiss_glm_16,
                   swiss_glm_17,
                   swiss_glm_18,
                   swiss_glm_19,
                   swiss_glm_20)

#Filter for SE < 5
swiss_glm_filter <- swiss_glm %>% filter(SE<5.5)

dim(swiss_glm)
dim(swiss_glm_filter)

#Export
write_csv(swiss_glm_filter,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_filter.csv")
write_csv(swiss_glm,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ab_allSE.csv")







#Discriptive Plots

#all_slopes
env_slope<- ggplot(swiss_glm_filter,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
#env_slope <- env_slope  + theme(
#  axis.text.x = element_text(size=12, face="bold"),
#  axis.text.y = element_text(size=12,face="bold"),
#  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
#  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_slope <- env_slope + facet_wrap(.~Site)
env_slope

