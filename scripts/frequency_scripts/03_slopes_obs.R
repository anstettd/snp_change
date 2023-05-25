##################################################################################
## Get pop specific SNP slopes over the timeseries for "observed" data set
## Done for snp set BF>10 in windows and BF>30 overall
## Does not include 1000 permuations of random sampling for BF<0
## Also done for unique slopes across all env (histPop)
## Author Daniel Anstett
## 
## 
## Last Modified April 24, 2023
###################################################################################
#Functions

# Run GLM

glm_obs <-function(df,env_input){
  freq_env1_slope<-data.frame()
  counter<-1
  for (i in 3:dim(df)[2]){
    for(j in 1:12){
      chr<-colnames(df)[i]
      popSNP <- df %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
      colnames(popSNP)[3]<-"snp_ID"
      
      if(all(is.na(popSNP$snp_ID))==FALSE & sum(is.na(popSNP$snp_ID)) !=(dim(popSNP)[1]-1)){
        rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
        freq_env1_slope[counter,1]<-unique(popSNP$Site)
        freq_env1_slope[counter,2]<-chr
        freq_env1_slope[counter,3]<-rSNP$coefficients[2] * 2
        sum_snp <- summary(rSNP)
        freq_env1_slope[counter,4]<-sum_snp$coefficients[2,2]
        freq_env1_slope[counter,5]<-env_input
        freq_env1_slope[counter,6]<-"obs"
        counter<-counter+1
      } else {
        freq_env1_slope[counter,1]<-unique(popSNP$Site)
        freq_env1_slope[counter,2]<-chr
        freq_env1_slope[counter,3]<-NA
        freq_env1_slope[counter,4]<-NA
        freq_env1_slope[counter,5]<-env_input
        freq_env1_slope[counter,6]<-"obs"
        counter<-counter+1
      }
      
    }
  }
  colnames(freq_env1_slope)<-c("Site","snp_ID","Slope","SE","ENV","Type")
  return(freq_env1_slope)
}



###################################################################################
#Import libraries
library(tidyverse)
library(boot)

#Import timeseries frequencies
freq_env1 <- read_csv("data/freq_env1.csv")
freq_env2 <- read_csv("data/freq_env2.csv")
freq_env3 <- read_csv("data/freq_env3.csv")
freq_env4 <- read_csv("data/freq_env4.csv")
freq_env5 <- read_csv("data/freq_env5.csv")
freq_env6 <- read_csv("data/freq_env6.csv")
freq_env7 <- read_csv("data/freq_env7.csv")
freq_env8 <- read_csv("data/freq_env8.csv")
freq_env9 <- read_csv("data/freq_env9.csv")


###################################################################################
#Single Case
#freq_env1_pop1 <- freq_env1 %>% filter(Site==1)
#rSNP_pop1 <- glm(CE10_chr1_865732 ~ Year, family = binomial, data = freq_env1_pop1)
#test_pop1 <- as.data.frame(summary(rSNP_pop1)$coefficients)
#test_pop1$Estimate[2]
#rSNP_pop1$coefficients[2]
#summary(rSNP_pop1)
#rSNP_pop1 <- glm(CE10_chr1_9110324 ~ Year, family = binomial, data = freq_env1_pop1)
#rSNP_pop1$coefficients[2]
#summary(rSNP_pop1)

###################################################################################

## Single Loop

# freq_env1
freq_env1_slope<-data.frame()
counter<-1
for (i in 3:dim(freq_env1)[2]){
  for(j in 1:12){
    chr<-colnames(freq_env1)[i]
    popSNP <- freq_env1 %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE & sum(is.na(popSNP$snp_ID)) !=(dim(popSNP)[1]-1)){
      rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
      freq_env1_slope[counter,1]<-unique(popSNP$Site)
      freq_env1_slope[counter,2]<-chr
      freq_env1_slope[counter,3]<-rSNP$coefficients[2] * 2
      sum_snp <- summary(rSNP)
      freq_env1_slope[counter,4]<-sum_snp$coefficients[2,2]
      freq_env1_slope[counter,5]<-"MAT"
      freq_env1_slope[counter,6]<-"obs"
      counter<-counter+1
    } else {
      freq_env1_slope[counter,1]<-unique(popSNP$Site)
      freq_env1_slope[counter,2]<-chr
      freq_env1_slope[counter,3]<-NA
      freq_env1_slope[counter,4]<-NA
      freq_env1_slope[counter,5]<-"MAT"
      freq_env1_slope[counter,6]<-"obs"
      
      counter<-counter+1
    }
    
  }
}
colnames(freq_env1_slope)<-c("Site","snp_ID","Slope","SE","ENV","Type")

plot(freq_env1_slope$Slope,freq_env1_slope$SE)

###################################################################################
###################################################################################
#Calc glm per env

slope_env1 <- glm_obs(freq_env1,"MAT")
slope_env2 <- glm_obs(freq_env2,"MAP")
slope_env3 <- glm_obs(freq_env3,"PAS")
slope_env4 <- glm_obs(freq_env4,"EXT")
slope_env5 <- glm_obs(freq_env5,"CMD")
slope_env6 <- glm_obs(freq_env6,"Tave_wt")
slope_env7 <- glm_obs(freq_env7,"Tave_sm")
slope_env8 <- glm_obs(freq_env8,"PPT_wt")
slope_env9 <- glm_obs(freq_env9,"PPT_sm")
slope_env_all <- rbind(slope_env1,
                       slope_env2,
                       slope_env3,
                       slope_env4,
                       slope_env5,
                       slope_env6,
                       slope_env7,
                       slope_env8,
                       slope_env9)

#Get unique SNP for merged env

slope_env_all_unique <- data.frame()
for (i in 1:dim(slope_env_all)){
  slope.temp <- slope_env_all %>% filter(Site==i)
  slope.unique <- slope.temp [!duplicated(slope.temp [ , "snp_ID"]), ]
  slope_env_all_unique <- rbind(slope_env_all_unique,slope.unique)
}


###################################################################################
#Export
write_csv(slope_env1, "data/slope_obs_env1.csv")
write_csv(slope_env2, "data/slope_obs_env2.csv")
write_csv(slope_env3, "data/slope_obs_env3.csv")
write_csv(slope_env4, "data/slope_obs_env4.csv")
write_csv(slope_env5, "data/slope_obs_env5.csv")
write_csv(slope_env6, "data/slope_obs_env6.csv")
write_csv(slope_env7, "data/slope_obs_env7.csv")
write_csv(slope_env8, "data/slope_obs_env8.csv")
write_csv(slope_env9, "data/slope_obs_env9.csv")
write_csv(slope_env_all, "data/slope_obs_all.csv")
write_csv(slope_env_all_unique, "data/slope_obs_all_unique.csv")




