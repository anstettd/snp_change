##################################################################################
## Get pop specific SNP slopes over the timeseries for "observed" data set
## Done for snp set BF>10 in windows and BF>30 overall
## Does not include 1000 permuations of random sampling for BF<0
## Author Daniel Anstett
## 
## 
## Last Modified April 24, 2023
###################################################################################


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
freq_env1_pop1 <- freq_env1 %>% filter(Site==1)
rSNP_pop1 <- glm(CE10_chr1_865732 ~ Year, family = binomial, data = freq_env1_pop1)
test_pop1 <- as.data.frame(summary(rSNP_pop1)$coefficients)
test_pop1$Estimate[2]
rSNP_pop1$coefficients[2]

rSNP_pop1 <- glm(CE10_chr1_9110324 ~ Year, family = binomial, data = freq_env1_pop1)
rSNP_pop1$coefficients[2]



###################################################################################

## Generate SNP slope tables per site

# freq_env1
freq_env1_slope<-data.frame()
counter<-1
for (i in 3:dim(freq_env1)[2]){
  for(j in 1:12){
    chr<-colnames(freq_env1)[i]
    popSNP <- freq_env1 %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
      freq_env1_slope[counter,1]<-unique(popSNP$Site)
      freq_env1_slope[counter,2]<-chr
      freq_env1_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else {
      freq_env1_slope[counter,1]<-unique(popSNP$Site)
      freq_env1_slope[counter,2]<-chr
      freq_env1_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(freq_env1_slope)<-c("Site","snp_ID","Slope")

###################################################################################
#Export
write_csv(freq_env1_slope, "data/slope_obs_env1.csv")













# freq_MAP
freq_MAP_slope<-data.frame()
counter<-1
for (i in 3:dim(freq_MAP)[2]){
  for(j in 1:12){
    chr<-colnames(freq_MAP)[i]
    popSNP <- freq_MAP %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
      freq_MAP_slope[counter,1]<-unique(popSNP$Site)
      freq_MAP_slope[counter,2]<-chr
      freq_MAP_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else {
      freq_MAP_slope[counter,1]<-unique(popSNP$Site)
      freq_MAP_slope[counter,2]<-chr
      freq_MAP_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(freq_MAP_slope)<-c("Site","snp_ID","Slope")


# freq_CMD
freq_CMD_slope<-data.frame()
counter<-1
for (i in 3:dim(freq_CMD)[2]){
  for(j in 1:12){
    chr<-colnames(freq_CMD)[i]
    popSNP <- freq_CMD %>% filter(Site==j) %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE){
      rSNP <- glm(snp_ID ~ Year, family = binomial, data = popSNP)
      freq_CMD_slope[counter,1]<-unique(popSNP$Site)
      freq_CMD_slope[counter,2]<-chr
      freq_CMD_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else {
      freq_CMD_slope[counter,1]<-unique(popSNP$Site)
      freq_CMD_slope[counter,2]<-chr
      freq_CMD_slope[counter,3]<-NA
      counter<-counter+1
    }
    
  }
}
colnames(freq_CMD_slope)<-c("Site","snp_ID","Slope")



