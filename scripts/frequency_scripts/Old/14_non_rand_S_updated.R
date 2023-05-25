##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Updated to count only difference between obs and upper CI rather than whole bin
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
###################################################################################
#Import libraries
library(tidyverse)


###################################################################################
#Functions

ns <- function(ci,env_name){
  #Filter for positive slope bin env
  env_pos <- ci %>% filter(S>0)
  env_obs <- env_pos %>% select(p1:p12) 
  env_rand <-env_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
  instruction_env_P <- env_obs - env_rand 
  instruction_env_P <- cbind (instruction_env_P,env_pos$S)
  colnames(instruction_env_P)[13] <- "S"
  
  #Get non-random slopes
  #env
  slope.cumul.env <-  as.data.frame(c(1,2,3,4,5,6,7,8,9,10,11,12))
  colnames(slope.cumul.env) <- "Site"
  
  for (j in 1:12){
    instruction_env <- instruction_env_P %>% select(eval(paste("p",j, sep="")),S)
    instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
    Evol <- instruction_env[,1]*instruction_env[,2]
    slope.cumul.env[j,2] <- sum(Evol)
  }
  colnames(slope.cumul.env) <- c("Site",env_name)
  return(slope.cumul.env )
}    
    

ns_all <- function(ci,env_name){
  #Filter for positive slope bin env
  env_pos <- ci
  env_obs <- env_pos %>% select(p1:p12) 
  env_rand <-env_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
  instruction_env_P <- env_obs - env_rand 
  instruction_env_P <- cbind (instruction_env_P,env_pos$S)
  colnames(instruction_env_P)[13] <- "S"
  
  #Get non-random slopes
  #env
  slope.cumul.env <-  as.data.frame(c(1,2,3,4,5,6,7,8,9,10,11,12))
  colnames(slope.cumul.env) <- "Site"
  
  for (j in 1:12){
    instruction_env <- instruction_env_P %>% select(eval(paste("p",j, sep="")),S)
    instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
    Evol <- instruction_env[,1]*instruction_env[,2]
    slope.cumul.env[j,2] <- sum(abs(Evol))
  }
  colnames(slope.cumul.env) <- c("Site",env_name)
  return(slope.cumul.env )
}    



###################################################################################
#Import observed slopes
env1_slope <- read_csv("data/slope_obs_env1.csv")
env2_slope <- read_csv("data/slope_obs_env2.csv")
env3_slope <- read_csv("data/slope_obs_env3.csv")
env4_slope <- read_csv("data/slope_obs_env4.csv")
env5_slope <- read_csv("data/slope_obs_env5.csv")
env6_slope <- read_csv("data/slope_obs_env6.csv")
env7_slope <- read_csv("data/slope_obs_env7.csv")
env8_slope <- read_csv("data/slope_obs_env8.csv")
env9_slope <- read_csv("data/slope_obs_env9.csv")

#Import files that give slopes
env_obs_ci <- read_csv("data/obs_ci_env.csv")
env_obs_ci[is.na(env_obs_ci)] <- 0
env1_obs_ci <- env_obs_ci %>% filter(env == "A MAT")
env2_obs_ci <- env_obs_ci %>% filter(env == "B MAP")
env3_obs_ci <- env_obs_ci %>% filter(env == "C PAS")
env4_obs_ci <- env_obs_ci %>% filter(env == "D EXT")
env5_obs_ci <- env_obs_ci %>% filter(env == "E CMD")
env6_obs_ci <- env_obs_ci %>% filter(env == "F Tave_wt")
env7_obs_ci <- env_obs_ci %>% filter(env == "G Tave_sm")
env8_obs_ci <- env_obs_ci %>% filter(env == "H PPT_wt")
env9_obs_ci <- env_obs_ci %>% filter(env == "I PPT_sm")



############################################################################################################
#Example out of function

#Filter for positive slope bin env1
env1_pos <- env1_obs_ci %>% filter(S>0)
env1_obs <- env1_pos %>% select(p1:p12) 
env1_rand <-env1_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
instruction_env1 <- env1_obs - env1_rand 
instruction_env1 <- cbind (instruction_env1,env1_pos$S)
colnames(instruction_env1)[13] <- "S"

#Get non-random slopes
#env1
slope.cumul.env <-  as.data.frame(c(1,2,3,4,5,6,7,8,9,10,11,12))
colnames(slope.cumul.env) <- "Site"

for (j in 1:12){
    instruction_env <- instruction_env1 %>% select(eval(paste("p",j, sep="")),S)
    instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
    Evol <- instruction_env[,1]*instruction_env[,2]
  slope.cumul.env[j,2] <- sum(Evol)
}
colnames(slope.cumul.env) <- c("Site","env1")
  
  

#Filter for all slope bin env1
env1_pos <- env1_obs_ci
env1_obs <- env1_pos %>% select(p1:p12) 
env1_rand <-env1_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
instruction_env1 <- env1_obs - env1_rand 
instruction_env1 <- cbind (instruction_env1,env1_pos$S)
colnames(instruction_env1)[13] <- "S"

#Get non-random slopes
#env1
slope.cumul.env_all <-  as.data.frame(c(1,2,3,4,5,6,7,8,9,10,11,12))
colnames(slope.cumul.env_all) <- "Site"

for (j in 1:12){
  instruction_env <- instruction_env1 %>% select(eval(paste("p",j, sep="")),S)
  instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
  Evol <- instruction_env[,1]*instruction_env[,2]
  slope.cumul.env_all[j,2] <- sum(abs(Evol))
}
colnames(slope.cumul.env_all) <- c("Site","env1")


  
  
############################################################################################################
#Run ns function to get cumulative score
slope_pos_env1 <-ns(env1_obs_ci,"env1")
slope_pos_env2 <-ns(env2_obs_ci,"env2") %>% select(env2)
slope_pos_env3 <-ns(env3_obs_ci,"env3") %>% select(env3)
slope_pos_env4 <-ns(env4_obs_ci,"env4") %>% select(env4)
slope_pos_env5 <-ns(env5_obs_ci,"env5") %>% select(env5)
slope_pos_env6 <-ns(env6_obs_ci,"env6") %>% select(env6)
slope_pos_env7 <-ns(env7_obs_ci,"env7") %>% select(env7)
slope_pos_env8 <-ns(env8_obs_ci,"env8") %>% select(env8)
slope_pos_env9 <-ns(env9_obs_ci,"env9") %>% select(env9)

slope_all_env1 <-ns_all(env1_obs_ci,"env1")
slope_all_env2 <-ns_all(env2_obs_ci,"env2") %>% select(env2)
slope_all_env3 <-ns_all(env3_obs_ci,"env3") %>% select(env3)
slope_all_env4 <-ns_all(env4_obs_ci,"env4") %>% select(env4)
slope_all_env5 <-ns_all(env5_obs_ci,"env5") %>% select(env5)
slope_all_env6 <-ns_all(env6_obs_ci,"env6") %>% select(env6)
slope_all_env7 <-ns_all(env7_obs_ci,"env7") %>% select(env7)
slope_all_env8 <-ns_all(env8_obs_ci,"env8") %>% select(env8)
slope_all_env9 <-ns_all(env9_obs_ci,"env9") %>% select(env9)


#Bind data frames
cumul_env_pos <- cbind(slope_pos_env1,
                  slope_pos_env2,
                  slope_pos_env3,
                  slope_pos_env4,
                  slope_pos_env5,
                  slope_pos_env6,
                  slope_pos_env7,
                  slope_pos_env8,
                  slope_pos_env9) 

cumul_env_all <- cbind(slope_all_env1,
                   slope_all_env2,
                   slope_all_env3,
                   slope_all_env4,
                   slope_all_env5,
                   slope_all_env6,
                   slope_all_env7,
                   slope_all_env8,
                   slope_all_env9) 

#Add up env
cumul_pos <- cumul_env_pos  %>% rowwise(Site) %>% mutate(cumul_pos = sum(c(env1,
                                                               env2,
                                                               env3,
                                                               env4,
                                                               env5,
                                                               env6,
                                                               env7,
                                                               env8,
                                                               env9)))


cumul_all <- cumul_env_all  %>% rowwise(Site) %>% mutate(cumul_all = sum(c(env1,
                                                                       env2,
                                                                       env3,
                                                                       env4,
                                                                       env5,
                                                                       env6,
                                                                       env7,
                                                                       env8,
                                                                       env9)))


#Merge cumul
cumul_pos_noENV <- cumul_pos %>% select(Site,cumul_pos)
cumul_all_noENV <- cumul_all %>% select(Site,cumul_all)

cumul_grand <- left_join(cumul_pos_noENV,cumul_all_noENV,by="Site")


##Integrate with timeseries
timeseries <- read_csv("data/offset_pop_timeseries_beagle.csv")
time_cumul <-left_join(timeseries,cumul_grand,by=c("Paper_ID"="Site"))
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
write_csv(cumul_pos,"data/cumul_pos.csv")
write_csv(cumul_all,"data/cumul_all.csv")
write_csv(time_cumul,"data/time_cumul_beagle.csv")










############################################################################################################
#Make unique snp_list across env

#Diagnostics
#snp_list_p1 <- snp_list %>% filter(Site==1)
#Not unique SNP ID
#unique(snp_list_p1$snp_ID[duplicated(snp_list_p1$snp_ID)]) 
#duplicated(snp_list_p1$snp_ID[duplicated(snp_list_p1$snp_ID)]) 


#Remove not unique SNPs
snp_list_filter <- data.frame()
for(i in 1:12){
  snp_list_p <- snp_list %>% filter(Site==i)
  snp_list_p_filtered <- snp_list_p %>% filter(duplicated(snp_ID) == FALSE)
  snp_list_filter <- rbind(snp_list_filter,snp_list_p_filtered)
}




#write_csv(snp_list,"Genomics_scripts/Data/snp_list.csv")
write_csv(snp_list_filter,"data/snp_list_unique.csv")




############################################################################################################


