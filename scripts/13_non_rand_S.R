##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
###################################################################################
#Import libraries
library(tidyverse)


###################################################################################
#Functions

ns <- function(ci,slope_df){
  #Filter for positive slope bin env
  env_pos <- ci %>% filter(S>0)
  env_obs <- env_pos %>% select(p1:p12) 
  env_rand <-env_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
  instruction_env_P <- env_obs - env_rand 
  instruction_env_P <- cbind (instruction_env_P,env_pos$S)
  colnames(instruction_env_P)[13] <- "S"
  
  #Get non-random slopes
  #env
  slope.out.env <- data.frame()
  for (j in 1:12){
    instruction_env <- instruction_env_P %>% select(eval(paste("p",j, sep="")),S)
    instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
    env_slope_pop <- slope_df %>% filter(Site==j)
    
    if(dim(instruction_env)[1]!=0){
      instruction_env <- instruction_env %>%  mutate(low_S=S-0.1,high_S=S+0.1)
      
      for (i in 1:dim(instruction_env)[1]){
        data.temp <- instruction_env[i,]
        slope.temp <- env_slope_pop %>% filter(Slope >= data.temp$low_S & Slope < data.temp$high_S)
        slope.out.env <- rbind(slope.out.env,slope.temp)
      }
    }
  }
  return(slope.out.env)
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
slope.out.env1 <- data.frame()
for (j in 1:12){
  instruction_env <- instruction_env1 %>% select(eval(paste("p",j, sep="")),S)
  instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
  env_slope_pop <- env1_slope %>% filter(Site==j)
  
  if(dim(instruction_env)[1]!=0){
  instruction_env <- instruction_env %>%  mutate(low_S=S-0.1,high_S=S+0.1)
  
  for (i in 1:dim(instruction_env)[1]){
    data.temp <- instruction_env[i,]
    slope.temp <- env_slope_pop %>% filter(Slope >= data.temp$low_S & Slope < data.temp$high_S) 
    slope.out.env1 <- rbind(slope.out.env1,slope.temp)
  }
}
}

############################################################################################################
#Run ns function
slope_out_env1 <-ns(env1_obs_ci,env1_slope)
slope_out_env2 <-ns(env2_obs_ci,env2_slope)
slope_out_env3 <-ns(env3_obs_ci,env3_slope)
slope_out_env4 <-ns(env4_obs_ci,env4_slope)
slope_out_env5 <-ns(env5_obs_ci,env5_slope)
slope_out_env6 <-ns(env6_obs_ci,env6_slope)
slope_out_env7 <-ns(env7_obs_ci,env7_slope)
slope_out_env8 <-ns(env8_obs_ci,env8_slope)
slope_out_env9 <-ns(env9_obs_ci,env9_slope)

snp_list <- rbind(slope_out_env1,
                  slope_out_env2,
                  slope_out_env3,
                  slope_out_env4,
                  slope_out_env5,
                  slope_out_env6,
                  slope_out_env7,
                  slope_out_env8,
                  slope_out_env9)

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


