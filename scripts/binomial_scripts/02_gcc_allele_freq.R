##################################################################################
## Generate SNP proportion matrix per population
## For peak window SNPS (WZA) with BF (BayPass) > 10 and all BF>30
## Done for all alleles positively associated with climate change
## Select if A or B is positively associated with climate change
## Author Daniel Anstett
## 
## For all 9 env variables
## Last Modified April 24, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
##setup functions

#Generate frequency matrix for prop A 
prop_A <- function(snp_table) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-tmp$A/(tmp$A + tmp$B) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  #snp_prop_A[is.na(snp_prop_A)] <- 0
  return(snp_prop_A)
}

#Generate frequency matrix for prop B
prop_B <- function(snp_table){
  #B Mat
  snp_prop_B<-snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="")
    snpB<-paste("V",i+1, sep="")
    P<-paste("P", pop_num, sep="")
    tmp<-snp_table %>% select(snpA,snpB)
    
    colnames(tmp)<-c("A", "B")
    
    snp_prop_B[,counter]<-tmp$B/(tmp$A + tmp$B)
    colnames (snp_prop_B)[counter]<-P
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  #snp_prop_B[is.na(snp_prop_B)] <- 0
  return(snp_prop_B)
}

###################################################################################
##set up frequency association table
#Positive assocation with climate
FAT_p <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-prop_A(snp_base) # call function
  snp_prop_B_in<-prop_B(snp_base) # call function
  
  snp_prop_A_in_time<-prop_A(snp_time)
  snp_prop_B_in_time<-prop_B(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
    env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
    colnames(env_pop)[3] <-"prop_A"
    colnames(env_pop)[4]<-"prop_B"
    
    lm.temp_A <- lm(env_pop$prop_A~env_pop[,2]) # save lm of cliamte predicting prop A
    lm.temp_B <- lm(env_pop$prop_B~env_pop[,2]) # save lm of cliamte predicting prop B
    
    #decides if A or B is positively associated 
    if (isTRUE(lm.temp_A$coefficents[2]>0) && isTRUE(lm.temp_B$coefficents[2]>0)){
      
      print(i)
      print(lm.temp_A) 
      print(lm.temp_B) 
      
    } else if(lm.temp_A$coefficients[2]>0){ 
      #print("A")
      tmp_in<-snp_prop_A_in_time %>% filter(chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (lm.temp_B$coefficients[2]>0){
      # print("B")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    }
  }
  return(freq.temp)
}


#Negative assocation with climate
FAT_n <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-prop_A(snp_base) # call function
  snp_prop_B_in<-prop_B(snp_base) # call function
  
  snp_prop_A_in_time<-prop_A(snp_time)
  snp_prop_B_in_time<-prop_B(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
    env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
    colnames(env_pop)[3] <-"prop_A"
    colnames(env_pop)[4]<-"prop_B"
    
    lm.temp_A <- lm(env_pop$prop_A~env_pop[,2]) # save lm of cliamte predicting prop A
    lm.temp_B <- lm(env_pop$prop_B~env_pop[,2]) # save lm of cliamte predicting prop B
    
    #decides if A or B is positively associated 
    if (isTRUE(lm.temp_A$coefficents[2]<0) && isTRUE(lm.temp_B$coefficents[2]<0)){
      
      print(i)
      print(lm.temp_A) 
      print(lm.temp_B) 
      
    } else if(lm.temp_A$coefficients[2]<0){ 
      #print("A")
      tmp_in<-snp_prop_A_in_time %>% filter( chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (lm.temp_B$coefficients[2]<0){
      # print("B")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
    }
  }
  return(freq.temp)
}
########################################################################################################
########################################################################################################

#Import Baseline Climate 
climate <- read_csv("data/climate_pop.csv")

#Import pop names (site_year names)
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

##SNP set abundances
#Timeseries SNP abundances
snp1_time <- read_csv("data/snp_set_time_env1.csv")
snp2_time <- read_csv("data/snp_set_time_env2.csv")
snp3_time <- read_csv("data/snp_set_time_env3.csv")
snp4_time <- read_csv("data/snp_set_time_env4.csv")
snp5_time <- read_csv("data/snp_set_time_env5.csv")
snp6_time <- read_csv("data/snp_set_time_env6.csv")
snp7_time <- read_csv("data/snp_set_time_env7.csv")
snp8_time <- read_csv("data/snp_set_time_env8.csv")
snp9_time <- read_csv("data/snp_set_time_env9.csv")

#Baseline SNP abundances
env1_base <- read_csv("data/snp_set_base_env1.csv")
env2_base <- read_csv("data/snp_set_base_env2.csv")
env3_base <- read_csv("data/snp_set_base_env3.csv")
env4_base <- read_csv("data/snp_set_base_env4.csv")
env5_base <- read_csv("data/snp_set_base_env5.csv")
env6_base <- read_csv("data/snp_set_base_env6.csv")
env7_base <- read_csv("data/snp_set_base_env7.csv")
env8_base <- read_csv("data/snp_set_base_env8.csv")
env9_base <- read_csv("data/snp_set_base_env9.csv")

###################################################################################

#Ensure baseline has same SNPs as timeseries and remove any not in the timeseries 
env1_base <- env1_base %>% filter (chr_snp %in% as.character(snp1_time$chr_snp))
env2_base <- env2_base %>% filter (chr_snp %in% as.character(snp2_time$chr_snp))
env3_base <- env3_base %>% filter (chr_snp %in% as.character(snp3_time$chr_snp))
env4_base <- env4_base %>% filter (chr_snp %in% as.character(snp4_time$chr_snp))
env5_base <- env5_base %>% filter (chr_snp %in% as.character(snp5_time$chr_snp))
env6_base <- env6_base %>% filter (chr_snp %in% as.character(snp6_time$chr_snp))
env7_base <- env7_base %>% filter (chr_snp %in% as.character(snp7_time$chr_snp))
env8_base <- env8_base %>% filter (chr_snp %in% as.character(snp8_time$chr_snp))
env9_base <- env9_base %>% filter (chr_snp %in% as.character(snp9_time$chr_snp))


########################################################################################################
########################################################################################################
##### Annual #####
# env 1 is MAT = Mean annual temperature (°C)
# env 2 is MAP = Mean annual precipitation (mm)
# env 3 is PAS = Precipitation as snow (mm) between August in previous year and July in current year
# env 4 is EXT = Extreme temperature over 30 years
# env 5 is CMD = Hargreaves climatic moisture deficit (mm)

##### Seasonal #####
# env 6 is Tave_wt = Winter mean temperature (°C)
# env 7 is Tave_sm = Summer mean temperature (°C)
# env 8 is PPT_wt = Winter precipitation (mm)
# env 9 is PPT_sm = Summer precipitation (mm)

## Make table with climate change associated SNPs for timeseries using baseline climate
freq_env1 <- FAT_p(env1_base,snp1_time,climate,"MAT")
freq_env2 <- FAT_n(env2_base,snp2_time,climate,"MAP")
freq_env3 <- FAT_n(env3_base,snp3_time,climate,"PAS")
freq_env4 <- FAT_p(env4_base,snp4_time,climate,"EXT")
freq_env5 <- FAT_p(env5_base,snp5_time,climate,"CMD")
freq_env6 <- FAT_p(env6_base,snp6_time,climate,"Tave_wt")
freq_env7 <- FAT_p(env7_base,snp7_time,climate,"Tave_sm")
freq_env8 <- FAT_n(env8_base,snp8_time,climate,"PPT_wt")
freq_env9 <- FAT_n(env9_base,snp9_time,climate,"PPT_sm")

#Make into dataframe from tibble
freq_env1 <- as.data.frame(freq_env1)
freq_env2 <- as.data.frame(freq_env2)
freq_env3 <- as.data.frame(freq_env3)
freq_env4 <- as.data.frame(freq_env4)
freq_env5 <- as.data.frame(freq_env5)
freq_env6 <- as.data.frame(freq_env6)
freq_env7 <- as.data.frame(freq_env7)
freq_env8 <- as.data.frame(freq_env8)
freq_env9 <- as.data.frame(freq_env9)

#Make first column a row name
rownames(freq_env1)<- as.vector(freq_env1$chr_snp)
rownames(freq_env2)<- as.vector(freq_env2$chr_snp)
rownames(freq_env3)<- as.vector(freq_env3$chr_snp)
rownames(freq_env4)<- as.vector(freq_env4$chr_snp)
rownames(freq_env5)<- as.vector(freq_env5$chr_snp)
rownames(freq_env6)<- as.vector(freq_env6$chr_snp)
rownames(freq_env7)<- as.vector(freq_env7$chr_snp)
rownames(freq_env8)<- as.vector(freq_env8$chr_snp)
rownames(freq_env9)<- as.vector(freq_env9$chr_snp)

#Remove 1st column
freq_env1 <- freq_env1 %>% select(-chr_snp)
freq_env2 <- freq_env2 %>% select(-chr_snp)
freq_env3 <- freq_env3 %>% select(-chr_snp)
freq_env4 <- freq_env4 %>% select(-chr_snp)
freq_env5 <- freq_env5 %>% select(-chr_snp)
freq_env6 <- freq_env6 %>% select(-chr_snp)
freq_env7 <- freq_env7 %>% select(-chr_snp)
freq_env8 <- freq_env8 %>% select(-chr_snp)
freq_env9 <- freq_env9 %>% select(-chr_snp)

# Add in row names to SNP frequency tables 
colnames(freq_env1)<- pop_order[,1] #name each pop/time combination
colnames(freq_env2)<- pop_order[,1] 
colnames(freq_env3)<- pop_order[,1] 
colnames(freq_env4)<- pop_order[,1] 
colnames(freq_env5)<- pop_order[,1] 
colnames(freq_env6)<- pop_order[,1] 
colnames(freq_env7)<- pop_order[,1] 
colnames(freq_env8)<- pop_order[,1] 
colnames(freq_env9)<- pop_order[,1] 

#Transpose and split up site_year
freq_env1_T <- as.data.frame(t(freq_env1)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env2_T <- as.data.frame(t(freq_env2)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env3_T <- as.data.frame(t(freq_env3)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env4_T <- as.data.frame(t(freq_env4)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env5_T <- as.data.frame(t(freq_env5)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env6_T <- as.data.frame(t(freq_env6)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env7_T <- as.data.frame(t(freq_env7)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env8_T <- as.data.frame(t(freq_env8)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env9_T <- as.data.frame(t(freq_env9)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))



# Write out climate change correlated timeseries frequency table 
write_csv(freq_env1_T, "data/freq_env1.csv")
write_csv(freq_env2_T, "data/freq_env2.csv")
write_csv(freq_env3_T, "data/freq_env3.csv")
write_csv(freq_env4_T, "data/freq_env4.csv")
write_csv(freq_env5_T, "data/freq_env5.csv")
write_csv(freq_env6_T, "data/freq_env6.csv")
write_csv(freq_env7_T, "data/freq_env7.csv")
write_csv(freq_env8_T, "data/freq_env8.csv")
write_csv(freq_env9_T, "data/freq_env9.csv")




  

