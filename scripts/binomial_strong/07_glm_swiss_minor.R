##################################################################################
## Calculate glm slope and SE for every BF>0 snp
## Computationally intensive
## Split into 100k chucks to make managable
## Author Daniel Anstett
## 
## 
## Last Modified Dec 18, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions

#Calibrate timeseries abundance so that A is always the minor allele

minor_A <- function(df,baseline){
  swiss_minor_1 <- data.frame()
  for(i in 1:dim(df)[1]){
    
    start_swiss <- baseline[i,] %>% select(-chr_snp)
    A_col <- sum(start_swiss[,seq(1, ncol(start_swiss), by = 2)])
    B_col <- sum(start_swiss[,seq(2, ncol(start_swiss), by = 2)])   
    swiss_minor_1[i,1]<- df[i,1]
    
    for(j in 0:61){
      if(A_col >= B_col){
        #B is the minor
        swiss_minor_1[i,2+(2*j)] <- df[i,3+(2*j)]
        swiss_minor_1[i,3+(2*j)] <- df[i,2+(2*j)]
      } else{
        swiss_minor_1[i,2+(2*j)] <- df[i,2+(2*j)]
        swiss_minor_1[i,3+(2*j)] <- df[i,3+(2*j)]
      }
    }
    #print(i)
  }
  colnames(swiss_minor_1) <- colnames(df)
  return(swiss_minor_1)
}


#Generate frequency matrix for prop A 
abA <- function(snp_table,pop_ID) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-as.numeric(tmp$A) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  
  colnames(snp_prop_A)<- pop_ID[,1] #name each pop/time combination
  rownames(snp_prop_A)<- snp_prop_A$chr_snp
  snp1A_T <- as.data.frame(t(snp_prop_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
  colnames(snp1A_T) <- snp1A_T[1,]
  snp1A_T <- snp1A_T[-1,]
  colnames(snp1A_T)[1]<- "Site"
  colnames(snp1A_T)[2]<- "Year"
  return(snp1A_T)
}



abB <- function(snp_table,pop_ID) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-as.numeric(tmp$B) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  
  colnames(snp_prop_A)<- pop_ID[,1] #name each pop/time combination
  rownames(snp_prop_A)<- snp_prop_A$chr_snp
  snp1A_T <- as.data.frame(t(snp_prop_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
  colnames(snp1A_T) <- snp1A_T[1,]
  snp1A_T <- snp1A_T[-1,]
  colnames(snp1A_T)[1]<- "Site"
  colnames(snp1A_T)[2]<- "Year"
  return(snp1A_T)
}



# Melt glm Josee
#Get Slope and SE from GLM
slope_melt <- function(dfA,dfB) {
  freq_1.melted <- reshape2::melt(dfA, id.vars = c("Site", "Year"))
  freq_2.melted <- reshape2::melt(dfB, id.vars = c("Site", "Year"))
  colnames(freq_1.melted)[3] <- "snp_ID"
  colnames(freq_1.melted)[4] <- "abA"
  colnames(freq_2.melted)[3] <- "snp_ID"
  colnames(freq_2.melted)[4] <- "abB"
  ab.melt <- cbind(freq_1.melted,freq_2.melted$abB)
  colnames(ab.melt)[5] <- "abB"
  rm(freq_1.melted)
  rm(freq_2.melted)

#Extract slope and SE from binomial glm removing all instances where there is not enough data
  freq_1.slope <- group_by(ab.melt, Site, snp_ID) %>%
    arrange(Site, Year) %>% 
    summarize(Slope = ifelse(sum(!is.na(as.numeric(abA) / as.numeric(abB))) >= 2,
        glm(cbind(as.numeric(abA), as.numeric(abB)) ~ as.numeric(Year), family = binomial)$coefficients[2],
        ifelse(sum(!is.na(as.numeric(abB) / as.numeric(abA))) >= 2,
          glm(cbind(as.numeric(abA), as.numeric(abB)) ~ as.numeric(Year), family = binomial)$coefficients[2],
          NA)),
        SE = ifelse(sum(!is.na(as.numeric(abA) / as.numeric(abB))) >= 2,
                    summary(glm(cbind(as.numeric(abA),as.numeric(abB)) ~ as.numeric(Year), family = binomial))$coefficients[2,2],
                    ifelse(sum(!is.na(as.numeric(abB) / as.numeric(abA))) >= 2,
                    summary(glm(cbind(as.numeric(abA),as.numeric(abB)) ~ as.numeric(Year), family = binomial))$coefficients[2,2],
                         NA)))
  return(freq_1.slope)
}


#################################################################################################
## Manupulate entire SNP datatset for timeseries

#Import BF>0 baseline SNPs
env1_bf0 <- read_csv("data/env1_BF0.csv")
env2_bf0 <- read_csv("data/env2_BF0.csv")
env3_bf0 <- read_csv("data/env3_BF0.csv")
env4_bf0 <- read_csv("data/env4_BF0.csv")
env5_bf0 <- read_csv("data/env5_BF0.csv")
env6_bf0 <- read_csv("data/env6_BF0.csv")
env7_bf0 <- read_csv("data/env7_BF0.csv")
env8_bf0 <- read_csv("data/env8_BF0.csv")
env9_bf0 <- read_csv("data/env9_BF0.csv")

#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Make pop order to organize site/year headers
pop_order_2 <- data.frame()
pop_order_2 [1,1] <- "chr_shp"
pop_order_2 <- rbind(pop_order_2,pop_order)

#Import Baseline data
#import full snp table
pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")
snp_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                     header=F, sep=" ")
loci_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                      header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united_base <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_base <-cbind(loci_united_base,snp_base) #add snp lables to rows

#Filter full snp table to remove climate associated SNPs
snp_swiss <-loci_snp %>% filter (!chr_snp %in% as.character(env1_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env2_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env3_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env4_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env5_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env6_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env7_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env8_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env9_bf0$chr_snp))

#Filter baseline by snp_swiss (reduced timeseries)
base_swiss <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_swiss$chr_snp))

#Filter snp_swiss (reduced timeseries) by baseline
snp_swiss <- snp_swiss %>% filter (chr_snp %in% as.character(base_swiss$chr_snp))


tail(base_swiss)[,1]
tail(snp_swiss)[,1]

#Clear memory
rm(loci)
rm(loci_snp)
rm(snp)
rm(loci_united)
rm(env1_bf0)
rm(env2_bf0)
rm(env3_bf0)
rm(env4_bf0)
rm(env5_bf0)
rm(env6_bf0)
rm(env7_bf0)
rm(env8_bf0)
rm(env9_bf0)
rm(pop_order)



#################################################################################################

#Toy example to get minor allele per snp/pop/time

#swiss_1 <- snp_swiss[1:10,]
#base_1 <- base_swiss[1:10,]

#swiss_minor_1 <- data.frame()

#for(i in 1:10){
  
  #Get minror allele 
#  start_swiss <- base_1[i,] %>% select(-chr_snp)
#  A_col <- sum(start_swiss[,seq(1, ncol(start_swiss), by = 2)])
#  B_col <- sum(start_swiss[,seq(2, ncol(start_swiss), by = 2)])   
#  swiss_minor_1[i,1]<- swiss_1[i,1]
  
#  for(j in 0:61){
#    if(A_col >= B_col){
#      #B is the minor
#      swiss_minor_1[i,2+(2*j)] <- swiss_1[i,3+(2*j)]
#      swiss_minor_1[i,3+(2*j)] <- swiss_1[i,2+(2*j)]
#    } else{
#      swiss_minor_1[i,2+(2*j)] <- swiss_1[i,2+(2*j)]
#      swiss_minor_1[i,3+(2*j)] <- swiss_1[i,3+(2*j)]
#    }
#  }
#}

#swiss_abA <-  abA(swiss_1,pop_order_2)
#swiss_abB <-  abB(swiss_1,pop_order_2)


#################################################################################################
#GLM for Swiss


swiss_1 <- snp_swiss[1:100000,] #Split data set into 100k ch
base_1 <- base_swiss[1:100000,]
#swiss_1 <- snp_swiss[1:1000,] #Split data set into 100k ch
#base_1 <- base_swiss[1:1000,]
swiss_1 <- minor_A(swiss_1,base_1)
write_csv(swiss_1,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_1.csv")
#swiss_abA_1 <-  abA(swiss_1,pop_order_2)
#swiss_abB_1 <-  abB(swiss_1,pop_order_2)
#swiss_glm_1 <- slope_melt(swiss_abA_1,swiss_abB_1) #Run glm function
#write_csv(swiss_glm_1,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_1.csv")

rm(swiss_1)
rm(base_1)
#rm(swiss_abA_1)
#rm(swiss_abB_1)
#rm(swiss_glm_1)


swiss_2 <- snp_swiss[100001:200000,] #Split data set into 100k ch
base_2 <- base_swiss[100001:200000,] #Split data set into 100k ch
swiss_2 <- minor_A(swiss_2,base_2)
write_csv(swiss_2,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_2.csv")
#swiss_abA_2 <-  abA(swiss_2,pop_order_2)
#swiss_abB_2 <-  abB(swiss_2,pop_order_2)
#swiss_glm_2 <- slope_melt(swiss_abA_2,swiss_abB_2) #Run glm function
#write_csv(swiss_glm_2,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_2.csv")


rm(swiss_2)
rm(base_2)
#rm(swiss_abA_2)
#rm(swiss_abB_2)
#rm(swiss_glm_2)


swiss_3 <- snp_swiss[200001:300000,] #Split data set into 100k ch
base_3 <- base_swiss[200001:300000,] #Split data set into 100k ch
swiss_3 <- minor_A(swiss_3,base_3)
write_csv(swiss_3,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_3.csv")
#swiss_abA_3 <-  abA(swiss_3,pop_order_2)
#swiss_abB_3 <-  abB(swiss_3,pop_order_2)
#swiss_glm_3 <- slope_melt(swiss_abA_3,swiss_abB_3) #Run glm function
#write_csv(swiss_glm_3,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_3.csv")


rm(swiss_3)
rm(base_3)
#rm(swiss_abA_3)
#rm(swiss_abB_3)
#rm(swiss_glm_3)

swiss_4 <- snp_swiss[300001:400000,] #Split data set into 100k ch
base_4 <- base_swiss[300001:400000,] #Split data set into 100k ch
swiss_4 <- minor_A(swiss_4,base_4)
write_csv(swiss_4,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_4.csv")
#swiss_abA_4 <-  abA(swiss_4,pop_order_2)
#swiss_abB_4 <-  abB(swiss_4,pop_order_2)
#swiss_glm_4 <- slope_melt(swiss_abA_4,swiss_abB_4) #Run glm function
#write_csv(swiss_glm_4,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_4.csv")


rm(swiss_4)
rm(base_4)
#rm(swiss_abA_4)
#rm(swiss_abB_4)
#rm(swiss_glm_4)

swiss_5 <- snp_swiss[400001:500000,] #Split data set into 100k ch
base_5 <- base_swiss[400001:500000,] #Split data set into 100k ch
swiss_5 <- minor_A(swiss_5,base_5)
write_csv(swiss_5,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_5.csv")
#swiss_abA_5 <-  abA(swiss_5,pop_order_2)
#swiss_abB_5 <-  abB(swiss_5,pop_order_2)
#swiss_glm_5 <- slope_melt(swiss_abA_5,swiss_abB_5) #Run glm function
#write_csv(swiss_glm_5,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_5.csv")


rm(swiss_5)
rm(base_5)
#rm(swiss_abA_5)
#rm(swiss_abB_5)
#rm(swiss_glm_5)

swiss_6 <- snp_swiss[500001:600000,] #Split data set into 100k ch
base_6 <- base_swiss[500001:600000,] #Split data set into 100k ch
swiss_6 <- minor_A(swiss_6,base_6)
write_csv(swiss_6,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_6.csv")
#swiss_abA_6 <-  abA(swiss_6,pop_order_2)
#swiss_abB_6 <-  abB(swiss_6,pop_order_2)
#swiss_glm_6 <- slope_melt(swiss_abA_6,swiss_abB_6) #Run glm function
#write_csv(swiss_glm_6,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_6.csv")


rm(swiss_6)
rm(base_6)
#rm(swiss_abA_6)
#rm(swiss_abB_6)
#rm(swiss_glm_6)

swiss_7 <- snp_swiss[600001:700000,] #Split data set into 100k ch
base_7 <- base_swiss[600001:700000,] #Split data set into 100k ch
swiss_7 <- minor_A(swiss_7,base_7)
write_csv(swiss_7,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_7.csv")
#swiss_abA_7 <-  abA(swiss_7,pop_order_2)
#swiss_abB_7 <-  abB(swiss_7,pop_order_2)
#swiss_glm_7 <- slope_melt(swiss_abA_7,swiss_abB_7) #Run glm function
#write_csv(swiss_glm_7,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_7.csv")


rm(swiss_7)
rm(base_7)
#rm(swiss_abA_7)
#rm(swiss_abB_7)
#rm(swiss_glm_7)

swiss_8 <- snp_swiss[700001:800000,] #Split data set into 100k ch
base_8 <- base_swiss[700001:800000,] #Split data set into 100k ch
swiss_8 <- minor_A(swiss_8,base_8)
write_csv(swiss_8,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_8.csv")
#swiss_abA_8 <-  abA(swiss_8,pop_order_2)
#swiss_abB_8 <-  abB(swiss_8,pop_order_2)
#swiss_glm_8 <- slope_melt(swiss_abA_8,swiss_abB_8) #Run glm function
#write_csv(swiss_glm_8,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_8.csv")


rm(swiss_8)
rm(base_8)
#rm(swiss_abA_8)
#rm(swiss_abB_8)
#rm(swiss_glm_8)

swiss_9 <- snp_swiss[800001:900000,] #Split data set into 100k ch
base_9 <- base_swiss[800001:900000,] #Split data set into 100k ch
swiss_9 <- minor_A(swiss_9,base_9)
write_csv(swiss_9,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_9.csv")
#swiss_abA_9 <-  abA(swiss_9,pop_order_2)
#swiss_abB_9 <-  abB(swiss_9,pop_order_2)
#swiss_glm_9 <- slope_melt(swiss_abA_9,swiss_abB_9) #Run glm function
#write_csv(swiss_glm_9,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_9.csv")


rm(swiss_9)
rm(base_9)
#rm(swiss_abA_9)
#rm(swiss_abB_9)
#rm(swiss_glm_9)

swiss_10 <- snp_swiss[900001:1000000,] #Split data set into 100k ch
base_10 <- base_swiss[900001:1000000,] #Split data set into 100k ch
swiss_10 <- minor_A(swiss_10,base_10)
write_csv(swiss_10,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_10.csv")
#swiss_abA_10 <-  abA(swiss_10,pop_order_2)
#swiss_abB_10 <-  abB(swiss_10,pop_order_2)
#swiss_glm_10 <- slope_melt(swiss_abA_10,swiss_abB_10) #Run glm function
#write_csv(swiss_glm_10,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_10.csv")


rm(swiss_10)
rm(base_10)
#rm(swiss_abA_10)
#rm(swiss_abB_10)
#rm(swiss_glm_10)


swiss_11 <- snp_swiss[1000001:1100000,] #Split data set into 100k ch
base_11 <- base_swiss[1000001:1100000,] #Split data set into 100k ch
swiss_11 <- minor_A(swiss_11,base_11)
write_csv(swiss_11,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_11.csv")
#swiss_abA_11 <-  abA(swiss_11,pop_order_2)
#swiss_abB_11 <-  abB(swiss_11,pop_order_2)
#swiss_glm_11 <- slope_melt(swiss_abA_11,swiss_abB_11) #Run glm function
#write_csv(swiss_glm_11,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_11.csv")


rm(swiss_11)
rm(base_11)
#rm(swiss_abA_11)
#rm(swiss_abB_11)
#rm(swiss_glm_11)

swiss_12 <- snp_swiss[1100001:1200000,] #Split data set into 100k ch
base_12 <- base_swiss[1100001:1200000,] #Split data set into 100k ch
swiss_12 <- minor_A(swiss_12,base_12)
write_csv(swiss_12,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_12.csv")
#swiss_abA_12 <-  abA(swiss_12,pop_order_2)
#swiss_abB_12 <-  abB(swiss_12,pop_order_2)
#swiss_glm_12 <- slope_melt(swiss_abA_12,swiss_abB_12) #Run glm function
#write_csv(swiss_glm_12,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_12.csv")


rm(swiss_12)
rm(base_12)
#rm(swiss_abA_12)
#rm(swiss_abB_12)
#rm(swiss_glm_12)

swiss_13 <- snp_swiss[1200001:1300000,] #Split data set into 100k ch
base_13 <- base_swiss[1200001:1300000,] #Split data set into 100k ch
swiss_13 <- minor_A(swiss_13,base_13)
write_csv(swiss_13,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_13.csv")
#swiss_abA_13 <-  abA(swiss_13,pop_order_2)
#swiss_abB_13 <-  abB(swiss_13,pop_order_2)
#swiss_glm_13 <- slope_melt(swiss_abA_13,swiss_abB_13) #Run glm function
#write_csv(swiss_glm_13,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_13.csv")


rm(swiss_13)
rm(base_13)
#rm(swiss_abA_13)
#rm(swiss_abB_13)
#rm(swiss_glm_13)

swiss_14 <- snp_swiss[1300001:1400000,] #Split data set into 100k ch
base_14 <- base_swiss[1300001:1400000,] #Split data set into 100k ch
swiss_14 <- minor_A(swiss_14,base_14)
write_csv(swiss_14,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_14.csv")
#swiss_abA_14 <-  abA(swiss_14,pop_order_2)
#swiss_abB_14 <-  abB(swiss_14,pop_order_2)
#swiss_glm_14 <- slope_melt(swiss_abA_14,swiss_abB_14) #Run glm function
#write_csv(swiss_glm_14,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_14.csv")


rm(swiss_14)
rm(base_14)
#rm(swiss_abA_14)
#rm(swiss_abB_14)
#rm(swiss_glm_14)

swiss_15 <- snp_swiss[1400001:1500000,] #Split data set into 100k ch
base_15 <- base_swiss[1400001:1500000,] #Split data set into 100k ch
swiss_15 <- minor_A(swiss_15,base_15)
write_csv(swiss_15,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_15.csv")
#swiss_abA_15 <-  abA(swiss_15,pop_order_2)
#swiss_abB_15 <-  abB(swiss_15,pop_order_2)
#swiss_glm_15 <- slope_melt(swiss_abA_15,swiss_abB_15) #Run glm function
#write_csv(swiss_glm_15,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_15.csv")


rm(swiss_15)
rm(base_15)
#rm(swiss_abA_15)
#rm(swiss_abB_15)
#rm(swiss_glm_15)

swiss_16 <- snp_swiss[1500001:1600000,] #Split data set into 100k ch
base_16 <- base_swiss[1500001:1600000,] #Split data set into 100k ch
swiss_16 <- minor_A(swiss_16,base_16)
write_csv(swiss_16,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_16.csv")
#swiss_abA_16 <-  abA(swiss_16,pop_order_2)
#swiss_abB_16 <-  abB(swiss_16,pop_order_2)
#swiss_glm_16 <- slope_melt(swiss_abA_16,swiss_abB_16) #Run glm function
#write_csv(swiss_glm_16,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_16.csv")


rm(swiss_16)
rm(base_16)
#rm(swiss_abA_16)
#rm(swiss_abB_16)
#rm(swiss_glm_16)

swiss_17 <- snp_swiss[1600001:1607008,] #Split data set into 100k ch
base_17 <- base_swiss[1600001:1607008,] #Split data set into 100k ch
swiss_17 <- minor_A(swiss_17,base_17)
write_csv(swiss_17,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_17.csv")
#swiss_abA_17 <-  abA(swiss_17,pop_order_2)
#swiss_abB_17 <-  abB(swiss_17,pop_order_2)
#swiss_glm_17 <- slope_melt(swiss_abA_17,swiss_abB_17) #Run glm function
#write_csv(swiss_glm_17,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_minor_17.csv")


rm(swiss_17)
rm(base_17)
#rm(swiss_abA_17)
#rm(swiss_abB_17)
#rm(swiss_glm_17)











