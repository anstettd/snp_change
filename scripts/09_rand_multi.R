##################################################################################
## Generate stratified distribution of random non-climate associated slopes
## Author Daniel Anstett
## 
## 
## Last Modified April 15, 2022
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions for generating stratigied random distribution
###################################################################################
#Large
#pop_snp = p1A
#freq_count = freq_count_env1_IDX
large <- function(pop_snp,freq_count){
  rand_pop <- data.frame()
  bin_fraction<-1/length(freq_count)
  bin_step<-seq(0,1,bin_fraction)
  
  for (i in 1: length(freq_count)){
    if(i==length(freq_count)){
      p1_n<-pop_snp %>% filter(snpA>=bin_step[i] & snpA <= bin_step[i+1])
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i],replace = FALSE)
      
    }else{
      p1_n<-pop_snp %>% filter(snpA>=bin_step[i] & snpA < bin_step[i+1])
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i],replace = FALSE)
      
    }
   
    for (j in 1:length(list_env_p1_n)){
      rand_pop <- rbind(rand_pop ,p1_n %>% 
                          filter(chr_snp==as.character(p1_n$chr_snp[list_env_p1_n[j]])))
    }
    
  }

  return(rand_pop)
}

##Freq Count functions

#Diagnostic tool
#Used to diagnose what bin size to assign to each population. 
#Use function on data to generate bin table (freq_count_env1).
#Diagnose bin size per pop by taking into account the number and placement of zeros
freq_bins <- function(basetime){
  freq_count_calc<-data.frame()
  bin_size<-10
  for (i in 1:12){
    
#    if(i==1 || i==7){
#      bin_size<-10
#    }else if (i==3 || i==4 || i==6 || i==11){
#      bin_size<-4
#    }else if (i==8){
#      bin_size<-2
#    }else if(i==2 || i==5 || i==9 || i==10 || i==12){
#      bin_size<-5
#    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}

#env_bin5
#Increase binwidtch at 0.2 increments throughout
freq_bin5 <- function(basetime){
  freq_count_calc<-data.frame()
  bin_size<-5
  for (i in 1:12){
    
    #    if(i==1 || i==7){
    #      bin_size<-10
    #    }else if (i==3 || i==4 || i==6 || i==11){
    #      bin_size<-4
    #    }else if (i==8){
    #      bin_size<-2
    #    }else if(i==2 || i==5 || i==9 || i==10 || i==12){
    #      bin_size<-5
    #    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}


#env1
freq_bins_env1 <- function(basetime){
  freq_count_calc<-data.frame()
  for (i in 1:12){
    
    if(i==4){
      bin_size<-4
    }else{
      bin_size<-5
    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)

    #counting the number of SNPs that fall within a given frequency bin for a given Pop   
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
    if (j==bin_size){
      freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
    }else{
      freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}

#env6
freq_bins_env6 <- function(basetime){
  freq_count_calc<-data.frame()
  for (i in 1:12){
    
    if(i==8 | i==11){
      bin_size<-4
    }else{
      bin_size<-5
    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    #counting the number of SNPs that fall within a given frequency bin for a given Pop   
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}


#env7
freq_bins_env7 <- function(basetime){
  freq_count_calc<-data.frame()
  for (i in 1:12){
    
    if(i==2 | i==12){
      bin_size<-4
    }else{
      bin_size<-5
    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    #counting the number of SNPs that fall within a given frequency bin for a given Pop   
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}


#env8
freq_bins_env8 <- function(basetime){
  freq_count_calc<-data.frame()
  for (i in 1:12){
    
    if(i==2){
      bin_size<-4
    }else{
      bin_size<-5
    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    #counting the number of SNPs that fall within a given frequency bin for a given Pop   
    for(j in 1:bin_size){
      test_ENV <- basetime %>% filter(Site==i) %>% select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}


#Generate frequency matrix for prop A 
prop_A <- function(snp_table,pop_ID) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-as.numeric(tmp$A/(tmp$A + tmp$B)) #calc proportion for A, assign to output
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

glm_slopes<-function(snp_popX){
  
  snp_popX_slope<-data.frame()
  counter<-1
  
  for (i in 3:dim(snp_popX)[2]){
    chr<-colnames(snp_popX)[i]
    popSNP <- snp_popX %>% select(Site,Year,all_of(chr))
    colnames(popSNP)[3]<-"snp_ID"
    
    if(all(is.na(popSNP$snp_ID))==FALSE && length(unique(as.numeric(popSNP$snp_ID)))>1 
       && sum(!is.na( as.numeric((popSNP$snp_ID))))>1){
      popSNP <- na.omit(popSNP)
      rSNP <- glm(as.numeric(snp_ID) ~ Year, family = binomial, data = popSNP)
      snp_popX_slope[counter,1]<-unique(popSNP$Site)
      snp_popX_slope[counter,2]<-chr
      snp_popX_slope[counter,3]<-rSNP$coefficients[2]
      counter<-counter+1
    } else if((all(is.na(popSNP$snp_ID))==FALSE && length(unique(as.numeric(popSNP$snp_ID)))==1)){
      snp_popX_slope[counter,1]<-unique(popSNP$Site)
      snp_popX_slope[counter,2]<-chr
      snp_popX_slope[counter,3]<-0
      counter<-counter+1
    }else {
      snp_popX_slope[counter,1]<-unique(popSNP$Site)
      snp_popX_slope[counter,2]<-chr
      snp_popX_slope[counter,3]<-NA
      counter<-counter+1
    }
  }
    colnames(snp_popX_slope)<-c("Site","snp_ID","Slope")
 return(snp_popX_slope) 
}


###################################################################################


## Manupulate entire SNP datatset for timeseries

#Import timeseries counts for A and B
#snp1_time <- read_csv("Genomics_scripts/Data/snp1_filter.csv")
#snp2_time <- read_csv("Genomics_scripts/Data/snp2_filter.csv")
#snp5_time <- read_csv("Genomics_scripts/Data/snp5_filter.csv")

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

#Add Chromsome ID to pop order
pop_order_2 <- data.frame()
pop_order_2 [1,1] <- "chr_shp"
pop_order_2 <- rbind(pop_order_2,pop_order)

#Import pop_order table with regional and V information
pop_order_V <- read_csv("data/pop_order_V.csv")

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


# Import non-climate associated slopes with low SE
swiss_glm <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_filter.csv")


###################################################################################
#Import transformed timeseries frequencies
basetime_env1 <- read_csv("data/env1_low.csv")
basetime_env2 <- read_csv("data/env2_low.csv")
basetime_env3 <- read_csv("data/env3_low.csv")
basetime_env4 <- read_csv("data/env4_low.csv")
basetime_env5 <- read_csv("data/env5_low.csv")
basetime_env6 <- read_csv("data/env6_low.csv")
basetime_env7 <- read_csv("data/env7_low.csv")
basetime_env8 <- read_csv("data/env8_low.csv")
basetime_env9 <- read_csv("data/env9_low.csv")
###################################################################################

#Diagnose bin size per pop per env variable
#freq_count_env1 <- freq_bins(basetime_env1)
#freq_count_env1
#freq_count_env2 <- freq_bins(basetime_env2)
#freq_count_env2
#freq_count_env3 <- freq_bins(basetime_env3)
#freq_count_env3
#freq_count_env4 <- freq_bins(basetime_env4)
#freq_count_env4
#freq_count_env5 <- freq_bins(basetime_env5)
#freq_count_env5
#freq_count_env6 <- freq_bins(basetime_env6)
#freq_count_env6
#freq_count_env7 <- freq_bins(basetime_env7)
#freq_count_env7
#freq_count_env8 <- freq_bins(basetime_env8)
#freq_count_env8
#freq_count_env9 <- freq_bins(basetime_env9)
#freq_count_env9



#Make frequency count table
freq_count_env1 <- freq_bins_env1(basetime_env1)
freq_count_env2 <- freq_bin5(basetime_env2)
freq_count_env3 <- freq_bin5(basetime_env3)
freq_count_env4 <- freq_bin5(basetime_env4)
freq_count_env5 <- freq_bin5(basetime_env5)
freq_count_env6 <- freq_bins_env6(basetime_env6)
freq_count_env7 <- freq_bins_env7(basetime_env7)
freq_count_env8 <- freq_bins_env8(basetime_env8)
freq_count_env9 <- freq_bin5(basetime_env9)

#Check to ensure no zeros are present
#freq_count_env1
#freq_count_env2
#freq_count_env3
#freq_count_env4
#freq_count_env5
#freq_count_env6
#freq_count_env7
#freq_count_env8
#freq_count_env9


######################################################################################################
#Break up freq count into individual vectors
######################################################################################################
freq_count_env1_1 <- as.vector(freq_count_env1[,1])
freq_count_env1_2 <- as.vector(freq_count_env1[,2])
freq_count_env1_3 <- as.vector(freq_count_env1[,3])
freq_count_env1_4 <- as.vector(freq_count_env1[,4])
freq_count_env1_5 <- as.vector(freq_count_env1[,5])
freq_count_env1_6 <- as.vector(freq_count_env1[,6])
freq_count_env1_7 <- as.vector(freq_count_env1[,7])
freq_count_env1_8 <- as.vector(freq_count_env1[,8])
freq_count_env1_9 <- as.vector(freq_count_env1[,9])
freq_count_env1_10 <- as.vector(freq_count_env1[,10])
freq_count_env1_11 <- as.vector(freq_count_env1[,11])
freq_count_env1_12 <- as.vector(freq_count_env1[,12])


freq_count_env1_1<-freq_count_env1_1[!is.na(freq_count_env1_1)]
freq_count_env1_2<-freq_count_env1_2[!is.na(freq_count_env1_2)]
freq_count_env1_3<-freq_count_env1_3[!is.na(freq_count_env1_3)]
freq_count_env1_4<-freq_count_env1_4[!is.na(freq_count_env1_4)]
freq_count_env1_5<-freq_count_env1_5[!is.na(freq_count_env1_5)]
freq_count_env1_6<-freq_count_env1_6[!is.na(freq_count_env1_6)]
freq_count_env1_7<-freq_count_env1_7[!is.na(freq_count_env1_7)]
freq_count_env1_8<-freq_count_env1_8[!is.na(freq_count_env1_8)]
freq_count_env1_9<-freq_count_env1_9[!is.na(freq_count_env1_9)]
freq_count_env1_10<-freq_count_env1_10[!is.na(freq_count_env1_10)]
freq_count_env1_11<-freq_count_env1_11[!is.na(freq_count_env1_11)]
freq_count_env1_12<-freq_count_env1_12[!is.na(freq_count_env1_12)]

freq_count_env2_1 <- as.vector(freq_count_env2[,1])
freq_count_env2_2 <- as.vector(freq_count_env2[,2])
freq_count_env2_3 <- as.vector(freq_count_env2[,3])
freq_count_env2_4 <- as.vector(freq_count_env2[,4])
freq_count_env2_5 <- as.vector(freq_count_env2[,5])
freq_count_env2_6 <- as.vector(freq_count_env2[,6])
freq_count_env2_7 <- as.vector(freq_count_env2[,7])
freq_count_env2_8 <- as.vector(freq_count_env2[,8])
freq_count_env2_9 <- as.vector(freq_count_env2[,9])
freq_count_env2_10 <- as.vector(freq_count_env2[,10])
freq_count_env2_11 <- as.vector(freq_count_env2[,11])
freq_count_env2_12 <- as.vector(freq_count_env2[,12])

freq_count_env2_1<-freq_count_env2_1[!is.na(freq_count_env2_1)]
freq_count_env2_2<-freq_count_env2_2[!is.na(freq_count_env2_2)]
freq_count_env2_3<-freq_count_env2_3[!is.na(freq_count_env2_3)]
freq_count_env2_4<-freq_count_env2_4[!is.na(freq_count_env2_4)]
freq_count_env2_5<-freq_count_env2_5[!is.na(freq_count_env2_5)]
freq_count_env2_6<-freq_count_env2_6[!is.na(freq_count_env2_6)]
freq_count_env2_7<-freq_count_env2_7[!is.na(freq_count_env2_7)]
freq_count_env2_8<-freq_count_env2_8[!is.na(freq_count_env2_8)]
freq_count_env2_9<-freq_count_env2_9[!is.na(freq_count_env2_9)]
freq_count_env2_10<-freq_count_env2_10[!is.na(freq_count_env2_10)]
freq_count_env2_11<-freq_count_env2_11[!is.na(freq_count_env2_11)]
freq_count_env2_12<-freq_count_env2_12[!is.na(freq_count_env2_12)]


freq_count_env5_1 <- as.vector(freq_count_env5[,1])
freq_count_env5_2 <- as.vector(freq_count_env5[,2])
freq_count_env5_3 <- as.vector(freq_count_env5[,3])
freq_count_env5_4 <- as.vector(freq_count_env5[,4])
freq_count_env5_5 <- as.vector(freq_count_env5[,5])
freq_count_env5_6 <- as.vector(freq_count_env5[,6])
freq_count_env5_7 <- as.vector(freq_count_env5[,7])
freq_count_env5_8 <- as.vector(freq_count_env5[,8])
freq_count_env5_9 <- as.vector(freq_count_env5[,9])
freq_count_env5_10 <- as.vector(freq_count_env5[,10])
freq_count_env5_11 <- as.vector(freq_count_env5[,11])
freq_count_env5_12 <- as.vector(freq_count_env5[,12])

freq_count_env5_1<-freq_count_env5_1[!is.na(freq_count_env5_1)]
freq_count_env5_2<-freq_count_env5_2[!is.na(freq_count_env5_2)]
freq_count_env5_3<-freq_count_env5_3[!is.na(freq_count_env5_3)]
freq_count_env5_4<-freq_count_env5_4[!is.na(freq_count_env5_4)]
freq_count_env5_5<-freq_count_env5_5[!is.na(freq_count_env5_5)]
freq_count_env5_6<-freq_count_env5_6[!is.na(freq_count_env5_6)]
freq_count_env5_7<-freq_count_env5_7[!is.na(freq_count_env5_7)]
freq_count_env5_8<-freq_count_env5_8[!is.na(freq_count_env5_8)]
freq_count_env5_9<-freq_count_env5_9[!is.na(freq_count_env5_9)]
freq_count_env5_10<-freq_count_env5_10[!is.na(freq_count_env5_10)]
freq_count_env5_11<-freq_count_env5_11[!is.na(freq_count_env5_11)]
freq_count_env5_12<-freq_count_env5_12[!is.na(freq_count_env5_12)]


#######################################################################################################
#Make SNP table for neutral loci

#Filter pop_order_V to call for 2010 regional basetime V column names
pop_V1 <- pop_order_V %>% filter(Year==2010 & Pop==1)
pop_V2 <- pop_order_V %>% filter(Year==2010 & Pop==2)
pop_V3 <- pop_order_V %>% filter(Year==2010 & Pop==3)
pop_V4 <- pop_order_V %>% filter(Year==2010 & Pop==4)
pop_V5 <- pop_order_V %>% filter(Year==2010 & Pop==5)
pop_V6 <- pop_order_V %>% filter(Year==2010 & Pop==6)

pop_V7 <- pop_order_V %>% filter(Year==2010 & Pop==7)
pop_V8 <- pop_order_V %>% filter(Year==2011 & Pop==8)
pop_V9 <- pop_order_V %>% filter(Year==2010 & Pop==9)
pop_V10 <- pop_order_V %>% filter(Year==2011 & Pop==10)
pop_V11 <- pop_order_V %>% filter(Year==2010 & Pop==11)
pop_V12 <- pop_order_V %>% filter(Year==2010 & Pop==11)


#Filter Baseline A and B numbers for each population
loci_base_1 <- snp_swiss %>% select(chr_snp,all_of(pop_V1$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_2 <- snp_swiss %>% select(chr_snp,all_of(pop_V2$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_3 <- snp_swiss %>% select(chr_snp,all_of(pop_V3$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_4 <- snp_swiss %>% select(chr_snp,all_of(pop_V4$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_5 <- snp_swiss %>% select(chr_snp,all_of(pop_V5$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_6 <- snp_swiss %>% select(chr_snp,all_of(pop_V6$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))

loci_base_7 <- snp_swiss %>% select(chr_snp,all_of(pop_V7$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_8 <- snp_swiss %>% select(chr_snp,all_of(pop_V8$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_9 <- snp_swiss %>% select(chr_snp,all_of(pop_V9$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_10 <- snp_swiss %>% select(chr_snp,all_of(pop_V10$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_11 <- snp_swiss %>% select(chr_snp,all_of(pop_V11$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))
loci_base_12 <- snp_swiss %>% select(chr_snp,all_of(pop_V12$V_ID)) %>% mutate(total_freq= rowSums(across(where(is.numeric))))


#Calculate frequency for SNP A for 12 basetime pops
#Gives all basetime frequencies per region
p1A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p2A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p3A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p4A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p5A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p6A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)

p7A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p8A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p9A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p10A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p11A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)
p12A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(total_freq)) %>% select(chr_snp,snpA)



#Filter full snp table to remove climate associated SNPs
rand_slope_env1_out<-data.frame()
rand_slope_env2_out<-data.frame()
rand_slope_env5_out<-data.frame()


#######################################################################################################
# Stratified random sampling SNPs for a given population
for (seed_num in 1:1000){
set.seed(seed_num)
#env1
rand_env1_pop1 <- large(p1A,freq_count_env1_1)
rand_env1_pop2 <- large(p2A,freq_count_env1_2)
rand_env1_pop3 <- large(p3A,freq_count_env1_3)
rand_env1_pop4 <- large(p4A,freq_count_env1_4)
rand_env1_pop5 <- large(p5A,freq_count_env1_5)
rand_env1_pop6 <- large(p6A,freq_count_env1_6)
rand_env1_pop7 <- large(p7A,freq_count_env1_7)
rand_env1_pop8 <- large(p8A,freq_count_env1_8)
rand_env1_pop9 <- large(p9A,freq_count_env1_9)
rand_env1_pop10 <- large(p10A,freq_count_env1_10)
rand_env1_pop11 <- large(p11A,freq_count_env1_11)
rand_env1_pop12 <- large(p12A,freq_count_env1_12)

rand_env2_pop1 <- large(p1A,freq_count_env2_1)
rand_env2_pop2 <- large(p2A,freq_count_env2_2)
rand_env2_pop3 <- large(p3A,freq_count_env2_3)
rand_env2_pop4 <- large(p4A,freq_count_env2_4)
rand_env2_pop5 <- large(p5A,freq_count_env2_5)
rand_env2_pop6 <- large(p6A,freq_count_env2_6)
rand_env2_pop7 <- large(p7A,freq_count_env2_7)
rand_env2_pop8 <- large(p8A,freq_count_env2_8)
rand_env2_pop9 <- large(p9A,freq_count_env2_9)
rand_env2_pop10 <- large(p10A,freq_count_env2_10)
rand_env2_pop11 <- large(p11A,freq_count_env2_11)
rand_env2_pop12 <- large(p12A,freq_count_env2_12)

rand_env5_pop1 <- large(p1A,freq_count_env5_1)
rand_env5_pop2 <- large(p2A,freq_count_env5_2)
rand_env5_pop3 <- large(p3A,freq_count_env5_3)
rand_env5_pop4 <- large(p4A,freq_count_env5_4)
rand_env5_pop5 <- large(p5A,freq_count_env5_5)
rand_env5_pop6 <- large(p6A,freq_count_env5_6)
rand_env5_pop7 <- large(p7A,freq_count_env5_7)
rand_env5_pop8 <- large(p8A,freq_count_env5_8)
rand_env5_pop9 <- large(p9A,freq_count_env5_9)
rand_env5_pop10 <- large(p10A,freq_count_env5_10)
rand_env5_pop11 <- large(p11A,freq_count_env5_11)
rand_env5_pop12 <- large(p12A,freq_count_env5_12)

#Get Full timeseires for each randomly selected neutral location
rand_time_AB_env1_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop1$chr_snp)) #env1
rand_time_AB_env1_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop2$chr_snp)) #env1
rand_time_AB_env1_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop3$chr_snp)) #env1
rand_time_AB_env1_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop4$chr_snp)) #env1
rand_time_AB_env1_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop5$chr_snp)) #env1
rand_time_AB_env1_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop6$chr_snp)) #env1
rand_time_AB_env1_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop7$chr_snp)) #env1
rand_time_AB_env1_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop8$chr_snp)) #env1
rand_time_AB_env1_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop9$chr_snp)) #env1
rand_time_AB_env1_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop10$chr_snp)) #env1
rand_time_AB_env1_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop11$chr_snp)) #env1
rand_time_AB_env1_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env1_pop12$chr_snp)) #env1

rand_time_AB_env2_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop1$chr_snp)) #env2
rand_time_AB_env2_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop2$chr_snp)) #env2
rand_time_AB_env2_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop3$chr_snp)) #env2
rand_time_AB_env2_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop4$chr_snp)) #env2
rand_time_AB_env2_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop5$chr_snp)) #env2
rand_time_AB_env2_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop6$chr_snp)) #env2
rand_time_AB_env2_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop7$chr_snp)) #env2
rand_time_AB_env2_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop8$chr_snp)) #env2
rand_time_AB_env2_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop9$chr_snp)) #env2
rand_time_AB_env2_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop10$chr_snp)) #env2
rand_time_AB_env2_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop11$chr_snp)) #env2
rand_time_AB_env2_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env2_pop12$chr_snp)) #env2

rand_time_AB_env5_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop1$chr_snp)) #env5
rand_time_AB_env5_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop2$chr_snp)) #env5
rand_time_AB_env5_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop3$chr_snp)) #env5
rand_time_AB_env5_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop4$chr_snp)) #env5
rand_time_AB_env5_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop5$chr_snp)) #env5
rand_time_AB_env5_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop6$chr_snp)) #env5
rand_time_AB_env5_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop7$chr_snp)) #env5
rand_time_AB_env5_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop8$chr_snp)) #env5
rand_time_AB_env5_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop9$chr_snp)) #env5
rand_time_AB_env5_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop10$chr_snp)) #env5
rand_time_AB_env5_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop11$chr_snp)) #env5
rand_time_AB_env5_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop12$chr_snp)) #env5

# Calc SNP A and transpose
rand_env1_pop1 <- prop_A(rand_time_AB_env1_1,pop_order_2)
rand_env1_pop2 <- prop_A(rand_time_AB_env1_2,pop_order_2)
rand_env1_pop3 <- prop_A(rand_time_AB_env1_3,pop_order_2)
rand_env1_pop4 <- prop_A(rand_time_AB_env1_4,pop_order_2)
rand_env1_pop5 <- prop_A(rand_time_AB_env1_5,pop_order_2)
rand_env1_pop6 <- prop_A(rand_time_AB_env1_6,pop_order_2)
rand_env1_pop7 <- prop_A(rand_time_AB_env1_7,pop_order_2)
rand_env1_pop8 <- prop_A(rand_time_AB_env1_8,pop_order_2)
rand_env1_pop9 <- prop_A(rand_time_AB_env1_9,pop_order_2)
rand_env1_pop10 <- prop_A(rand_time_AB_env1_10,pop_order_2)
rand_env1_pop11 <- prop_A(rand_time_AB_env1_11,pop_order_2)
rand_env1_pop12 <- prop_A(rand_time_AB_env1_12,pop_order_2)

rand_env2_pop1 <- prop_A(rand_time_AB_env2_1,pop_order_2)
rand_env2_pop2 <- prop_A(rand_time_AB_env2_2,pop_order_2)
rand_env2_pop3 <- prop_A(rand_time_AB_env2_3,pop_order_2)
rand_env2_pop4 <- prop_A(rand_time_AB_env2_4,pop_order_2)
rand_env2_pop5 <- prop_A(rand_time_AB_env2_5,pop_order_2)
rand_env2_pop6 <- prop_A(rand_time_AB_env2_6,pop_order_2)
rand_env2_pop7 <- prop_A(rand_time_AB_env2_7,pop_order_2)
rand_env2_pop8 <- prop_A(rand_time_AB_env2_8,pop_order_2)
rand_env2_pop9 <- prop_A(rand_time_AB_env2_9,pop_order_2)
rand_env2_pop10 <- prop_A(rand_time_AB_env2_10,pop_order_2)
rand_env2_pop11 <- prop_A(rand_time_AB_env2_11,pop_order_2)
rand_env2_pop12 <- prop_A(rand_time_AB_env2_12,pop_order_2)

rand_env5_pop1 <- prop_A(rand_time_AB_env5_1,pop_order_2)
rand_env5_pop2 <- prop_A(rand_time_AB_env5_2,pop_order_2)
rand_env5_pop3 <- prop_A(rand_time_AB_env5_3,pop_order_2)
rand_env5_pop4 <- prop_A(rand_time_AB_env5_4,pop_order_2)
rand_env5_pop5 <- prop_A(rand_time_AB_env5_5,pop_order_2)
rand_env5_pop6 <- prop_A(rand_time_AB_env5_6,pop_order_2)
rand_env5_pop7 <- prop_A(rand_time_AB_env5_7,pop_order_2)
rand_env5_pop8 <- prop_A(rand_time_AB_env5_8,pop_order_2)
rand_env5_pop9 <- prop_A(rand_time_AB_env5_9,pop_order_2)
rand_env5_pop10 <- prop_A(rand_time_AB_env5_10,pop_order_2)
rand_env5_pop11 <- prop_A(rand_time_AB_env5_11,pop_order_2)
rand_env5_pop12 <- prop_A(rand_time_AB_env5_12,pop_order_2)

# Filter for site (its currently duplicated since timeseries filtering was on 12 pop datatset)
rand_env1_pop1 <- rand_env1_pop1 %>% filter(Site==1)
rand_env1_pop2 <- rand_env1_pop2 %>% filter(Site==2)
rand_env1_pop3 <- rand_env1_pop3 %>% filter(Site==3)
rand_env1_pop4 <- rand_env1_pop4 %>% filter(Site==4)
rand_env1_pop5 <- rand_env1_pop5 %>% filter(Site==5)
rand_env1_pop6 <- rand_env1_pop6 %>% filter(Site==6)
rand_env1_pop7 <- rand_env1_pop7 %>% filter(Site==7)
rand_env1_pop8 <- rand_env1_pop8 %>% filter(Site==8)
rand_env1_pop9 <- rand_env1_pop9 %>% filter(Site==9)
rand_env1_pop10 <- rand_env1_pop10 %>% filter(Site==10)
rand_env1_pop11 <- rand_env1_pop11 %>% filter(Site==11)
rand_env1_pop12 <- rand_env1_pop12 %>% filter(Site==12)

rand_env2_pop1 <- rand_env2_pop1 %>% filter(Site==1)
rand_env2_pop2 <- rand_env2_pop2 %>% filter(Site==2)
rand_env2_pop3 <- rand_env2_pop3 %>% filter(Site==3)
rand_env2_pop4 <- rand_env2_pop4 %>% filter(Site==4)
rand_env2_pop5 <- rand_env2_pop5 %>% filter(Site==5)
rand_env2_pop6 <- rand_env2_pop6 %>% filter(Site==6)
rand_env2_pop7 <- rand_env2_pop7 %>% filter(Site==7)
rand_env2_pop8 <- rand_env2_pop8 %>% filter(Site==8)
rand_env2_pop9 <- rand_env2_pop9 %>% filter(Site==9)
rand_env2_pop10 <- rand_env2_pop10 %>% filter(Site==10)
rand_env2_pop11 <- rand_env2_pop11 %>% filter(Site==11)
rand_env2_pop12 <- rand_env2_pop12 %>% filter(Site==12)

rand_env5_pop1 <- rand_env5_pop1 %>% filter(Site==1)
rand_env5_pop2 <- rand_env5_pop2 %>% filter(Site==2)
rand_env5_pop3 <- rand_env5_pop3 %>% filter(Site==3)
rand_env5_pop4 <- rand_env5_pop4 %>% filter(Site==4)
rand_env5_pop5 <- rand_env5_pop5 %>% filter(Site==5)
rand_env5_pop6 <- rand_env5_pop6 %>% filter(Site==6)
rand_env5_pop7 <- rand_env5_pop7 %>% filter(Site==7)
rand_env5_pop8 <- rand_env5_pop8 %>% filter(Site==8)
rand_env5_pop9 <- rand_env5_pop9 %>% filter(Site==9)
rand_env5_pop10 <- rand_env5_pop10 %>% filter(Site==10)
rand_env5_pop11 <- rand_env5_pop11 %>% filter(Site==11)
rand_env5_pop12 <- rand_env5_pop12 %>% filter(Site==12)

#Get slopes
rand_slope_env1_pop1 <- glm_slopes(rand_env1_pop1)
rand_slope_env1_pop2 <- glm_slopes(rand_env1_pop2)
rand_slope_env1_pop3 <- glm_slopes(rand_env1_pop3)
rand_slope_env1_pop4 <- glm_slopes(rand_env1_pop4)
rand_slope_env1_pop5 <- glm_slopes(rand_env1_pop5)
rand_slope_env1_pop6 <- glm_slopes(rand_env1_pop6)
rand_slope_env1_pop7 <- glm_slopes(rand_env1_pop7)
rand_slope_env1_pop8 <- glm_slopes(rand_env1_pop8)
rand_slope_env1_pop9 <- glm_slopes(rand_env1_pop9)
rand_slope_env1_pop10 <- glm_slopes(rand_env1_pop10)
rand_slope_env1_pop11 <- glm_slopes(rand_env1_pop11)
rand_slope_env1_pop12 <- glm_slopes(rand_env1_pop12)

rand_slope_env2_pop1 <- glm_slopes(rand_env2_pop1)
rand_slope_env2_pop2 <- glm_slopes(rand_env2_pop2)
rand_slope_env2_pop3 <- glm_slopes(rand_env2_pop3)
rand_slope_env2_pop4 <- glm_slopes(rand_env2_pop4)
rand_slope_env2_pop5 <- glm_slopes(rand_env2_pop5)
rand_slope_env2_pop6 <- glm_slopes(rand_env2_pop6)
rand_slope_env2_pop7 <- glm_slopes(rand_env2_pop7)
rand_slope_env2_pop8 <- glm_slopes(rand_env2_pop8)
rand_slope_env2_pop9 <- glm_slopes(rand_env2_pop9)
rand_slope_env2_pop10 <- glm_slopes(rand_env2_pop10)
rand_slope_env2_pop11 <- glm_slopes(rand_env2_pop11)
rand_slope_env2_pop12 <- glm_slopes(rand_env2_pop12)

rand_slope_env5_pop1 <- glm_slopes(rand_env5_pop1)
rand_slope_env5_pop2 <- glm_slopes(rand_env5_pop2)
rand_slope_env5_pop3 <- glm_slopes(rand_env5_pop3)
rand_slope_env5_pop4 <- glm_slopes(rand_env5_pop4)
rand_slope_env5_pop5 <- glm_slopes(rand_env5_pop5)
rand_slope_env5_pop6 <- glm_slopes(rand_env5_pop6)
rand_slope_env5_pop7 <- glm_slopes(rand_env5_pop7)
rand_slope_env5_pop8 <- glm_slopes(rand_env5_pop8)
rand_slope_env5_pop9 <- glm_slopes(rand_env5_pop9)
rand_slope_env5_pop10 <- glm_slopes(rand_env5_pop10)
rand_slope_env5_pop11 <- glm_slopes(rand_env5_pop11)
rand_slope_env5_pop12 <- glm_slopes(rand_env5_pop12)


#Bind populations for each env
rand_slope_env1 <- rbind(rand_slope_env1_pop1,
                      rand_slope_env1_pop2,
                      rand_slope_env1_pop3,
                      rand_slope_env1_pop4,
                      rand_slope_env1_pop5,
                      rand_slope_env1_pop6,
                      rand_slope_env1_pop7,
                      rand_slope_env1_pop8,
                      rand_slope_env1_pop9,
                      rand_slope_env1_pop10,
                      rand_slope_env1_pop11,
                      rand_slope_env1_pop12)

rand_slope_env2 <- rbind(rand_slope_env2_pop1,
                      rand_slope_env2_pop2,
                      rand_slope_env2_pop3,
                      rand_slope_env2_pop4,
                      rand_slope_env2_pop5,
                      rand_slope_env2_pop6,
                      rand_slope_env2_pop7,
                      rand_slope_env2_pop8,
                      rand_slope_env2_pop9,
                      rand_slope_env2_pop10,
                      rand_slope_env2_pop11,
                      rand_slope_env2_pop12)

rand_slope_env5 <- rbind(rand_slope_env5_pop1,
                      rand_slope_env5_pop2,
                      rand_slope_env5_pop3,
                      rand_slope_env5_pop4,
                      rand_slope_env5_pop5,
                      rand_slope_env5_pop6,
                      rand_slope_env5_pop7,
                      rand_slope_env5_pop8,
                      rand_slope_env5_pop9,
                      rand_slope_env5_pop10,
                      rand_slope_env5_pop11,
                      rand_slope_env5_pop12)
rand_slope_env1 <- rand_slope_env1 %>% mutate (Seed_ID = seed_num)
rand_slope_env2 <- rand_slope_env2 %>% mutate (Seed_ID = seed_num)
rand_slope_env5 <- rand_slope_env5 %>% mutate (Seed_ID = seed_num)
#Export each joint df

rand_slope_env1_out<-rbind(rand_slope_env1_out,rand_slope_env1)
rand_slope_env2_out<-rbind(rand_slope_env2_out,rand_slope_env2)
rand_slope_env5_out<-rbind(rand_slope_env5_out,rand_slope_env5)

print(seed_num)

}

#Save large files in folder outside of github
setwd("~/Dropbox/AM_Workshop/Large_files")

write_csv(rand_slope_env1_out, "rand_slope_env1.csv")
write_csv(rand_slope_env2_out, "rand_slope_env2.csv")
write_csv(rand_slope_env5_out, "rand_slope_env5.csv")

setwd("~/Dropbox/AM_Workshop/AM_Workshop/")























rand_slope_env5_out<-data.frame()

#env5 ONLY

#######################################################################################################
# Stratified random sampling SNPs for a given population
for (seed_num in 1:1000){
  set.seed(seed_num)
  
  rand_env5_pop1 <- large(p1A,freq_count_env5_1)
  rand_env5_pop2 <- large(p2A,freq_count_env5_2)
  rand_env5_pop3 <- large(p3A,freq_count_env5_3)
  rand_env5_pop4 <- large(p4A,freq_count_env5_4)
  rand_env5_pop5 <- large(p5A,freq_count_env5_5)
  rand_env5_pop6 <- large(p6A,freq_count_env5_6)
  rand_env5_pop7 <- large(p7A,freq_count_env5_7)
  rand_env5_pop8 <- large(p8A,freq_count_env5_8)
  rand_env5_pop9 <- large(p9A,freq_count_env5_9)
  rand_env5_pop10 <- large(p10A,freq_count_env5_10)
  rand_env5_pop11 <- large(p11A,freq_count_env5_11)
  rand_env5_pop12 <- large(p12A,freq_count_env5_12)
  
  #Get Full timeseires for each randomly selected neutral location
  rand_time_AB_env5_1 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop1$chr_snp)) #env5
  rand_time_AB_env5_2 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop2$chr_snp)) #env5
  rand_time_AB_env5_3 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop3$chr_snp)) #env5
  rand_time_AB_env5_4 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop4$chr_snp)) #env5
  rand_time_AB_env5_5 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop5$chr_snp)) #env5
  rand_time_AB_env5_6 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop6$chr_snp)) #env5
  rand_time_AB_env5_7 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop7$chr_snp)) #env5
  rand_time_AB_env5_8 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop8$chr_snp)) #env5
  rand_time_AB_env5_9 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop9$chr_snp)) #env5
  rand_time_AB_env5_10 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop10$chr_snp)) #env5
  rand_time_AB_env5_11 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop11$chr_snp)) #env5
  rand_time_AB_env5_12 <-loci_snp %>% filter (chr_snp %in% as.character(rand_env5_pop12$chr_snp)) #env5
  
  # Calc SNP A and transpose
  rand_env5_pop1 <- prop_A(rand_time_AB_env5_1,pop_order_2)
  rand_env5_pop2 <- prop_A(rand_time_AB_env5_2,pop_order_2)
  rand_env5_pop3 <- prop_A(rand_time_AB_env5_3,pop_order_2)
  rand_env5_pop4 <- prop_A(rand_time_AB_env5_4,pop_order_2)
  rand_env5_pop5 <- prop_A(rand_time_AB_env5_5,pop_order_2)
  rand_env5_pop6 <- prop_A(rand_time_AB_env5_6,pop_order_2)
  rand_env5_pop7 <- prop_A(rand_time_AB_env5_7,pop_order_2)
  rand_env5_pop8 <- prop_A(rand_time_AB_env5_8,pop_order_2)
  rand_env5_pop9 <- prop_A(rand_time_AB_env5_9,pop_order_2)
  rand_env5_pop10 <- prop_A(rand_time_AB_env5_10,pop_order_2)
  rand_env5_pop11 <- prop_A(rand_time_AB_env5_11,pop_order_2)
  rand_env5_pop12 <- prop_A(rand_time_AB_env5_12,pop_order_2)
  
  # Filter for site (its currently duplicated since timeseries filtering was on 12 pop datatset)
  rand_env5_pop1 <- rand_env5_pop1 %>% filter(Site==1)
  rand_env5_pop2 <- rand_env5_pop2 %>% filter(Site==2)
  rand_env5_pop3 <- rand_env5_pop3 %>% filter(Site==3)
  rand_env5_pop4 <- rand_env5_pop4 %>% filter(Site==4)
  rand_env5_pop5 <- rand_env5_pop5 %>% filter(Site==5)
  rand_env5_pop6 <- rand_env5_pop6 %>% filter(Site==6)
  rand_env5_pop7 <- rand_env5_pop7 %>% filter(Site==7)
  rand_env5_pop8 <- rand_env5_pop8 %>% filter(Site==8)
  rand_env5_pop9 <- rand_env5_pop9 %>% filter(Site==9)
  rand_env5_pop10 <- rand_env5_pop10 %>% filter(Site==10)
  rand_env5_pop11 <- rand_env5_pop11 %>% filter(Site==11)
  rand_env5_pop12 <- rand_env5_pop12 %>% filter(Site==12)
  
  #Get slopes
  rand_slope_env5_pop1 <- glm_slopes(rand_env5_pop1)
  rand_slope_env5_pop2 <- glm_slopes(rand_env5_pop2)
  rand_slope_env5_pop3 <- glm_slopes(rand_env5_pop3)
  rand_slope_env5_pop4 <- glm_slopes(rand_env5_pop4)
  rand_slope_env5_pop5 <- glm_slopes(rand_env5_pop5)
  rand_slope_env5_pop6 <- glm_slopes(rand_env5_pop6)
  rand_slope_env5_pop7 <- glm_slopes(rand_env5_pop7)
  rand_slope_env5_pop8 <- glm_slopes(rand_env5_pop8)
  rand_slope_env5_pop9 <- glm_slopes(rand_env5_pop9)
  rand_slope_env5_pop10 <- glm_slopes(rand_env5_pop10)
  rand_slope_env5_pop11 <- glm_slopes(rand_env5_pop11)
  rand_slope_env5_pop12 <- glm_slopes(rand_env5_pop12)
  
  
  #Bind populations for each env
  rand_slope_env5 <- rbind(rand_slope_env5_pop1,
                          rand_slope_env5_pop2,
                          rand_slope_env5_pop3,
                          rand_slope_env5_pop4,
                          rand_slope_env5_pop5,
                          rand_slope_env5_pop6,
                          rand_slope_env5_pop7,
                          rand_slope_env5_pop8,
                          rand_slope_env5_pop9,
                          rand_slope_env5_pop10,
                          rand_slope_env5_pop11,
                          rand_slope_env5_pop12)
  rand_slope_env5 <- rand_slope_env5 %>% mutate (Seed_ID = seed_num)
  #Export each joint df
  rand_slope_env5_out<-rbind(rand_slope_env5_out,rand_slope_env5)
  
  print(seed_num)
  
}









