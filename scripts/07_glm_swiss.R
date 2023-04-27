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

#Melt example
slope_melt <- function(df) {
  freq_1.melted <- reshape2::melt(df, id.vars = c("Site", "Year"))
  colnames(freq_1.melted)[3] <- "snp_ID"
  colnames(freq_1.melted)[4] <- "freq"
  
  freq_1.melted_mod <- na.omit(freq_1.melted)
  
  freq_1.slope <- group_by(freq_1.melted_mod, Site, snp_ID) %>%
    arrange(Site, Year) %>%
    summarize(Slope = glm(as.numeric(freq) ~ as.numeric(Year), family = binomial)$coefficients[2]*2,
     SE = summary(glm(as.numeric(freq) ~ as.numeric(Year), family = binomial))$coefficients[2,2])
  
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

#Clear memory
rm(loci)
rm(loci_snp)
rm(snp)
rm(loci_united)



#################################################################################################
#Calculate frequencies from abundance

swiss_propA <- prop_A(snp_swiss,pop_order_2)
#rm(snp_swiss)

#Split data set into 100k ch
swiss_1 <- swiss_propA[1:100000,]




















#################################################################################################
##Filter by site & calc slopes and SE

swiss_propA_p1 <- swiss_propA %>% filter(Site==1)
swiss_glm_p1 <- slope_melt(swiss_propA_p1) #Run glm function
#rm(swiss_propA_p1)

swiss_propA_p2 <- swiss_propA %>% filter(Site==2)
swiss_propA_p3 <- swiss_propA %>% filter(Site==3)
swiss_propA_p4 <- swiss_propA %>% filter(Site==4)
swiss_propA_p5 <- swiss_propA %>% filter(Site==5)
swiss_propA_p6 <- swiss_propA %>% filter(Site==6)
swiss_propA_p7 <- swiss_propA %>% filter(Site==7)
swiss_propA_p8 <- swiss_propA %>% filter(Site==8)
swiss_propA_p9 <- swiss_propA %>% filter(Site==9)
swiss_propA_p10 <- swiss_propA %>% filter(Site==10)
swiss_propA_p11 <- swiss_propA %>% filter(Site==11)
swiss_propA_p12 <- swiss_propA %>% filter(Site==12)

###################################################################################




#write_csv(swiss_glm,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm.csv")

