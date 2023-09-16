##################################################################################
## Plot SNP frequency for all baseline populations
## Done for all SNP set that are in the timeseries
## 
## 
## Last Modified Sept 16, 2023
###################################################################################
#Import libraries
library(tidyverse)
library(egg)

###################################################################################
theme_spaghetti <- function(){ 
  
  theme_classic() %+replace%    #replace elements we want to change
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
          axis.title = element_text(size =16, face = "bold"), 
          axis.text.y = element_text(size = 14, face = "bold"),legend.title = element_blank(),
          legend.text = element_text(size=12,face="bold"),
          strip.background = element_blank(), 
          strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
}


###################################################################################
#Import SNP binomial data
abund_clim <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/data/abund_table_baseline.csv",
                       col_names = T)

###################################################################################

abund_clim$ENV <- as.factor(abund_clim$ENV)
abund_clim$ENV <- factor(abund_clim$ENV, levels = c("MAT","MAP","PAS","EXT","CMD",
                                                     "Tave_wt","Tave_sm","PPT_wt","PPT_sm"))

#Plot frequency vs latitude vs
ggplot(data=abund_clim,aes(Latitude,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.15,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Latitude") +  facet_wrap(.~ENV) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_baseline_lat.pdf",width=8, height = 7, units = "in")


#Plot frequency vs climate
plot_env_1<-ggplot(data=abund_clim,aes(MAT,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="MAT")  + theme_spaghetti()

plot_env_2<-ggplot(data=abund_clim,aes(MAP,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="MAP")  + theme_spaghetti()

plot_env_3<-ggplot(data=abund_clim,aes(PAS,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="PAS")  + theme_spaghetti()

plot_env_4<-ggplot(data=abund_clim,aes(EXT,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="EXT")  + theme_spaghetti()

plot_env_5<-ggplot(data=abund_clim,aes(CMD,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="CMD")  + theme_spaghetti()

plot_env_6<-ggplot(data=abund_clim,aes(Tave_wt,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Tave_wt")  + theme_spaghetti()

plot_env_7<-ggplot(data=abund_clim,aes(Tave_sm,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Tave_sm")  + theme_spaghetti()

plot_env_8<-ggplot(data=abund_clim,aes(PPT_wt,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="PPT_wt")  + theme_spaghetti()

plot_env_9<-ggplot(data=abund_clim,aes(PPT_sm,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="PPT_sm")  + theme_spaghetti()



arrg_1<-ggarrange(plot_env_1, 
                  plot_env_2 +theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank()),
                  plot_env_3 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()),                  
                  plot_env_4, 
                  plot_env_5 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()), 
                  plot_env_6 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()),
                  plot_env_7, 
                  plot_env_8 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()), 
                  plot_env_9 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()), 
                  nrow = 3, ncol = 3)

#Export 8 X 7

#ggsave(arrg_1,"graphs/spaghetii/spaghetii_baseline_env.pdf",width=8, height = 7, units = "in")







