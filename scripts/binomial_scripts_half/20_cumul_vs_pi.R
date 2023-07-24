##################################################################################
## Plot cumulative slope against demography
## Author Daniel Anstett
## 
## 
## Last Modified April 19, 2023
###################################################################################

#Library install and import
library(tidyverse)
library(car)
library(RColorBrewer)

#Import data
offset_pop <- read_csv("data/binomial_data/time_cumul_beagle.csv")
pi.df <- read_csv("data/binomial_data/raw_pi.csv")

pi_cumul <- left_join(offset_pop,pi.df,by=c("Paper_ID"="Site")) 

#Models
lm1 <- lm(cumul_pos~pi_snp_set,data=pi_cumul)
lm2 <- lm(cumul_all~pi_snp_set,data=pi_cumul)
lm3 <- lm(cumul_pos~pi_all_snps,data=pi_cumul)
lm4 <- lm(cumul_all~pi_all_snps,data=pi_cumul)

#Anova

Anova(lm1,type="III")
summary(lm1)
Anova(lm2,type="III")
summary(lm2)
Anova(lm3,type="III")
summary(lm3)
Anova(lm4,type="III")
summary(lm4)




###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(offset_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#cumul slope plotted against 1215 offset
ggplot(pi_cumul, aes(x=pi_snp_set, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="PI Climate Associated")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plot_pi/1_pi_cumul1.pdf",width=8, height = 6, units = "in")


#cumul slope plotted against 1215 offset
ggplot(pi_cumul, aes(x=pi_snp_set, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="PI Climate Associated")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plot_pi/2_pi_cumul2.pdf",width=8, height = 6, units = "in")


###########################################################################################################


#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
ggplot(pi_cumul, aes(x=pi_all_snps, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="PI Genome-Wide")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plot_pi/3_pi_cumul3.pdf",width=8, height = 6, units = "in")

#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
ggplot(pi_cumul, aes(x=pi_all_snps, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="PI Genome-Wide")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plot_pi/4_pi_cumul4.pdf",width=8, height = 6, units = "in")
