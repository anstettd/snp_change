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
#Import data
offset_pop <- read_csv("data/time_cumul.csv")
offset_pop_3 <- offset_pop %>% filter(Paper_ID!=3)


lm.cumul.1215 <- lm(cumul_slope~offset_1215,data=offset_pop)
lm.cumul.ssp245 <- lm(cumul_slope~offset_SSP245,data=offset_pop)
lm.cumul.ssp585 <- lm(cumul_slope~offset_SSP585,data=offset_pop)
lm.cumul.clim <- lm(cumul_slope~offset_climate,data=offset_pop)

lm.cumul.1215_no3 <- lm(cumul_slope~offset_1215,data=offset_pop_3)
lm.cumul.clim_no3 <- lm(cumul_slope~offset_climate,data=offset_pop_3)


Anova(lm.cumul.1215 ,type="III")
Anova(lm.cumul.ssp245 ,type="III")
Anova(lm.cumul.ssp585 ,type="III")
Anova(lm.cumul.clim ,type="III")

Anova(lm.cumul.1215_no3 ,type="III")
Anova(lm.cumul.clim_no3 ,type="III")


###########################################################################################################

#cumul slope plotted against 1215 offset
ggplot(offset_pop, aes(x=offset_1215, y=cumul_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Rapid Evolution")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/3_rapid_offset_1215.pdf",width=7, height = 6, units = "in")

#cumul slope plotted against 1215 offset
ggplot(offset_pop_3, aes(x=offset_1215, y=cumul_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Rapid Evolution")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/4_rapid_offset_1215.pdf",width=7, height = 6, units = "in")




###########################################################################################################

#cumul slope plotted against 2012-2015 Climate Distance
ggplot(offset_pop, aes(x=offset_climate, y=cumul_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Rapid Evolution")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/5_rapid_offset_climate_distance.pdf",width=7, height = 6, units = "in")


#cumul slope plotted against 2012-2015 Climate Distance
ggplot(offset_pop_3, aes(x=offset_climate, y=cumul_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Rapid Evolution")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/6_rapid_offset_climate_distance.pdf",width=7, height = 6, units = "in")


###########################################################################################################


#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP245, y=cumul_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Rapid Evolution")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/7_rapid_offset_ssp245.pdf",width=7, height = 6, units = "in")


#cumul slope plotted against 2040-2070 SSP585 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP585, y=cumul_slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Rapid Evolution")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/8_rapid_offset_ssp585.pdf",width=7, height = 6, units = "in")



