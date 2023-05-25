##################################################################################
## Plot cumulative slope against demography
## Author Daniel Anstett
## 
## 
## Last Modified April 19, 2023
###################################################################################

#Library install and import
library(tidyverse) 
#Import data
offset_pop <- read_csv("data/demo_cumul.csv")



#stats
lm.cumul.slope <- lm(lambda.slope~cumul_slope,data=offset_pop)
lm.cumul.mean <- lm(lambda.mean~cumul_slope,data=offset_pop)


Anova(lm.cumul.slope ,type="III")
Anova(lm.cumul.mean,type="III")



###########################################################################################################

#cumul slope plotted against lambda slope
ggplot(offset_pop, aes(x=cumul_slope, y=lambda.slope, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_hline(yintercept=c(0,0), linetype="dotted")+
  #geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Lambda Slope")+
  scale_x_continuous(name="Cumulative Rapid Evolution")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/1_rapid_slope.pdf",width=7, height = 6, units = "in")

#cumul slope plotted against lambda slope
ggplot(offset_pop, aes(x=cumul_slope, y=lambda.mean, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Mean Lambda")+
  scale_x_continuous(name="Cumulative Rapid Evolution")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.85, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/2_rapid_mean.pdf",width=7, height = 6, units = "in")

