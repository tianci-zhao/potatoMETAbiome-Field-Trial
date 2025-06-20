# PotatoMETAbiom Field trial experiment
# NL - plant development - Week5

## load libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(readxl)


#Figure 1 ####
#load dataset for analysis
#biomass_w5.xlsx
df2 <- as.data.frame(read_excel("Biomass_w5.xlsx"))
df2$Treatment <- factor(df2$Treatment, levels=c("Control", "Consortium B", "Consortium BP", "Fertilizers", "Pesticides", "Fertilizers-Pesticides"))
df2$Variety<- factor(df2$Variety, levels=c("ATOL", "DESIREE", "JELLY", "KRAB", "PASJA POMORSKA", "RUDAWA", "SALTO"))

##varieties/treatment####
fig1 <- ggplot(df2) + 
  geom_boxplot(aes(x = Variety, y = Height, color = Variety,fill = Variety), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Variety, y = Height, fill = Variety, color = Variety), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Treatment),ncol = 100) +
  xlab("") +
  ylab("Height / cm")+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),  
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))
fig1

fig2 <- ggplot(df2) + 
  geom_boxplot(aes(x = Variety, y = Aboveground_biomass, color = Variety,fill = Variety), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Variety, y =  Aboveground_biomass , fill = Variety, color = Variety), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Treatment),ncol = 100) +
  xlab("") +
  ylab("Aboveground biomass / g")+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),   
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))
fig2

fig3 <- ggplot(df2) +
  geom_boxplot(aes(x = Variety, y =  Belowground_biomass, color = Variety,fill = Variety), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Variety, y =  Belowground_biomass, fill = Variety, color = Variety), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Treatment),ncol = 100) +
  xlab(" ") +
  ylab("Belowground biomass / g")+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),   
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

fig3


fig4 <- ggplot(df2) +
  geom_boxplot(aes(x = Variety, y =  root_shoot, color = Variety,fill = Variety), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Variety, y =  root_shoot, fill = Variety, color = Variety), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Treatment),ncol = 100) +
  xlab("Control Consortium B Consortium BP Fertilizers Pesticides Fertilizers-Pesticides") +
  ylab("root/shoot ratio")+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),   
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

fig4


fig1 + fig2 + fig3 + fig4+ plot_layout(nrow = 4, guides= "collect")&
  theme(legend.position = 'bottom')

#ggsave("biomass_root_shoot_var.pdf",fig1 + fig2 + fig3 + fig4 + plot_layout(nrow = 4, guides= "collect")&theme(legend.position = 'bottom'), width=15,height=12)


##Two-way ANOVA
library(agricolae)
library(tidyverse)
library(dplyr)
library(readxl)


res <- aov(df2$root_shoot ~ df2$Treatment*df2$Variety, data = df2)

summary(res)


#post hoc - treatment
df3 <- df2 %>% filter(Treatment == "Fertilizers-Pesticides")

re.lm <- lm(root_shoot ~ Variety, data = df3)
re.av <- aov(re.lm) #betweenness
summary(re.av)  

re.test <- duncan.test(re.av, trt = 'Variety')
re.test


#Supplementary Figure 2####
##treatments in same variety####
fig1 <- ggplot(df2) + 
  geom_boxplot(aes(x = Treatment, y = Height, color = Treatment,fill = Treatment), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Treatment, y = Height, fill = Treatment, color = Treatment), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Variety),ncol = 100) +
  xlab("") +
  ylab("Height / cm")+
  scale_color_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),   
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))
fig1

fig2 <- ggplot(df2) + 
  geom_boxplot(aes(x = Treatment, y = Aboveground_biomass, color = Treatment,fill = Treatment), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Treatment, y =  Aboveground_biomass , fill = Treatment, color = Treatment), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Variety),ncol = 100) +
  xlab("") +
  ylab("Aboveground biomass / g")+
  scale_color_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),   
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))
fig2

fig3 <- ggplot(df2) +
  geom_boxplot(aes(x = Treatment, y =  Belowground_biomass, color = Treatment,fill = Treatment), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Treatment, y =  Belowground_biomass, fill = Treatment, color = Treatment), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Variety),ncol = 100) +
  xlab(" ") +
  ylab("Belowground biomass / g")+
  scale_color_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),   
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

fig3


fig4 <- ggplot(df2) +
  geom_boxplot(aes(x = Treatment, y =  root_shoot, color = Treatment,fill = Treatment), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"),alpha = 0.7) + 
  geom_point(aes(x = Treatment, y =  root_shoot, fill = Treatment, color = Treatment), 
             size = 2, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  facet_wrap(. ~ factor(Variety),ncol = 100) +
  xlab("ATOL   DESIREE   JELLY   KRAB   PASJA POMORSKA   RUDAWA   SALTO") +
  ylab("root/shoot ratio")+
  scale_color_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  theme(axis.title.x = element_text(family="Arial", colour = "black", size=14),
        axis.title.y = element_text(family="Arial", colour = "black", size=14),
        axis.text.x  = element_blank(),  
        #axis.text.x  = element_text(family="Arial", colour = "black", size=14),
        axis.text.y  = element_text(family="Arial", colour = "black",size=14),
        axis.ticks.x = element_blank(),   
        legend.text=element_text(family="Arial", colour = "black",size=14),
        legend.title=element_text(family="Arial", colour = "black",size=14),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(family="Arial",colour = "black",size=14),
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  #theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1))+
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

fig4

#(p <- ggarrange(fig1, fig2, fig3,fig4, common.legend = TRUE, legend = "right", ncol = 1, nrow = 4, align = "hv"))

fig1 + fig2 + fig3 + fig4+ plot_layout(nrow = 4, guides= "collect")&
  theme(legend.position = 'bottom')

#ggsave("biomass_root_shoot_tre_1.pdf",fig1 + fig2 + fig3 + fig4 + plot_layout(nrow = 4, guides= "collect")&theme(legend.position = 'bottom'), width=15,height=12)

library(agricolae)
library(tidyverse)
library(dplyr)
library(readxl)
##Two-way ANOVA

res <- aov(df2$root_shoot ~ df2$Treatment * df2$Variety, data = df2)

summary(res)


#post hoc - Variety
df3 <- df2 %>% filter(Variety == "ATOL")

re.lm <- lm(root_shoot ~ Treatment, data = df3)
re.av <- aov(re.lm) #betweenness
summary(re.av)  

re.test <- duncan.test(re.av, trt = 'Treatment')
re.test

#Supplementary Figure 3####
##microbial alpha diversity####
#same code above Sulpmentary Figure 2, but use the alpha diversity as input data
df2 <- as.data.frame(read_excel("alpha_diversity_rhizo_w5_bacteria.xls"))
df2 <- as.data.frame(read_excel("alpha_diversity_rhizo_w5_fungi.xls"))

df2$Treatment <- factor(df2$Treatment, levels=c("Control", "Consortium B", "Consortium BP", "Fertiliser", "Pesticide", "Fertiliser-Pesticide"))
df2$Variety<- factor(df2$Variety, levels=c("ATOL", "DESIREE", "JELLY", "KRAB", "PASJA_POMORSKA", "RUDAWA", "SALTO"))


