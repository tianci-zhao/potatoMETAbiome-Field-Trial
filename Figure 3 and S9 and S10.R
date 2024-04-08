# PotatoMETAbiom Field trial experiment
# NL - microbial community composition - Week5

##Initiate libraries####
rm(list=ls())

library(vegan)
library(multcompView)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

#Figure 3####

### import bacterial distance data (Bray Curtis)####
mydata <- data.table::fread("bac_dist_T1.csv")

# analysis of variance; 1-DISTANCE = similarity
anova <- aov(1-Distance ~ Management, data = mydata)
summary(anova)    #<2e-16 ***

# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)
# compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# table with factors and 3rd quantile
Tk <- group_by(mydata, Management) %>%
  summarise(mean=mean(1-Distance), quant = quantile(1-Distance, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$Management)
Tk$cld <- cld$Letters

print(Tk)

# boxplot # similarity = (1 - Bray_Curtis dissimilarity)
p1 <- ggplot(mydata, aes(x=Management, y=1-Distance, fill=Management)) +
  geom_jitter(aes(color=Management),shape=16, alpha = 0.7, size = 2,
              position=position_jitter(0.2))+
  geom_boxplot(alpha = 0.9, outlier.shape = NA)+
  ylim(0, 1)+
  ylab("Similarity to Control")+ 
  xlab("")+
  labs(title = "Bacteria")+ #Bacterial community distance
  # geom_text(data = Tk, aes(x = Management, y = mean, label = cld),
  #size = 5, vjust=-8, hjust = 0.3)+
  scale_color_manual(values=c("#a6c39d","#ea8e94")) +
  scale_fill_manual(values=c("#a6c39d","#ea8e94"))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x =element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(size = 12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12),
        axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5))
p1

#ggsave("distance_T1_bac_man_1.jpeg", dpi=600, width=15, height=14,  units="cm")


### import fungal distance data####
mydata <- data.table::fread("fun_dist_T1.csv")

# analysis of variance; 1-DISTANCE = similarity
anova <- aov(1-Distance ~ Management, data = mydata)
summary(anova)    #<2e-16 ***

# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)
# compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# table with factors and 3rd quantile
Tk <- group_by(mydata, Management) %>%
  summarise(mean=mean(1-Distance), quant = quantile(1-Distance, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$Management)
Tk$cld <- cld$Letters

print(Tk)

# boxplot 
p2 <- ggplot(mydata, aes(x=Management, y=1-Distance, fill=Management)) +
  geom_jitter(aes(color=Management),shape=16, alpha = 0.7, size = 2,
              position=position_jitter(0.2))+
  geom_boxplot(alpha = 0.9, outlier.shape = NA)+
  ylim(0, 1)+
  ylab("Similarity to Control")+
  xlab("")+
  labs(title = "Fungi")+
  geom_text(data = Tk, aes(x = Management, y = mean, label = cld),
            size = 5, vjust=-6, hjust = 0.5)+
  scale_color_manual(values=c("#a6c39d","#ea8e94")) +
  scale_fill_manual(values=c("#a6c39d","#ea8e94"))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x =element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(size = 12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12),
        axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5))
p2

#ggsave("distance_T1_fun_man.jpeg", dpi=600, width=15, height=14,  units="cm")


Fig3 <- p1 +p2 + plot_layout(nrow = 1, guides= "collect")&
  theme(legend.position = 'bottom')

ggsave("distance_control_man_1.pdf",p1 +p2+ plot_layout(nrow = 1, guides= "collect")&
         theme(legend.position = 'bottom'),width=7,height=5)

##Figure S9####

### import bacterial distance data####
mydata <- data.table::fread("bac_dist_T1.csv")

mydata$Treatment <- factor(mydata$Treatment, 
                           levels =c("1", "2", "3", "4", "5", "6"), 
                           labels = c("Control", "Consortium B", "Consortium BP", "Fertiliser", "Pesticide",  "Fertiliser_Pesticide"))

# analysis of variance; 1-DISTANCE = similarity
anova <- aov(1-Distance ~ Treatment, data = mydata)
summary(anova)    #<2e-16 ***

# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)
# compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# table with factors and 3rd quantile
Tk <- group_by(mydata, Treatment) %>%
  summarise(mean=mean(1-Distance), quant = quantile(1-Distance, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$Treatment)
Tk$cld <- cld$Letters

print(Tk)

# boxplot similarity = (1 - Bray_Curtis dissimilarity)
p3 <- ggplot(mydata, aes(x=Treatment, y=1-Distance, fill=Treatment)) +
  geom_jitter(aes(color=Treatment),shape=16, alpha = 0.7, size = 2,
              position=position_jitter(0.2))+
  geom_boxplot(alpha = 0.9, outlier.shape = NA)+
  ylim(0, 1)+
  ylab("Similarity to Control")+ 
  xlab("")+
  labs(title = "Bacteria")+ #Bacterial community distance
  geom_text(data = Tk, aes(x = Treatment, y = mean, label = cld),
            size = 5, vjust=-4, hjust = 0.3)+
  scale_color_manual(values=c("#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x =element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(size = 12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
p3

### import fungal distance data####
mydata <- data.table::fread("fun_dist_T1.csv")

mydata$Treatment <- factor(mydata$Treatment, 
                           levels =c("1", "2", "3", "4", "5", "6"), 
                           labels = c("Control", "Consortium B", "Consortium BP", "Fertiliser", "Pesticide",  "Fertiliser_Pesticide"))

# analysis of variance; 1-DISTANCE = similarity
anova <- aov(1-Distance ~ Treatment, data = mydata)
summary(anova)    #<2e-16 ***

# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)
# compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# table with factors and 3rd quantile
Tk <- group_by(mydata, Treatment) %>%
  summarise(mean=mean(1-Distance), quant = quantile(1-Distance, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$Treatment)
Tk$cld <- cld$Letters

print(Tk)

# boxplot 
p4 <- ggplot(mydata, aes(x=Treatment, y=1-Distance, fill=Treatment)) +
  geom_jitter(aes(color=Treatment),shape=16, alpha = 0.7, size = 2,
              position=position_jitter(0.2))+
  geom_boxplot(alpha = 0.9, outlier.shape = NA)+
  ylim(0, 1)+
  ylab("Similarity to Control")+
  xlab("")+
  labs(title = "Fungi")+
  geom_text(data = Tk, aes(x = Treatment, y = mean, label = cld),
            size = 5, vjust=-4, hjust = 0.5)+#v 调节高度，h调节横向
  scale_color_manual(values=c("#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x =element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(size = 12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
p4

# merge bac and fun

p3 +p4 + plot_layout(nrow = 1, guides= "collect")&
  theme(legend.position = 'right')

ggsave("distance_control_4.pdf",p3 + p4 + plot_layout(nrow = 1, guides= "collect")&
         theme(legend.position = 'right'),width=9,height=5)


##Figure S10####

### import bacterial distance data####
mydata <- data.table::fread("bac_dist_all_variety.csv")

mydata$Variety <- factor(mydata$Variety,
                         levels=c("ATOL", "DESIREE", "JELLY",  "KRAB", "PASJA_POMORSKA", "RUDAWA", "SALTO"),
                         labels = c("ATOL", "DESIREE", "JELLY", "KRAB", "PASJA POMORSKA", "RUDAWA", "SALTO"))

# analysis of variance; 1-DISTANCE = similarity,DISTANCE = Dissimilarity
anova <- aov(Distance ~ Variety, data = mydata)
summary(anova)    #8.41e-07 ***

# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)

# compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# table with factors and 3rd quantile
Tk <- group_by(mydata, Variety) %>%
  summarise(mean=mean(Distance), quant = quantile(Distance, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$Variety)
Tk$cld <- cld$Letters

print(Tk)

# boxplot 
p5 <- ggplot(mydata, aes(x=Variety, y=Distance, fill=Variety)) +
  geom_jitter(aes(color=Variety),shape=16, alpha = 0.7, size = 2,
              position=position_jitter(0.2))+
  geom_boxplot(alpha = 0.9, outlier.shape = NA)+
  ylim(0, 1)+
  ylab("Dissimilarity")+
  xlab("")+
  labs(title = "Bacteria")+ #Bacterial community distance within varieties
  geom_text(data = Tk, aes(x = Variety, y = mean, label = cld),
            size = 5, vjust=-7, hjust = 0.5)+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x =element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(size = 12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
p5

### import fungal distance data####

mydata <- data.table::fread("fun_dist_all_variety.csv")

mydata$Variety <- factor(mydata$Variety, 
                         levels=c("ATOL", "DESIREE", "JELLY",  "KRAB", "PASJA_POMORSKA", "RUDAWA", "SALTO"),
                         labels = c("ATOL", "DESIREE", "JELLY", "KRAB", "PASJA POMORSKA", "RUDAWA", "SALTO"))

# analysis of variance; 1-DISTANCE = similarity, DISTANCE = Dissimilarity
anova <- aov(Distance ~ Variety, data = mydata)
summary(anova)    #<2e-16 ****

# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)
# compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# table with factors and 3rd quantile
Tk <- group_by(mydata, Variety) %>%
  summarise(mean=mean(Distance), quant = quantile(Distance, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$Variety)
Tk$cld <- cld$Letters

print(Tk)

# boxplot 
p6 <- ggplot(mydata, aes(x=Variety, y=Distance, fill=Variety)) +
  geom_jitter(aes(color=Variety),shape=16, alpha = 0.7, size = 2,
              position=position_jitter(0.2))+
  geom_boxplot(alpha = 0.9, outlier.shape = NA)+
  ylim(0, 1)+
  ylab("Dissimilarity")+
  xlab("")+
  labs(title = "Fungi")+
  geom_text(data = Tk, aes(x = Variety, y = mean, label = cld),
            size = 5, vjust=-7, hjust = 0.5)+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x =element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(size = 12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
p6

p5 + p6 + plot_layout(nrow = 1, guides= "collect")&
  theme(legend.position = 'right')

ggsave("distance_varieties_dissimi_1.pdf",p5 + p6 + plot_layout(nrow = 1, guides= "collect")&
         theme(legend.position = 'right'),width=10,height=4.5)



