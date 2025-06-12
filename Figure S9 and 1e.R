# PotatoMETAbiom Field trial experiment
# NL - microbial interactive traits (MITs) in previous controlled greenhouse experiment

##Initiate libraries####
rm(list=ls())

library(ggplot2)
library(dplyr)


#Figure S9 Performance of MIT-selected cultivars in a previous controlled greenhouse experiment.
#load metadata in vitro condition
df <- read.csv("metadata_for_zscore_GH.csv", header = 1, row.names = 1, sep = ",")

means <- colMeans(df, na.rm = TRUE)
sds <- apply(df, 2, sd, na.rm = TRUE)

#z-score
df_z <- scale(df, center = TRUE, scale = TRUE)

#save z-score data, manually add metadata, calculate mean value as the z-score of MITs
write.csv(df_z, "metadata_zscore_GH.csv")

#reload z-score data with replicates
data <- read.csv("metadata_zscore_GH.csv", header = 1, row.names = 1, sep = ",")
data$Cultivar <- reorder(data$Cultivar, data$Mean)
levels(data$Cultivar)

#average
mean_value <- mean(data$Mean, na.rm = TRUE)

#plot
fig1 <- ggplot(data) + 
  geom_boxplot(aes(x = Cultivar, y = Mean, color = Color), 
               lwd = 0.3, outlier.size = -1, 
               position = position_dodge(preserve = "single"), alpha = 0.7) + 
  geom_point(aes(x = Cultivar, y = Mean, fill = Color, color = Color), 
             size = 3, shape = 16, alpha = 0.9, 
             position = position_identity()) + 
  geom_hline(yintercept = mean_value, linetype = "dashed", color = "gray", lwd = 0.5) + 
  xlab("") +
  ylab("Rank of MITs (z score)")+ 
  scale_color_manual(values = c("#349839","#349839","#349839","#349839", "gray")) +  
  scale_fill_manual(values = c("#349839","#349839","#349839","#349839", "gray")) +
  theme_bw() +
  theme(text = element_text(size=12),
        plot.title = element_text(size=12, face="bold"),
        strip.background = element_blank(),
        strip.text=element_text(face="bold", size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

fig1


ggsave("FigureS1.pdf", width = 15, height = 8, units = "cm", fig1, scale = 1.6)


#Figure 1e The relationship between below-ground biomass in field trial and MIT scores of selected cultivars
library(ggpubr)
library(tidyverse)
library(ggpmisc)

mytheme <- 
  theme_bw() +
  theme(text = element_text(size=14),
        plot.title = element_text(size=14, face="bold"),
        strip.background = element_blank(),
        strip.text=element_text(face="bold", size = 14),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


df_all <- read.csv("MIT_belowground_FT.csv", header = 1, sep = ",")
df_all$Treatment <- factor(df_all$Treatment, levels=c("Control", "Consortium B", "Consortium BP", "Fertilizers", "Pesticides", "Fertilizers-Pesticides"))

(p <- ggplot(df_all, aes(x = MIT_ric_sha_ratio, y = Belowground_biomass )) +
    #geom_text(aes(label = WHC), vjust = 1, hjust = -0.1, size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", color= "#08519c",alpha=0.1)+
    geom_ribbon(stat = "smooth", method = "lm", alpha=0.1, color="#08519c", linetype = 2)+
    #spearman/pearson method
    stat_cor(method = "pearson",label.y = 8,size = 5,label.x=-0.5,
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    stat_regline_equation( label.y = 7, label.x=-0.5,size = 5)+ 
    #geom_point(aes(fill = Treatment), size = 3, shape = 19, color="#3182bd", alpha = 0.8) +
    geom_point(aes(color = Treatment), size = 5, shape = 19, alpha = 0.8) +
    scale_color_manual(values = c("Control" = "#bfd2ec", "Consortium B" = "#a6c39d", "Consortium BP" = "#2a9c90","Fertilizers" = "#eab68e","Pesticides" = "#e2db47","Fertilizers-Pesticides" = "#ea8e94" )) +
    #facet_wrap(facets= "Time_point", nrow=1,ncol=3)+
    labs(x = "MIT (z score)", y = "Belowground biomass (g)", title = "") +
    mytheme)


ggsave("Figure1e.pdf", width = 12, height = 6, p)



