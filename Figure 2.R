# PotatoMETAbiom Field trial experiment
# NL - microbial community composition - Week5

##Initiate libraries####
rm(list=ls())

library(tidyverse)
library(vegan)
library(ape) 
library(doBy)
library(reshape2)
library(patchwork)

mytheme <- 
  theme_bw() +
  theme(text = element_text(size=10),
        plot.title = element_text(size=11, face="bold"),
        strip.background = element_blank(),
        strip.text=element_text(face="bold", size = 10),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#Figure 2 The effect of treatment and cultivar on the composition of microbial communities in rhizosphere soil####
###Bacteria_import datasets####
# feature table
otu_table <- read.csv("otu_table_16S_rhi_W5.csv",header=T,row.names=1,sep = ";")

#metadata
sample_info <- read.csv("metadata_16S_rhi_W5.csv",header=T,row.names=1, sep = ",")

###beta-diversity####
all.dist <- vegdist(t(otu_table), method = "bray", binary = F)
all.pcoa <- cmdscale(all.dist, k=3, eig = T)
all_pcoa_points <- as.data.frame(all.pcoa$points)
sum_eig <- sum(all.pcoa$eig)
eig_percent <- round(all.pcoa$eig/sum_eig*100, 1)

colnames(all_pcoa_points) <- paste0("PCoA", 1:3)

## combine pcoa and sample info
all_pcoa_result <- cbind(all_pcoa_points, sample_info)
head(all_pcoa_result)

### adonis statistic####
set.seed(1)
all.div <- adonis2(all.dist ~ Treatment, data = all_pcoa_result, permutations = 999, method = "bray")
all.div # R2=0.1561,p=0.001 ***

all.div <- adonis2(all.dist ~ Variety, data = all_pcoa_result, permutations = 999, method = "bray")
all.div # ns R2=0.05133 , p= 0.169

## legend order
all_pcoa_result$Variety <- factor(all_pcoa_result$Variety,
                                  levels=c("ATOL", "DESIREE", "JELLY",  "KRAB", "PASJA_POMORSKA", "RUDAWA", "SALTO"),
                                  labels = c("ATOL", "DESIREE", "JELLY", "KRAB", "PASJA POMORSKA", "RUDAWA", "SALTO"))
all_pcoa_result$Treatment <- factor(all_pcoa_result$Treatment, 
                                    levels =c("1", "2", "3", "4", "5", "6"), 
                                    labels = c("Control", "Bacteria", "Bacteria-Protists", "Fertilizers", "Pesticides",  "Fertilizers-Pesticides"))

#write.csv(all_pcoa_result, file = "pcoa_rhi_bac_W5_result.csv")

Bac1 <- ggplot(all_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Treatment, group=Treatment)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep = ""), 
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
  geom_point(size=4, alpha=1.0) +
  #geom_encircle(aes(fill=Treatment),alpha=0.4, show.legend = T) +
  annotate("text", x = -0.06, y = 0.2, label = "")+
  scale_color_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  #scale_shape_manual(values = c(0, 1, 2, 3,7,8,10)) +
  mytheme + stat_ellipse(aes(group = Treatment, color = Treatment, fill = Treatment), level = 0.8, geom = "polygon", alpha = 0.1)
Bac1

Bac2 <- ggplot(all_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Variety, group=Variety,fill=Variety)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep = ""), 
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
  geom_point(size=4, alpha=1.0) +
  annotate("text", x = -0.06, y = 0.2, label = "")+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  scale_shape_manual(values = c(0, 1, 2, 3,7,8,10)) +
  mytheme+
  stat_ellipse(aes(group = Variety, color = Variety, fill =Variety), level = 0.8, geom = "polygon", alpha = 0.1)
Bac2


###Fungi_import datasets####
# feature table
otu_table <- read.csv("otu_table_ITS_rhi_W5.csv",header=T,row.names=1,sep = ";")

#metadata
sample_info <- read.csv("metadata_ITS_rhi_W5.csv",header=T,row.names=1, sep = ",")

###beta-diversity####
all.dist <- vegdist(t(otu_table), method = "bray", binary = F)
all.pcoa <- cmdscale(all.dist, k=3, eig = T)
all_pcoa_points <- as.data.frame(all.pcoa$points)
sum_eig <- sum(all.pcoa$eig)
eig_percent <- round(all.pcoa$eig/sum_eig*100, 1)

colnames(all_pcoa_points) <- paste0("PCoA", 1:3)

## combine pcoa and sample info
all_pcoa_result <- cbind(all_pcoa_points, sample_info)
head(all_pcoa_result)

all_pcoa_result$Treatment <- factor(all_pcoa_result$Treatment, 
                                    levels =c("1", "2", "3", "4", "5", "6"), 
                                    labels = c("Control", "Bacteria", "Bacteria-Protists", "Fertilizers", "Pesticides",  "Fertilizers-Pesticides"))

### adonis statistic####
set.seed(1)
all.div <- adonis2(all.dist ~ Treatment, data = all_pcoa_result, permutations = 999, method = "bray")
all.div # R2=0.24568,p=0.001 ***

all.div <- adonis2(all.dist ~ Variety, data = all_pcoa_result, permutations = 999, method = "bray")
all.div # ns R2=0.03991, p= 0.978

## legend order
all_pcoa_result$Variety <- factor(all_pcoa_result$Variety,
                                  levels=c("ATOL", "DESIREE", "JELLY",  "KRAB", "PASJA_POMORSKA", "RUDAWA", "SALTO"),
                                  labels = c("ATOL", "DESIREE", "JELLY", "KRAB", "PASJA POMORSKA", "RUDAWA", "SALTO"))

all_pcoa_result$Treatment <- factor(all_pcoa_result$Treatment, levels=c("Control","Bacteria", "Bacteria-Protists", "Fertilizers", "Pesticides", "Fertilizers-Pesticides"))

Fun1<-ggplot(all_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Treatment, group=Treatment)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep = ""), 
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
  geom_point(size=4, alpha=1.0) +
  #geom_encircle(aes(fill=Treatment),alpha=0.4, show.legend = T) +
  annotate("text", x = -0.06, y = 0.2, label = "")+
  scale_color_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94")) +
  scale_fill_manual(values=c("#bfd2ec","#a6c39d","#2a9c90","#eab68e","#e2db47","#ea8e94"))+
  #scale_shape_manual(values = c(0, 1, 2, 3,7,8,10)) +
  stat_ellipse(aes(group = Treatment, color = Treatment, fill = Treatment), level = 0.8, geom = "polygon", alpha = 0.1)+
  mytheme
Fun1


Fun2 <- ggplot(all_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Variety, group=Variety,fill=Variety)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep = ""), 
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
  geom_point(size=4, alpha=1.0) +
  annotate("text", x = -0.06, y = 0.2, label = "")+
  scale_color_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969")) +
  scale_fill_manual(values=c("#7799bc","#f6e0a5","#c98bb5","#a4757d","#e5977d","#8580a3","#535969"))+
  scale_shape_manual(values = c(0, 1, 2, 3,7,8,10)) +
  mytheme+
  stat_ellipse(aes(group = Variety, color = Variety, fill =Variety), level = 0.8, geom = "polygon", alpha = 0.1)
Fun2

## merge figures
Fig2 <- Bac1 + Fun1 + Bac2 + Fun2 + plot_layout(nrow = 2, guides= "collect")&
  theme(legend.position = 'right')

Fig2

ggsave("Figure2.pdf", width = 15, height = 12, units = "cm", Fig2, scale = 1.6)
