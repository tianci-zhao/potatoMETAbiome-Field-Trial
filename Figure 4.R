# PotatoMETAbiom Field trial experiment
# NL - microbial commmunity network - Week5

##Initiate libraries####
rm(list=ls())

library(BiocManager)
library(edgeR)
library(indicspecies)
library(igraph)
library(Hmisc)
library(sciplot)
library(reshape2)
library(ggpmisc)
library(readxl)
source("CorrDF.R")


###Soil meta co-occurrence network creation and analysis and defining keystone OTUs under biological management
#Biological management####
###bacterial####
otu <- read.table("otu_table_16S_rhi_W5.csv",row.names=1,sep=";",header=T, blank.lines.skip=F,check.names=F)
otu <- as.matrix(otu)

#metadata,according different management
design_16S<- read_excel("metadata_16S_bio.xls",sheet = 1,col_names = T)

#match otu table to metadata
otu_16S_t <- t(otu)

otu_16S_bio_t <- otu_16S_t[match(design_16S$`Sampleid`,rownames(otu_16S_t)),]

otu_16S_bio <- t(otu_16S_bio_t)

#delete rows containing 0 that more than 4 columns (40 samples-4)
otu_16S_re <- otu_16S_bio[rowSums(otu_16S_bio == 0) <= 36, ]

####indicator species, find significant difference OTUs between groups####
edgeR_16S_t <- DGEList(counts=otu_16S_re, 
                       group=design_16S$Treatment)
#CPM normalization
otu_norm_16S_t <- cpm(edgeR_16S_t, normalized.lib.sizes=T, log=F)

#write.table(otu_norm_16S_t,"otu_norm_16S_FP.txt",sep="\t",row.names=T,col.names=T,quote=F)
#indicator species
indic_t_16S <- as.data.frame(t(otu_norm_16S_t))
indic_t_groups_16S <- design_16S$Treatment


set.seed(8046)
# Identify each group of indicator species
indicatorsp_t_16S <- multipatt(indic_t_16S,indic_t_groups_16S,func = "r.g",control=how(nperm=1000))
indic_t_df_16S <- indicatorsp_t_16S$sign

#select p.value < 0.05
Con_16S <- as.matrix(indic_t_df_16S[which(indic_t_df_16S$s.2 == 1 & indic_t_df_16S$p.value < 0.05),])
Eco_16S <- as.matrix(indic_t_df_16S[which(indic_t_df_16S$s.3 == 1 & indic_t_df_16S$p.value < 0.05),])

#merge data
t_r_values_16S <- rbind(Con_16S,Eco_16S)

#deleted redundant "s in colnames
colnames(t_r_values_16S)[1:2] <- gsub("s.","",colnames(t_r_values_16S)[1:2])

####iedgeR, find significant difference OTUs between groups####
#Get Group responsive rhizosphere soil bacteria OTUs with likelihood ratio testing in edgeR
head(design_16S)
model_matt_16S <- model.matrix(~Treatment, data=design_16S)

edgeR_16S_t_Group <- DGEList(counts=otu_16S_re, group=design_16S$Treatment)
edgeR_16S_t_Group <- calcNormFactors(edgeR_16S_t_Group)
dge_Group_16S <- estimateGLMRobustDisp(edgeR_16S_t_Group, design=model_matt_16S)
fit_Group_16S <- glmFit(dge_Group_16S, design=model_matt_16S)

#comparision
lrt_Group_16S <- glmLRT(fit_Group_16S, coef=2)
Group_t_16S <- topTags(lrt_Group_16S, n=Inf, p.value=0.05)
Group_t_16S <- Group_t_16S$table

#intersection of significant OTUs based on indicator species and edgeR
indic_edge_16S_t <- intersect(rownames(t_r_values_16S),rownames(Group_t_16S))


###fungi####
otu_its <- read.table("otu_table_ITS_rhi_w5.csv",row.names=1,sep=";",header=T, blank.lines.skip=F,check.names=F)
otu <- as.matrix(otu_its)

#design_its<- read_excel("metadata_its_bio.xls",sheet = 1,col_names = T)# same number of samples, can use same data with bacteria 
####match otu table to metadata
otu_its_t <- t(otu)

otu_its_bio_t <- otu_its_t[match(design_16S$`Sampleid`,rownames(otu_its_t)),]

otu_its_bio <- t(otu_its_bio_t)
#otu_its_bio <- otu_its_bio[, colSums(is.na(otu_its_bio)) == 0]

otu_its_re <- otu_its_bio[rowSums(otu_its_bio == 0) <= 36, ]

edgeR_its_t <- DGEList(counts=otu_its_re, 
                       group=design_16S$Treatment)

otu_norm_its_t <- cpm(edgeR_its_t, normalized.lib.sizes=T, log=F)

indic_t_its <- as.data.frame(t(otu_norm_its_t))
indic_t_groups_its <- design_16S$Treatment

set.seed(8046)

indicatorsp_t_its <- multipatt(indic_t_its,indic_t_groups_its,func = "r.g",control=how(nperm=1000))
indic_t_df_its <- indicatorsp_t_its$sign

Con_its <- as.matrix(indic_t_df_its[which(indic_t_df_its$s.2 == 1 & indic_t_df_its$p.value < 0.05),])
Eco_its <- as.matrix(indic_t_df_its[which(indic_t_df_its$s.3 == 1 & indic_t_df_its$p.value < 0.05),])

t_r_values_its <- rbind(Con_its,Eco_its)

colnames(t_r_values_its)[1:2] <- gsub("s.","",colnames(t_r_values_its)[1:2])

head(design_16S)
model_matt_its <- model.matrix(~Treatment, data=design_16S)

edgeR_its_t_Group <- DGEList(counts=otu_its_re, group=design_16S$Treatment)
edgeR_its_t_Group <- calcNormFactors(edgeR_its_t_Group)
dge_Group_its <- estimateGLMRobustDisp(edgeR_its_t_Group, design=model_matt_its)
fit_Group_its <- glmFit(dge_Group_its, design=model_matt_its)

lrt_Group_its <- glmLRT(fit_Group_its, coef=2)
Group_t_its <- topTags(lrt_Group_its, n=Inf, p.value=0.05)
Group_t_its <- Group_t_its$table

indic_edge_its_t <- intersect(rownames(t_r_values_its),rownames(Group_t_its))

###merge data####
## Combine OTU counts of both kingdoms together
otu_norm_soil_combine <- rbind(otu_norm_16S_t, otu_norm_its_t)

## Perform Spearman correlation of all OTU pairs
all_soil_cor <- rcorr(t(otu_norm_soil_combine), type=c("spearman"))

## Create data frame of co-occurring OTUs
all_cor_soil_df <- CorrDF(all_soil_cor$r, all_soil_cor$P)
all_cor_soil_df$padj <- p.adjust(all_cor_soil_df$p, method="none")

all_cor_soil_df_padj <- all_cor_soil_df[which(abs(all_cor_soil_df$cor) > 0.6),]
all_cor_soil_df_padj <- all_cor_soil_df_padj[which(all_cor_soil_df_padj$padj < 0.01),]

## Make node attribute table
indic_edge_soil_combine <- c(indic_edge_16S_t, indic_edge_its_t)

soil_r_values_combine <- rbind(t_r_values_16S, t_r_values_its)

nodeattrib_soil_combine <- data.frame(node=union(all_cor_soil_df_padj$from,all_cor_soil_df_padj$to))
nodeattrib_soil_combine$indicgroup <- 0

for (i in as.character(nodeattrib_soil_combine$node))
{
  if (i %in% indic_edge_soil_combine == TRUE)
  {nodeattrib_soil_combine[nodeattrib_soil_combine$node==i,"indicgroup"] <- paste(colnames(soil_r_values_combine)[which(soil_r_values_combine[i,1:2]==1)],collapse = "_")}
  else
  {nodeattrib_soil_combine[nodeattrib_soil_combine$node==i,"indicgroup"]<- "NA"}
}

rownames(nodeattrib_soil_combine) <- as.character(nodeattrib_soil_combine$node)

#all_soil_net <- graph_from_data_frame(all_cor_soil_df_padj,direct=F,vertices=nodeattrib_soil_combine)

all_cor_soil_df_padj$from.class <- sapply(all_cor_soil_df_padj[,1], function(x) nodeattrib_soil_combine[nodeattrib_soil_combine[,1] == x, 2])

all_cor_soil_df_padj$to.class <- sapply(all_cor_soil_df_padj[,2], function(x) nodeattrib_soil_combine[nodeattrib_soil_combine[,1] == x, 2])

#write.table(nodeattrib_soil_combine,"nodeattrib_combine_bio.txt",sep="\t",row.names=T,col.names=T,quote=F)

#igraph to achieve co-occurence network
all_soil_net <- graph_from_data_frame(all_cor_soil_df_padj,direct=F,vertices = nodeattrib_soil_combine)

t_ra_all <- apply(otu_norm_soil_combine,1,mean)
t_ra_all <- t_ra_all[V(all_soil_net)$name]

## Set node shape
V(all_soil_net)$shape <- V(all_soil_net)$name

#V(all_soil_net)$shape[V(all_soil_net)$shape %in% names(soil_all_keystone)] <- "star"
V(all_soil_net)$shape[V(all_soil_net)$shape %in% rownames(otu_norm_16S_t)] <- "circle"
V(all_soil_net)$shape[V(all_soil_net)$shape %in% rownames(otu_norm_its_t)] <- "triangle"
node_shape <- c("circle","triangle")

cs <- c("2","3")#2, 3 is the treatment number

unique(V(all_soil_net)$indicgroup) # get significant information

V(all_soil_net)$color <- V(all_soil_net)$indicgroup
V(all_soil_net)$color[!V(all_soil_net)$color %in% cs] <- "gray30"
V(all_soil_net)$color[V(all_soil_net)$color == "2"] <- "#35A585"
#V(all_soil_net)$color[V(all_soil_net)$color == "2_3"] <- "#4b69cd"
V(all_soil_net)$color[V(all_soil_net)$color == "3"] <- "#4baacd"
V(all_soil_net)$frame.color <- V(all_soil_net)$color

node_cols <- c("gray30","#35A585","#4baacd")

V(all_soil_net)$size <- degree(all_soil_net)*0.1+2 # size of node, 
E(all_soil_net)$color[E(all_soil_net)$cor >= 0.6] <- "black" # color of edge
E(all_soil_net)$color[E(all_soil_net)$cor <= -0.6] <- "grey" # color of edge
E(all_soil_net)$width <- abs(E(all_soil_net)$cor)*0.8 # width of edge


soil_nodes <- rownames(nodeattrib_soil_combine[nodeattrib_soil_combine$indicgroup %in% cs,])

soil_nodesizes <- as.numeric(V(all_soil_net)$size)


#Vectorize each type of significant inter-group difference OTUs for subsequent calculations
Con_nodes <- rownames(nodeattrib_soil_combine[nodeattrib_soil_combine$indicgroup=="2",])
#Con_Eco_nodes <- rownames(nodeattrib_soil_combine[nodeattrib_soil_combine$indicgroup=="2_3",])
Eco_nodes <- rownames(nodeattrib_soil_combine[nodeattrib_soil_combine$indicgroup=="3",])

cs_nodes_t_all <- c(Con_nodes,Eco_nodes)


#set layout, Fruchterman & Reingold
set.seed(8051)

#takes time, can save this file
coords_soil_all <- layout_(all_soil_net,with_fr(niter=9999, grid="nogrid"))
#write.table(coords_soil_all,"coords_all_bio.txt",sep="\t",row.names=F,col.names=F,quote=F)
#coords_soil_all <- as.matrix(read.table("coords_all_bio.txt"))

dimnames(coords_soil_all) <- NULL


### plot figure4####

set.seed(8051)
pdf(paste0("bio_BF.pdf"),width=7,height=5)

par(mfrow=c(1,1), mar=c(4,0,0,0))

plot(all_soil_net,vertex.label=NA,vertex.size=soil_nodesizes, layout=coords_soil_all,)
legend("bottom",legend=c("NS","Comsortium","EcoStyle"), #"Non-significant","Enriched in Comsortium"
       col=node_cols, bty="n", fill=node_cols, horiz=TRUE, ncol=1, y.intersp=-1, x.intersp = 0.5, yjust = 0)
legend("bottom",legend=c("Bacteria", "Fungi"), 
       pch= c(16, 17), bty="n", col="gray30",horiz=TRUE, ncol=1, y.intersp=-4, x.intersp = 0.5, yjust = 0)

dev.off()



#Chemical management####
#Same method with biological management (showe above), just replace the data and slightly adjust
design_16S<- read_excel("metadata_16S_che.xls",sheet = 1,col_names = T)
design_its<- read_excel("metadata_its_che.xls",sheet = 1,col_names = T) # the samples have different number for fungi and bacteria

#delete rows containing 0 that more than 4 columns (No.samples-4)
otu_16S_re <- otu_16S_bio[rowSums(otu_16S_bio == 0) <= 36, ]# 36 should be replaced



#Links information####
library(ggplot2)
library(data.table)
library(ggpubr)

### bar of proportion of different links,all_cor_soil_df_padj all links
fun_rows <- grepl("^OTU", all_cor_soil_df_padj$from) & 
  grepl("^OTU", all_cor_soil_df_padj$to) # filter fungi to fungi links
fun_df <- subset(all_cor_soil_df_padj, fun_rows)
count_fun <- c(sum(fun_df$cor > 0), sum(fun_df$cor < 0)) #positive and negative links
#count_fun <- nrow(fun_df)

bac_rows <- !grepl("^OTU", all_cor_soil_df_padj$from) & 
  !grepl("^OTU", all_cor_soil_df_padj$to)# filter bacteria to bacteria links
bac_df <- subset(all_cor_soil_df_padj, bac_rows)
count_bac <- c(sum(bac_df$cor > 0), sum(bac_df$cor < 0))

fb_rows <- grepl("^OTU", all_cor_soil_df_padj$from) & 
  !grepl("^OTU", all_cor_soil_df_padj$to)# filter fungi to bacteria links
fb_df <- subset(all_cor_soil_df_padj, fb_rows)
count_fb <- c(sum(fb_df$cor > 0), sum(fb_df$cor < 0))


pro_adge_bac <- count_bac/381
pro_adge_fun <- count_fun/682
pro_adge_fb <- count_fb/700

adge_bac <- data.table(pro_adge_bac, labels=c("Positive", "Negative"))
adge_bac$mic <- "bac" #add a column, whcih mean from bac to bac

adge_fun<- data.table(pro_adge_fun, labels=c("Positive", "Negative"))
adge_fun$mic <- "fun" 

adge_fb<- data.table(pro_adge_fb, labels=c("Positive", "Negative"))
adge_fb$mic <- "fun to bac"

adge_plot <- rbind(adge_bac, adge_fun,adge_fb,use.names=FALSE)

#write.csv(adge_plot, file = "adge_bio.csv") #add link number info
#same for chemical management, merge adge_bio and adge_che as adge_all

#figure
adge_plot<- data.table::fread("adge_all.csv")

mic_order <- c("bac-bacteria","fungi-fungi", "fungi-bacteria")
adge_plot$mic <- factor(adge_plot$mic, levels = mic_order)


p <- ggplot(data=adge_plot, aes(x=mic, y=link_percent,fill=labels)) +
  geom_bar(width=0.6, stat="identity") +
  xlab("") +
  ylab("Proportion of positive and negative correlations") +
  #ggtitle("Proportion of Values in Column 3") +
  scale_fill_manual(values=c("#9ECAE1", "#F8766D")) +
  scale_y_continuous(limits = c(0, 0.5), labels = scales::percent)+
  facet_grid(cols=vars(management)) +
  theme_pubr()+
  theme(axis.title.y = element_text(size = 12, color = "black"),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, color = "black"),
        legend.position = "top", legend.title=element_blank(),legend.key.width=unit(1, "lines"), 
        legend.key.height=unit(1, "lines"),
        legend.text = element_text(size = 12), legend.justification = "left")
p

#ggsave("Proportion of positive and negative correlations_all.pdf", p, width = 9, height = 6, units = "cm", scale = 1.5)


