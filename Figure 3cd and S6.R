# PotatoMETAbiom Field trial experiment
# NL - heatmap-significant different genera among management - Week5

##Initiate libraries####
rm(list=ls())
library(dplyr)
library(ggplot2)
library(FSA)
library(rstatix)
library(tidyr)
library(tibble)
library(reshape2)
library(vegan)
library(ggplot2)

##Figure 3c heatmap illustrating significant differential genera between biological and chemical management for bacteria 
#bacterial genus table
bac <- read.csv("otu_bac_rhi_W5_t_genus.csv", header = 1, row.names = 1, sep = ",")

#Group by Genus and calculate the total number of OTUs for each Genus
Genus_col <- grep("Genus", colnames(bac))
Genus_otu <- aggregate(bac[, -Genus_col], by = list(bac[, Genus_col]), sum)

#Save file
write.csv(Genus_otu, file = "Genus_16S.csv", row.names = FALSE) #remove uncultured genus

bac_genus <- read.csv("Genus_16S.csv", header = 1, row.names = 1, sep = ",")

# metadata
sample_metadata <- read.csv("metadata_16S_rhi_W5.csv", header = 1, row.names = 1, sep = ",")

#treatments as factor
sample_metadata$Management <- as.factor(sample_metadata$Management)


##relative abundance select
otu_table <- as.data.frame(bac_genus)

otu_table <- as.data.frame(lapply(otu_table, as.numeric))

total_abundance <- colSums(otu_table)

relative_abundance <- bac_genus / 5700 #5700 is the cutoff of bacterial otu table

#Filter out those OTUs whose relative abundance in most samples is less than 0.1%
filtered_otu <- relative_abundance[rowSums(relative_abundance > 0.001) > 0, ]

filtered_otu <- filtered_otu %>%
  rownames_to_column("OTU") %>%
  pivot_longer(cols = -OTU, names_to = "SampleID", values_to = "Abundance")

#Merge sample metadata
otu_data <- filtered_otu %>%
  left_join(sample_metadata, by = "SampleID")

## Kruskal-Wallis test for management

kruskal_results <- otu_data %>%
  group_by(OTU) %>%
  summarise(p.value = kruskal.test(Abundance ~ Management)$p.value)

# significantly different genus
significant_otus <- kruskal_results %>%
  filter(p.value < 0.05) %>%
  pull(OTU)

# Dunn's test for pairwise comparisons, p adjust
dunn_results <- list()

for (otu in significant_otus) {
  subset_data <- subset(otu_data, OTU == otu)
  subset_data$Management <- as.factor(subset_data$Management) 
  dunn_test <- dunnTest(Abundance ~ Management, data = subset_data, method = "bh")#bonferroni
  dunn_results[[otu]] <- as.data.frame(dunn_test$res) 
}

# add genus onfoinfor
dunn_results <- lapply(names(dunn_results), function(otu) {
  result <- dunn_results[[otu]]
  result$OTU <- otu
  return(result)
})

# as dataframe
dunn_results_df <- bind_rows(dunn_results)

# Filter out rows with a P.unadj value less than 0.05
significant_results_man <- dunn_results_df %>%
  filter(P.adj < 0.05)

#keep comparison between biological -chemical
high_sig_man <- significant_results_man %>%
  filter(Comparison == "Biological - Chemical")%>% 
  filter(P.adj < 0.05)#%>%filter(abs(Z) > 4)
  
unique_importance_long <- unique(high_sig_man$OTU)
print(unique_importance_long)

significant_results_sorted_man <- high_sig_man[order(abs(high_sig_man$Z), decreasing = TRUE), ]

#write.csv(significant_results_sorted_man, "sig_features_man_bac_relabu0.01_genus.csv") #rownames

#select top genus otu table
df <- read.csv("sig_features_man_bac_relabu0.01_genus.csv", header = 1, row.names = 1, sep = ",")
bac_filtered <- bac_genus[rownames(bac_genus) %in% rownames(df), ]

write.csv(bac_filtered, "dif_man_bac_otu_table_relabu0.01_genus.csv") 

#use the dif genus table to plot the heatmap
library(ComplexHeatmap)
library(dplyr)

df <- read.csv("dif_man_bac_otu_table_relabu0.001_genus.csv", header = 1, row.names = 1, sep = ",") #as rownames

df_mean <- df %>%group_by(Treatment) %>%summarise(across(1:39, mean, na.rm = TRUE), .groups = "drop")

write.csv(df_mean, file = "top25_bac_0.001.csv") #Transformation for heatmap

#heatmap
mat <- read.csv("top25_bac_0.001.csv", header = 1, row.names = 1, sep = ",")

for (i in 1:nrow(mat)) mat[i, ] <- scale(log(unlist(mat[i, ])+0.001 , 2))#2 mean log2

mat <- as.matrix(mat)

#group
sample_group <- as.matrix(read.delim('group_man.txt', row.names = 1))

pdf("Figure 3c.pdf", width = 6, height = 6)

p3 <- Heatmap(
  mat, name = 'value', #functional abundance, Standardized\nvalue
  col = colorRampPalette(c("#7569a6", "#f5edf0", "#ffb800"))(100), 
  cluster_rows = T, cluster_columns = T,
  row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 7),  
  column_split = sample_group,  
  row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7))

p3
dev.off()



##Figure 3d heatmap illustrating significant differential genera between biological and chemical management for fungi 
#fun
# fun
fun <- read.csv("otu_fun_rhi_W5_t_genus.csv", header = 1, row.names = 1, sep = ",")

Genus_col <- grep("Genus", colnames(fun))
Genus_otu <- aggregate(fun[, -Genus_col], by = list(fun[, Genus_col]), sum)

write.csv(Genus_otu, file = "Genus_ITS.csv", row.names = FALSE) # remove uncultured genus

fun_genus <- read.csv("Genus_ITS.csv", header = 1, row.names = 1, sep = ",")

# metadata
sample_metadata <- read.csv("metadata_ITS_rhi_W5.csv", header = 1, row.names = 1, sep = ",")

#treatments as factor
sample_metadata$Management <- as.factor(sample_metadata$Management)


##relative abundance

otu_table <- as.data.frame(fun_genus)
otu_table <- as.data.frame(lapply(otu_table, as.numeric))

total_abundance <- colSums(otu_table)

relative_abundance <- fun_genus / 9800 #9800 is the cutoff of fungal otu table

filtered_otu <- relative_abundance[rowSums(relative_abundance > 0.001) > 0, ]

filtered_otu <- filtered_otu %>%
  rownames_to_column("OTU") %>%
  pivot_longer(cols = -OTU, names_to = "SampleID", values_to = "Abundance")

otu_data <- filtered_otu %>%
  left_join(sample_metadata, by = "SampleID")

# Kruskal-Wallis test for management

kruskal_results <- otu_data %>%
  group_by(OTU) %>%
  summarise(p.value = kruskal.test(Abundance ~ Management)$p.value)


significant_otus <- kruskal_results %>%
  filter(p.value < 0.05) %>%
  pull(OTU)

# Dunn
dunn_results <- list()

for (otu in significant_otus) {
  subset_data <- subset(otu_data, OTU == otu)
  subset_data$Management <- as.factor(subset_data$Management)
  dunn_test <- dunnTest(Abundance ~ Management, data = subset_data, method = "bh")#bonferroni
  dunn_results[[otu]] <- as.data.frame(dunn_test$res) 
}

dunn_results <- lapply(names(dunn_results), function(otu) {
  result <- dunn_results[[otu]]
  result$OTU <- otu
  return(result)
})

dunn_results_df <- bind_rows(dunn_results)

print(dunn_results_df)

write.csv(dunn_results_df, file = "dunn_results_man_fun_reab.csv", row.names = FALSE)

significant_results <- dunn_results_df %>%
  filter(P.adj < 0.05)


#keep comparison between biological -chemical
high_sig <- significant_results %>%
  filter(Comparison == "Biological - Chemical")%>% 
  filter(P.adj < 0.05)#%>% filter(abs(Z) > 3.5)
  
unique_importance_long <- unique(high_sig$OTU)
print(unique_importance_long)

significant_results_sorted <- high_sig[order(abs(high_sig$Z), decreasing = TRUE), ]

#write.csv(significant_results_sorted, "sig_features_man_fun_relabu_0.001_genus.csv") #change rownames as OTU list

#select top genus table
df <- read.csv("sig_features_man_fun_relabu_0.001_genus.csv", header = 1, row.names = 1, sep = ",")
fun_filtered <- fun_genus[rownames(fun_genus) %in% rownames(df), ]

write.csv(fun_filtered, "dif_man_fun_otu_table_relabu0.001_genus.csv") 


#FUN
#read table

df <- read.csv("dif_man_fun_otu_table_relabu0.001_genus.csv", header = 1, row.names = 1, sep = ",") #as rownames


df_mean <- df %>%group_by(Treatment) %>%summarise(across(1:27, mean, na.rm = TRUE), .groups = "drop")

#write.csv(df_mean, file = "top25_fun_0.001.csv")

mat <- read.csv("top25_fun_0.001.csv", header = 1, row.names = 1, sep = ",")

for (i in 1:nrow(mat)) mat[i, ] <- scale(log(unlist(mat[i, ])+0.001 , 2))#2 mean log2

mat <- as.matrix(mat)


sample_group <- as.matrix(read.delim('group_man.txt', row.names = 1))


pdf("Figure 3d", width = 6, height = 6)

p4<-Heatmap(
    mat, name = 'Relative abundance (log2-scaled)', #functional abundance, Standardized\nvalue
    col = colorRampPalette(c("#7569a6", "#f5edf0", "#ffb800"))(100), 
    cluster_rows = T, cluster_columns = T, 
    row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 7),  
    column_split = sample_group,  
    row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),  )

p4
dev.off()


##Figure 6, same code, but change the comparision among different levels of MIT

