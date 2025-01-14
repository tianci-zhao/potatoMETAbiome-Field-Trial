# PotatoMETAbiom Field trial experiment
# NL - heatmap - Treatment-responsive ASVs/OTUs in biological/chemical management

##Initiate libraries####
#library
library(ComplexHeatmap)
library(dplyr)


##Figure S6 Treatment-responsive ASVs/OTUs in biological management###
#got indicator species table
#bio indicator otu table
df <- read.csv("indicgroup_bio.csv", header = 1, row.names = 1, sep = ",") #as rownames
otu <- read.csv("otu_bac_fun_total_119.csv", header = 1, row.names = 1, sep = ",") #as rownames

sig_rows <- rownames(df)

otu_filtered <- otu[sig_rows,]

print(otu_filtered)

write.csv(otu_filtered, "indicgroup_sig_otu_bio.csv") #transpose

#read table

df <- read.csv("indicgroup_sig_otu_bio.csv", header = 1, row.names = 1, sep = ",") #as rownames

df_mean <- df %>%group_by(Treatment) %>%summarise(across(1:131, mean, na.rm = TRUE), .groups = "drop")

#write.csv(df_mean, file = "bio_indicator_otu_table_tre_mean.csv") #transpose for heatmap

#heatmap
mat <- read.csv("bio_indicator_otu_table_tre_mean_family_rm_unidentified.csv", header = 1, row.names = 1, sep = ",")

mat <- mat[,1:6]

for (i in 1:nrow(mat)) mat[i, ] <- scale(log(unlist(mat[i, ])+0.001 , 2))# 2 mean log2

mat <- as.matrix(mat)

sample_group <- as.matrix(read.delim('group_tre.txt', row.names = 1))

#heatmap

pdf("Figure S6-left.pdf", width = 8, height = 10)

ht<- Heatmap(
    mat, name = 'value', #functional abundance, Standardized\nvalue
    col = colorRampPalette(c("#7569a6", "#f5edf0", "#ffb800"))(100), 
    cluster_rows = T, cluster_columns = T,  
    row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
    column_split = sample_group,
    row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),)
dev.off()



#Extract row names and sort them in cluster order
ordered_rows <- rownames(mat)[row_order(ht)]
print(ordered_rows)
write.csv(ordered_rows, file = "ordered_rows_bio.csv")


#bio node lollipop
library(ggplot2)
library(ggbreak)

df <- read.csv("bio_family_degree.csv", header = TRUE, sep = ",")

df$Family <- reorder(df$Family, df$order)

# lollipop
p1 <- ggplot(df, aes(x = degree, y = Family, color = Taxa)) +
  geom_segment(aes(x = 0, xend = degree, y = Family, yend = Family)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#3e413a","#aaaea5")) +
  theme_minimal(base_size = 12) +theme(
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
  ) +
  labs(
    title = "      ",
    x = " ",
    y = " "
  ) #+scale_x_break(c(10,20), scales = 1)#X-axis cutoff: cut off between 1000 and 80000

print(p1)

ggsave("Figure S6-right.pdf", plot = p1, width = 5, height = 10)



#Figure S7 Treatment-responsive ASVs/OTUs in chemical management

#che indicator otu table
df <- read.csv("indicgroup_che.csv", header = 1, row.names = 1, sep = ",") #as rownames

sig_rows <- rownames(df)

otu_filtered <- otu[sig_rows,]

print(otu_filtered)

write.csv(otu_filtered, "indicgroup_sig_otu_che.csv")#transpose

df <- read.csv("indicgroup_sig_otu_che.csv", header = 1, row.names = 1, sep = ",") #as rownames

df_mean <- df %>%group_by(Treatment) %>%summarise(across(1:73, mean, na.rm = TRUE), .groups = "drop")

#write.csv(df_mean, file = "che_indicator_otu_table_tre_mean.csv") ##transpose for heatmap

mat <- read.csv("che_indicator_otu_table_tre_mean_family_rm_unidentified.csv", header = 1, row.names = 1, sep = ",")

mat <- mat[,1:6]

for (i in 1:nrow(mat)) mat[i, ] <- scale(log(unlist(mat[i, ])+0.001 , 2))#2 mean log2

mat <- as.matrix(mat)

sample_group <- as.matrix(read.delim('group_tre.txt', row.names = 1))




#heatmap
pdf("Figure S7-left.pdf", width = 8, height = 10)

ht<-Heatmap(
    mat, name = 'value', #functional abundance, Standardized\nvalue
    col = colorRampPalette(c("#7569a6", "#f5edf0", "#ffb800"))(100),  
    cluster_rows = T, cluster_columns = T,  
    row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),  
    column_split = sample_group, 
    row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),)
dev.off()

ordered_rows <- rownames(mat)[row_order(ht)]
print(ordered_rows)
write.csv(ordered_rows, file = "ordered_rows_che.csv")


#che node lollipop
library(ggplot2)
library(ggbreak)

df <- read.csv("che_family_degree.csv", header = TRUE, sep = ",")

df$Family <- reorder(df$Family, df$order)

# lollipop
p1 <- ggplot(df, aes(x = degree, y = Family, color = Taxa)) +
  geom_segment(aes(x = 0, xend = degree, y = Family, yend = Family)) +
  geom_point(size = 3) +）
  scale_color_manual(values = c("#3e413a","#aaaea5")) +
  theme_minimal(base_size = 12) +theme(
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
  ) +
  labs(
    title = "      ",
    x = " ",
    y = " "
  ) #+scale_x_break(c(6000,10000), scales = 1)


print(p1)

ggsave("Figure S7-right.pdf", plot = p1, width = 5, height = 10)




