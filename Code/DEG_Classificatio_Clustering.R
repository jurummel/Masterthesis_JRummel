library(tidyverse)
library(ComplexHeatmap)
library(dendextend)

#-------------------get Cluster--------------------

trans_cts <- read.csv("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Clustering_genes/DEGS_Clustering_csv.csv", sep = ";", dec=".")
sample_info <- read.csv("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Clustering_genes/metadata.csv", sep = ";", dec=".")

trans_cts_mean_pre <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = F2_001CC:BL6_18, names_to = "sample", values_to = "cts") %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>% 
  # for each gene
  group_by(X) %>% 
  # scale the cts column (z-score)
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(X, strain, condition) %>%
  # calculate the mean z-score
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()

# get df for clustering 
trans_cts_mean_pre$sample_new = paste(trans_cts_mean_pre$strain, trans_cts_mean_pre$condition, sep="_")

# get df for clustering
trans_cts_mean <- trans_cts_mean_pre[c("X", "sample_new", "mean_cts_scaled")]
trans_cts_mean <- pivot_wider(trans_cts_mean, id_cols=(c(X)), values_from = mean_cts_scaled,  names_from=sample_new)
trans_cts_mean_pre <- trans_cts_mean_pre[c("X", "strain", "condition", "mean_cts_scaled","nrep")]

# Create a matrix
hclust_matrix <- trans_cts_mean%>% 
  select(-X) %>% 
  as.matrix()


# assign rownames
rownames(hclust_matrix) <- trans_cts_mean$X 

# calculate distance (euclidean)
gene_dist <- dist(hclust_matrix)

# hierarchical cluster
gene_hclust <- hclust(gene_dist, method = "complete")


# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 1.8, col = "brown", lwd = 2)

# save cluster information
gene_cluster <- cutree(gene_hclust, k = 10) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(X = name, cluster = value)

# test
head(gene_cluster)

# join cluster info + expression info 
trans_cts_cluster <- trans_cts_mean_pre %>% 
  inner_join(gene_cluster, by = "X")

# test
head(trans_cts_cluster)

# Duplicate data frame
data_new2 <- trans_cts_cluster 

# change cluster to match with heatmap
#1=1, 2=2, 30=5, 40=8, 50=10, 60=4, 70=6, 80=9, 90=7, 100=3
data_new2["cluster"] <- sapply(data_new2["cluster"], function(x) replace(x,  x %in% c(100), 10)) 
data_new2["cluster"] <- sapply(data_new2["cluster"], function(x) replace(x,  x %in% c(3), 10))  
data_new2["cluster"] <- sapply(data_new2["cluster"],function(x) replace(x,  x %in% c(99), 3))  
trans_cts_cluster <- data_new2

# print diagram for every cluster
trans_cts_cluster %>% 
  ggplot(aes(factor(condition, level=c("untreated", "treated")), mean_cts_scaled)) +
  geom_line(aes(group = X)) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(strain), cols = vars(cluster)) +
  xlab("Condition") + ylab("Z-score")

# get complex heatmap 
column_new=c("cc_untreated", "cc_treated", "bl6_untreated", "bl6_treated")
dend = as.dendrogram(gene_hclust)
dend = color_branches(dend, k = 10)
ht = Heatmap(name="Z-score",hclust_matrix, show_row_names = FALSE, cluster_rows=dend, row_split = 10, 
             border = TRUE, row_dend_reorder = FALSE, column_order = column_new, cluster_columns = FALSE, 
             row_gap = unit(2, "mm"), column_names_side = "bottom", column_names_rot = 45, row_title = "Genes",
             left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:1),
                                                              labels = c("1", "2", "3","4", "5", "6","7", "8", "9","10"), 
                                                              labels_gp = gpar(col = "white", fontsize = 8))))
# define values for later use
ht = draw(ht)
orders = row_order(ht)
list_components()

# decorate heatmap
for (i in 1:10){
  decorate_heatmap_body("Z-score", {
    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 4))
  },slice = c(i))
}

# loop to extract genes for each cluster.
for (i in 1:length(orders)){
  if (i == 1) {
    clu <- t(t(row.names(hclust_matrix[orders[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
    } else {
      clu <- t(t(row.names(hclust_matrix[orders[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
      }
  }

# define number of genes in each cluster
out = data.frame(out)
cluster_1_n = count(out[out$Cluster=="cluster1",]) #167
cluster_2_n = count(out[out$Cluster=="cluster2",]) #625
cluster_3_n = count(out[out$Cluster=="cluster3",]) #112
cluster_4_n = count(out[out$Cluster=="cluster4",]) #78
cluster_5_n = count(out[out$Cluster=="cluster5",]) #39
cluster_6_n = count(out[out$Cluster=="cluster6",]) #68
cluster_7_n = count(out[out$Cluster=="cluster7",]) #73
cluster_8_n = count(out[out$Cluster=="cluster8",]) #21
cluster_9_n = count(out[out$Cluster=="cluster9",]) #90
cluster_10_n = count(out[out$Cluster=="cluster10",]) #363

# as vector
cluster_num_n=c(cluster_1_n,cluster_2_n,cluster_3_n,cluster_4_n,cluster_5_n,cluster_6_n,cluster_7_n,cluster_8_n,cluster_9_n,cluster_10_n)
cluster_num_n

# compare cluster from heatmap to cluster of hierarchical clustering  
#1=1, 2=2, 3=5, 4=8, 5=10, 6=4, 7=6, 8=9, 9=7, 10=3
cluster_1 = count(trans_cts_cluster[trans_cts_cluster$cluster==1,])/4 #167
cluster_2 = count(trans_cts_cluster[trans_cts_cluster$cluster==2,])/4 #625
cluster_3 = count(trans_cts_cluster[trans_cts_cluster$cluster==3,])/4 #363
cluster_4 = count(trans_cts_cluster[trans_cts_cluster$cluster==4,])/4 #68
cluster_5 = count(trans_cts_cluster[trans_cts_cluster$cluster==5,])/4 #112
cluster_6 = count(trans_cts_cluster[trans_cts_cluster$cluster==6,])/4 #73
cluster_7 = count(trans_cts_cluster[trans_cts_cluster$cluster==7,])/4 #90
cluster_8 = count(trans_cts_cluster[trans_cts_cluster$cluster==8,])/4 #78
cluster_9 = count(trans_cts_cluster[trans_cts_cluster$cluster==9,])/4 #21
cluster_10 = count(trans_cts_cluster[trans_cts_cluster$cluster==10,])/4 #39

# as vector
cluster_num=c(cluster_1,cluster_2,cluster_3,cluster_4,cluster_5,cluster_6,cluster_7,cluster_8,cluster_9,cluster_10)
cluster_num

# check if identical
out_1 <- lapply(list(out$GeneID[out$Cluster=="cluster1"]), sort)
trans_1 <- lapply(list(unique(trans_cts_cluster$X[trans_cts_cluster$cluster==1])), sort)
identical(out_1, trans_1)

# get genes of each cluster
genes_cluster_1 <- lapply(list(out$GeneID[out$Cluster=="cluster1"]), sort)
genes_cluster_2 <- lapply(list(out$GeneID[out$Cluster=="cluster2"]), sort)
genes_cluster_3 <- lapply(list(out$GeneID[out$Cluster=="cluster3"]), sort)
genes_cluster_4 <- lapply(list(out$GeneID[out$Cluster=="cluster4"]), sort)
genes_cluster_5 <- lapply(list(out$GeneID[out$Cluster=="cluster5"]), sort)
genes_cluster_6 <- lapply(list(out$GeneID[out$Cluster=="cluster6"]), sort)
genes_cluster_7 <- lapply(list(out$GeneID[out$Cluster=="cluster7"]), sort)
genes_cluster_8 <- lapply(list(out$GeneID[out$Cluster=="cluster8"]), sort)
genes_cluster_9 <- lapply(list(out$GeneID[out$Cluster=="cluster9"]), sort)
genes_cluster_10 <- lapply(list(out$GeneID[out$Cluster=="cluster10"]), sort)

#-------------gprofiler2-------------------
library(gprofiler2)

# perform enrichment analysis for all cluster 
for (i in 1:10){
  nam <- paste("enrichment_cluster_", i, sep ="")
  variable <- lapply(list(out$GeneID[out$Cluster==paste("cluster",i, sep = "")]), sort)
  assign(nam, gost(query = c(variable),organism = "mmusculus", significant=TRUE, 
                               correction_method = "gSCS"))
}

# get # of all significant terms for all cluster
count_rows = list(nrow(enrichment_cluster_1$result), nrow(enrichment_cluster_2$result), nrow(enrichment_cluster_3$result), 
                  nrow(enrichment_cluster_4$result), nrow(enrichment_cluster_5$result), nrow(enrichment_cluster_6$result),
                  nrow(enrichment_cluster_7$result), nrow(enrichment_cluster_8$result), nrow(enrichment_cluster_9$result),
                  nrow(enrichment_cluster_10$result))

# plot overview of significant terms and save output to PNG
p <- gostplot(enrichment_cluster_1, capped = TRUE, interactive = TRUE)
publish_gostplot(p, width = NA, height = NA, filename ="/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Clustering_genes/cluster_10_enrichment.png")

# create data frame with the last column (parents) NULL and write to CSV
cluster_results <- data.frame(enrichment_cluster_10$result)
cluster_results[,c(14)] <- NULL
write.csv(cluster_results,
          "/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Clustering_genes/cluster_10_enrichment.csv")

# highlight significant terms of plor
pp <- publish_gostplot(p, highlight_terms = c("GO:0006631", "GO:0044255", "KEGG:00062", "KEGG:01040", "WP:WP4351", "WP:WP4350"), 
                       width = NA, height = NA, filename = NULL )


#-------------get Gene names-------------------

# read DEGs of BL6, CC and Overlap
DEG_BL6 <- read.csv("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/GitHub/Results/mRNA-seq/DEG/DEG_BL6_sum.csv", sep = ",", dec=".")
DEG_CC <- read.csv("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/GitHub/Results/mRNA-seq/DEG/DEG_CC_sum.csv", sep = ",", dec=".")
DEG_Overlapp <- merge(DEG_BL6, DEG_CC, by="X")

# subset BL6 and CC w/o overlapping genes (redundancy)
DEG_BL6_sub <- subset(DEG_BL6, !(X %in% c(DEG_Overlapp$X)))
DEG_CC_sub <- subset(DEG_CC, !(X %in% c(DEG_Overlapp$X)))

# change gene IDs to gene name with the ENSEMBL data base
mart <- useMart('ENSEMBL_MART_ENSEMBL', host = 'https://useast.ensembl.org')
mart <- useDataset('mmusculus_gene_ensembl', mart)
G_list_BL6 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol","gene_biotype"),values=c(DEG_BL6_sub$X),mart= mart)
G_list_CC <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol","gene_biotype"),values=c(DEG_CC_sub$X),mart= mart)
G_list_Lap <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol","gene_biotype"),values=c(DEG_Overlapp$X),mart= mart)

# safe gene names to file
fileBL6<-file("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/GitHub/Results/mRNA-seq/DEG/deg_BL6.txt")
writeLines(c(G_list_BL6$ensembl_gene_id), fileBL6)
close(fileBL6)

# safe gene names to file
fileCC<-file("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/GitHub/Results/mRNA-seq/DEG/deg_CC.txt")
writeLines(c(G_list_CC$ensembl_gene_id), fileCC)
close(fileCC)

# safe gene names to file
fileLAP<-file("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/GitHub/Results/mRNA-seq/DEG/deg_LAP.txt")
writeLines(c(G_list_Lap$ensembl_gene_id), fileLAP)
close(fileLAP)

# safe gene names to file
fileLAP<-file("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/GitHub/Results/mRNA-seq/DEG/degs_name.txt")
writeLines(c(c(G_list_BL6$mgi_symbol),c(G_list_CC$mgi_symbol), c(G_list_Lap$mgi_symbol)), fileLAP)
close(fileLAP)

# get gene names for different cluster
numerator <- 1
for (i in c(genes_cluster_1,genes_cluster_2,genes_cluster_3,genes_cluster_4,genes_cluster_5,genes_cluster_6,genes_cluster_7,genes_cluster_8,genes_cluster_9,genes_cluster_10)){
  nam <- paste("genes_cluster_",numerator, "_names", sep="")
  variable <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol","gene_biotype"),values=c(i),mart= mart)
  variable$mgi_symbol
  assign(nam, variable$mgi_symbol)
  numerator = numerator + 1
}

# check in which cluster the gene is present
sapply(list.files("/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/GitHub/Results/mRNA-seq/Test", full.names=TRUE), FUN=function(x){
  grep("Cdkn1a", readLines(x))
})
