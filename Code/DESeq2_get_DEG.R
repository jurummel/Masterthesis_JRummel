library(stringr)
library(DESeq2)
library(VennDiagram)
library(VennDetail)
library(RColorBrewer)
library(apeglm)
library(ggfortify)

# data from snakepipes, RNA-seq Analysis, feature counts, gene expression
feature_counts = "/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/featureCounts_all.csv"  
data_feature_counts = read.csv(feature_counts, header = TRUE, sep = ";", dec=".", row.names = 1)
df_feature_counts = as.data.frame(data_feature_counts)
column_names_vector = c(colnames(df_feature_counts))

#metadata for sum 
meta_data_BL6_gen1 = "/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/meta_data_BL6_sum.csv"
data_meta_data_BL6_gen1 = read.csv(meta_data_BL6_gen1, header = TRUE, sep = ";", dec=".", row.names = 1)
meta_data_CC_gen1 = "/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/meta_data_CC_sum.csv"
data_meta_data_CC_gen1 = read.csv(meta_data_CC_gen1, header = TRUE, sep = ";", dec=".", row.names = 1)
meta_data_XB_gen1 = "/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/meta_data_XB_sum.csv"
data_meta_data_XB_gen1 = read.csv(meta_data_XB_gen1, header = TRUE, sep = ";", dec=".", row.names = 1)
meta_data_YB_gen1 = "/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/meta_data_YB_sum.csv"
data_meta_data_YB_gen1 = read.csv(meta_data_YB_gen1, header = TRUE, sep = ";", dec=".", row.names = 1)

# get table for different species
# C57Bl6/J animals bought from Janvier
df_fc_BL6_gen1 = df_feature_counts[, str_subset(column_names_vector, "^BL6\\w*")]
meta_BL6_gen1 = as.data.frame(data_meta_data_BL6_gen1)

# add columns with sum of genome1, genome2 and all 
df_fc_BL6_gen1$SUM_BL6_11 <- rowSums(df_fc_BL6_gen1[ , c(1,2,3)], na.rm=TRUE)
df_fc_BL6_gen1$SUM_BL6_12 <- rowSums(df_fc_BL6_gen1[ , c(4,5,6)], na.rm=TRUE)
df_fc_BL6_gen1$SUM_BL6_13 <- rowSums(df_fc_BL6_gen1[ , c(7,8,9)], na.rm=TRUE)
df_fc_BL6_gen1$SUM_BL6_14 <- rowSums(df_fc_BL6_gen1[ , c(10,11,12)], na.rm=TRUE)
df_fc_BL6_gen1$SUM_BL6_15 <- rowSums(df_fc_BL6_gen1[ , c(13,14,15)], na.rm=TRUE)
df_fc_BL6_gen1$SUM_BL6_16 <- rowSums(df_fc_BL6_gen1[ , c(16,17,18)], na.rm=TRUE)
df_fc_BL6_gen1$SUM_BL6_17 <- rowSums(df_fc_BL6_gen1[ , c(19,20,21)], na.rm=TRUE)
df_fc_BL6_gen1$SUM_BL6_18 <- rowSums(df_fc_BL6_gen1[ , c(22,23,24)], na.rm=TRUE)

# create table for DEG Analysis
new_cols_BL6 = c(colnames(df_fc_BL6_gen1))
df_fc_BL6_sum = df_fc_BL6_gen1[, str_subset(new_cols_BL6, "^SUM\\w*")]
write.csv(df_fc_BL6_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/BL6_SUM_FC.csv", row.names = TRUE)

# F2 Generation, Cast-EiJ father, Cast-EiJ mother
df_fc_CC_gen1 = df_feature_counts[, str_subset(column_names_vector, "^F2\\w*CC")]
meta_CC_gen1 = as.data.frame(data_meta_data_CC_gen1)

# add columns with sum of genome1, genome2 and all 
df_fc_CC_gen1$SUM_F2_001CC <- rowSums(df_fc_CC_gen1[ , c(1,2,3)], na.rm=TRUE)
df_fc_CC_gen1$SUM_F2_002CC <- rowSums(df_fc_CC_gen1[ , c(4,5,6)], na.rm=TRUE)
df_fc_CC_gen1$SUM_F2_003CC <- rowSums(df_fc_CC_gen1[ , c(7,8,9)], na.rm=TRUE)
df_fc_CC_gen1$SUM_F2_004CC <- rowSums(df_fc_CC_gen1[ , c(10,11,12)], na.rm=TRUE)
df_fc_CC_gen1$SUM_F2_006CC <- rowSums(df_fc_CC_gen1[ , c(13,14,15)], na.rm=TRUE)
df_fc_CC_gen1$SUM_F2_007CC <- rowSums(df_fc_CC_gen1[ , c(16,17,18)], na.rm=TRUE)

# create table for DEG Analysis
new_cols_CC = c(colnames(df_fc_CC_gen1))
df_fc_CC_sum = df_fc_CC_gen1[, str_subset(new_cols_CC, "^SUM\\w*")]
head(df_fc_CC_sum)
write.csv(df_fc_CC_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/CC_SUM_FC.csv", row.names = TRUE)

# F2 Generation, Cast-EiJ mother, Bl6 father 
df_fc_YB_gen1 = df_feature_counts[, str_subset(column_names_vector, "^F2\\w*YB")]
meta_YB_gen1 = as.data.frame(data_meta_data_YB_gen1)

# add columns with sum of genome1, genome2 and all 
df_fc_YB_gen1$SUM_F2_008YB <- rowSums(df_fc_YB_gen1[ , c(1,2,3)], na.rm=TRUE)
df_fc_YB_gen1$SUM_F2_009YB <- rowSums(df_fc_YB_gen1[ , c(4,5,6)], na.rm=TRUE)
df_fc_YB_gen1$SUM_F2_010YB <- rowSums(df_fc_YB_gen1[ , c(7,8,9)], na.rm=TRUE)
df_fc_YB_gen1$SUM_F2_012YB <- rowSums(df_fc_YB_gen1[ , c(10,11,12)], na.rm=TRUE)
df_fc_YB_gen1$SUM_F2_013YB <- rowSums(df_fc_YB_gen1[ , c(13,14,15)], na.rm=TRUE)
df_fc_YB_gen1$SUM_F2_022YB <- rowSums(df_fc_YB_gen1[ , c(16,17,18)], na.rm=TRUE)
df_fc_YB_gen1$SUM_F2_025YB <- rowSums(df_fc_YB_gen1[ , c(19,20,21)], na.rm=TRUE)

# create table for DEG Analysis
new_cols_YB = c(colnames(df_fc_YB_gen1))
df_fc_YB_sum = df_fc_YB_gen1[, str_subset(new_cols_YB, "^SUM\\w*")]
write.csv(df_fc_YB_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/YB_SUM_FC.csv", row.names = TRUE)

# F2 Generation, Bl6 mother, Cast-EiJ father
df_fc_XB_gen1 = df_feature_counts[, str_subset(column_names_vector, "^F2\\w*XB")]
meta_XB_gen1 = as.data.frame(data_meta_data_XB_gen1)

# add columns with sum of genome1, genome2 and all 
df_fc_XB_gen1$SUM_F2_016XB <- rowSums(df_fc_XB_gen1[ , c(1,2,3)], na.rm=TRUE)
df_fc_XB_gen1$SUM_F2_017XB <- rowSums(df_fc_XB_gen1[ , c(4,5,6)], na.rm=TRUE)
df_fc_XB_gen1$SUM_F2_018XB <- rowSums(df_fc_XB_gen1[ , c(7,8,9)], na.rm=TRUE)
df_fc_XB_gen1$SUM_F2_019XB <- rowSums(df_fc_XB_gen1[ , c(10,11,12)], na.rm=TRUE)
df_fc_XB_gen1$SUM_F2_020XB <- rowSums(df_fc_XB_gen1[ , c(13,14,15)], na.rm=TRUE)
df_fc_XB_gen1$SUM_F2_027XB <- rowSums(df_fc_XB_gen1[ , c(16,17,18)], na.rm=TRUE)
df_fc_XB_gen1$SUM_F2_028XB <- rowSums(df_fc_XB_gen1[ , c(19,20,21)], na.rm=TRUE)

# create table for DEG Analysis
new_cols_XB = c(colnames(df_fc_XB_gen1))
df_fc_XB_sum = df_fc_XB_gen1[, str_subset(new_cols_XB, "^SUM\\w*")]
write.csv(df_fc_XB_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/XB_SUM_FC.csv", row.names = TRUE)

# DEG analysis with DESeq 2
# BL6 sum
dds_BL6_sum = DESeqDataSetFromMatrix(countData = df_fc_BL6_sum, colData = meta_BL6_gen1, design = ~ condition)
dds_BL6_sum <- DESeq(dds_BL6_sum)
res_BL6_sum <- results(dds_BL6_sum)
res_LFC_BL6_sum <- lfcShrink(dds_BL6_sum, coef="condition_untreated_vs_treated", type="apeglm")
sum(res_BL6_sum$padj < 0.05, na.rm=TRUE)
sub_res_BL6_sum = subset(res_BL6_sum, res_BL6_sum$padj < 0.05)
write.csv(sub_res_BL6_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/DEG_BL6_sum.csv", row.names = TRUE)

# CC sum
dds_CC_sum = DESeqDataSetFromMatrix(countData = df_fc_CC_sum, colData = meta_CC_gen1, design = ~ condition)
dds_CC_sum<- DESeq(dds_CC_sum)
res_CC_sum <- results(dds_CC_sum)
res_LFC_CC_sum <- lfcShrink(dds_CC_sum, coef="condition_untreated_vs_treated", type="apeglm")
sum(res_CC_sum$padj < 0.05, na.rm=TRUE)
sub_res_CC_sum = subset(res_CC_sum, res_CC_sum$padj < 0.05)
write.csv(sub_res_CC_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/DEG_CC_sum.csv", row.names = TRUE)

# YB sum
dds_YB_sum = DESeqDataSetFromMatrix(countData = df_fc_YB_sum, colData = meta_YB_gen1, design = ~ condition)
dds_YB_sum <- DESeq(dds_YB_sum)
res_YB_sum <- results(dds_YB_sum)
res_LFC_YB_sum <- lfcShrink(dds_YB_sum, coef="condition_untreated_vs_treated", type="apeglm")
sum(res_YB_sum$padj < 0.05, na.rm=TRUE)
sub_res_YB_sum = subset(res_YB_sum, res_YB_sum$padj < 0.05)
write.csv(sub_res_YB_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/DEG_YB_sum.csv", row.names = TRUE)

#XB sum
dds_XB_sum = DESeqDataSetFromMatrix(countData = df_fc_XB_sum, colData = meta_XB_gen1, design = ~ condition)
dds_XB_sum <- DESeq(dds_XB_sum)
res_XB_sum <- results(dds_XB_sum)
res_LFC_XB_sum <- lfcShrink(dds_XB_sum, coef="condition_untreated_vs_treated", type="apeglm")
sum(res_XB_sum$padj < 0.05, na.rm=TRUE)
sub_res_XB_sum = subset(res_XB_sum, res_XB_sum$padj < 0.05)
write.csv(sub_res_XB_sum,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/DEG_XB_sum.csv", row.names = TRUE)


#VennDiagram
myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(rownames(sub_res_BL6_sum), rownames(sub_res_CC_sum), rownames(sub_res_YB_sum), rownames(sub_res_XB_sum)),
  category.names = c("BL6","CC", "YB", "XB"),
  filename = '/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/DEG_venn_diagramm_all.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans"
)

                
#MAplots
plotMA(res_LFC_BL6_sum, ylim=c(-2,2))
plotMA(res_LFC_CC_sum, ylim=c(-2,2))
plotMA(res_LFC_YB_sum, ylim=c(-2,2))
plotMA(res_LFC_XB_sum, ylim=c(-2,2))

#PCA Plot Samples
df_BL6_pca <- data.frame(t(df_fc_BL6_sum), rownames=TRUE)
df_BL6_pca <- df_BL6_pca[ , which(apply(df_BL6_pca, 2, var) != 0)]
meta_pca_BL6 = data.frame("samples"=c("11","12","13","14","15","16","17","18"), "tr_untr"=c("treated","treated","treated","treated","untreated","untreated","untreated","untreated"))
pca_res_BL6 <- stats::prcomp(df_BL6_pca, scale. = TRUE, center = TRUE)
PCA_BL6 <- autoplot(pca_res_BL6, data=meta_pca_BL6, shape='tr_untr', colour='samples', xlim = c(-1,1), ylim = c(-1,1)) 
summary(pca_res_BL6)

df_CC_pca <- data.frame(t(df_fc_CC_sum), rownames=TRUE)
df_CC_pca <- df_CC_pca[ , which(apply(df_CC_pca, 2, var) != 0)]
meta_pca_CC = data.frame("samples"=c("01","02","03","04","06","07"), "tr_untr"=c("treated","treated","treated","treated","untreated","untreated"))
pca_res_CC <- stats::prcomp(df_CC_pca, scale. = TRUE, center = TRUE)
PCA_CC <- autoplot(pca_res_CC, data=meta_pca_CC, shape='tr_untr', colour='samples', xlim = c(-1,1), ylim = c(-1,1))
meta_pca_CC_new = data.frame("samples"=c(colnames(df_fc_CC_sum)), "condition"=c("treated","treated","treated","treated","untreated","untreated"))
write.csv(meta_pca_CC_new,"/Users/julian/Documents/Bioinfo/Master/3_Semester/Praktikum/RNA_seq_final/Analysis/Sum_treated_untreated/META_CC_SUM_FC.csv", row.names = TRUE)


df_XB_pca <- data.frame(t(df_fc_XB_sum), rownames=TRUE)
df_XB_pca <- df_XB_pca[ , which(apply(df_XB_pca, 2, var) != 0)]
meta_pca_XB = data.frame("samples"=c("16","17","18","19","20","27","28"), "tr_untr"=c("treated","treated","untreated","untreated","untreated","treated","treated"))
pca_res_XB <- stats::prcomp(df_XB_pca, scale. = TRUE, center = TRUE)
PCA_XB <- autoplot(pca_res_XB, data=meta_pca_XB, shape='tr_untr', colour='samples', xlim = c(-1,1), ylim = c(-1,1))

df_YB_pca <- data.frame(t(df_fc_YB_sum), rownames=TRUE)
df_YB_pca <- df_YB_pca[ , which(apply(df_YB_pca, 2, var) != 0)]
meta_pca_YB = data.frame("samples"=c("08","09","10","12","13","22","25"), "tr_untr"=c("untreated","treated","untreated","untreated","treated","treated","treated"))
pca_res_YB <- stats::prcomp(df_YB_pca, scale. = TRUE, center = TRUE)
PCA_YB <- autoplot(pca_res_YB, data=meta_pca_YB, shape='tr_untr', colour='samples', xlim = c(-1,1), ylim = c(-1,1))

PCA_BL6
PCA_CC
PCA_XB
PCA_YB
