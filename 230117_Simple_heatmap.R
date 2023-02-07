## Libraries Used In Script
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(tibble)

## Before you start, set the working directory to the KK_NRL-L75Pfs/rnaseq directory
## You can do this by selecting from the menu above Session > Set Working Directory > Choose Directory... 
## The neavigating to KK_NRL-L75Pfs/rnaseq
## Alternatively, you can edit the file paths below for loading/saving data

## Load the log2(CPM+1) counts for each gene (filtered at >1CPM in at least one sample)
logCPM <- read.csv("processed/Gene_filtlog2CPM_MSTR_1CPM.tsv", sep='\t')

# Heatmap colors
col <- c("#023858","#3690c0","#e0f3f8","#fff7bc","#fec44f","#ec7014","#fc4e2a","#bd0026","#800026")
col <- colorRampPalette(col)(n=1000)

## Prep data for heatmap
## Pick genes to show EITHER by ensembl Gene ID (ensembl v106)
gene_list <- c("ENSG00000000419","ENSG00000000457","ENSG00000001036")
gene_mat <- logCPM %>% filter(X %in% gene_list) %>% 
              select(X, Control_1:NRL_L75Pfs_3) %>% 
              column_to_rownames(var = "X") %>% 
              as.matrix()
## OR by gene_name
gene_list <- c("TSPAN6","CFH","DPM1")
gene_mat <- logCPM %>% filter(external_gene_name %in% gene_list) %>% 
  select(external_gene_name, Control_1:NRL_L75Pfs_3) %>% 
  column_to_rownames(var = "external_gene_name") %>% 
  as.matrix()


###### Heatmap by log2(CPM+1)
p_cpm <- Heatmap(gene_mat, name = "log2(CPM+1)", col = col, cluster_columns = F, cluster_rows = T)
p_cpm

## If you want to save the log2CPM+1 heatmap, run all three of these lines in a row
pdf(file = paste0("plots/",format(Sys.time(), "%Y_%m_%d-%H%M_"),"CPM_heatmap.pdf"), 
    width = 5, height = 10 , onefile = T)
p_cpm
dev.off()

## Matrix in Z scale format
gene_mat_Z <- t(scale(t(gene_mat)))


###### Heatmap by Z score
p_Z <- Heatmap(gene_mat_Z, name = "Z score", col = col, cluster_columns = F, cluster_rows = T)
p_Z

## If you want to save the Z heatmap, run all three of these lines in a row
pdf(file = paste0("plots/",format(Sys.time(), "%Y_%m_%d-%H%M_"),"Z_score_heatmap.pdf"), 
    width = 5, height = 10 , onefile = T)
p_Z
dev.off()

