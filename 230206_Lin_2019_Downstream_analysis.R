## 02-06-2023
## Downstream analysis of gene expression changes in USH1C Organoids -- external data -- Lin 2019 PMID:31652328
## And figure generation

## Load the upstream processed data
load("data/processed/Lin_2019_Processed_RNAseq_w_CPM.Rdata")

## PCA
library(PCAtools)
rownames(meta) <- c(paste0(meta$Group,"_",meta$Replicate))
p <- pca(gene.dge.filt$lcpm, metadata = meta, removeVar = 0.1)


## Variance per PC (Scree Plot)
screeplot(p,
          hline=75,
          vline = 5)

## PCA Plot
biplot(p, 
       colby = "Group", 
       colkey = c("#61ADF1", "#f76a06"),
                           encircle = T)+
  theme_classic()

## Eigencor plot - Pearson's correlation of factor to PCs with FDR corrected p-value
eigencorplot(p,
             col = colorRampPalette(c("white", "red3"))(20),
             components = getComponents(p, 1:5),
             metavars = c('Replicate','Group'),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'fdr')




library(gprofiler2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(dplyr)

## 


## Set focal comparison 
df <- DE$Control_v_NMDA  ## Control vs NDMA damaged

# ## Get log2 fold change 
# gene_list <- df$logFC
# 
# ## name the vector
# names(gene_list) <- df$ensembl_ID
# 
# ## Sort the list in decreasing order (required for clusterProfiler)
# gene_list = sort(gene_list, decreasing = TRUE)
# 
# 
# ## Run GSEA on GO terms
# gse <- gseGO(geneList=gene_list, 
#              ont ="ALL", 
#              keyType = "ENSEMBL", 
#              minGSSize = 3, 
#              maxGSSize = 3000, 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = org.Mm.eg.db, 
#              pAdjustMethod = "BH")
# 
# 
# 
# ## GO term Dotplot
# dotplot(gse, showCategory=40, split=".sign") + facet_grid(.~.sign)
# 
# 
# ## GO term Ridgeplot
# ridgeplot(gse, showCategory = 40) + labs(x = "enrichment distribution") + xlim(-3,3)
# 
# ## GO term GSEA plot
# ## Look at the list of GSEA gene sets and pick one to plot by ID
# ## You can check all available sets by uncommenting the next line
# # gse@result$Description 
# 
# genesets <- grep(pattern="glia", ignore.case = T, x = gse@result$Description)
# geneset = genesets[1]
# gseaplot(gse, by = "all", title = gse@result$Description[geneset], geneSetID = geneset)




## Analysis of just DE genes
DE_df <- df %>% filter(FDR < 0.05) %>% filter(logFC>1) ## Up 
DE_df <- df %>% filter(FDR < 0.05) %>% filter(logFC < -1) ## Down

## gProfiler
gprof_organism <- "mmusculus"
gostres <- gost(query = DE_df$ensembl_ID,
                organism = gprof_organism,
                user_threshold = 0.05,
                correction_method = "g_SCS")

## Interactive gProfiler plot
gostplot(gostres, capped = TRUE, interactive = TRUE)


## If there are a lot of significant terms in gProfiler output, we can simplify with revigo
library(rrvgo)

## Calculate similarity between terms
simMatrix <- calculateSimMatrix(gostres$result$term_id,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

## Reduce terms and set names for plotting
scores <- setNames(-log10(gostres$result$p_value), gostres$result$term_id)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mm.eg.db")

## Tree plot (grouped boxes)
treemapPlot(reducedTerms)


