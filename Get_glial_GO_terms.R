## Script to pull glia related GO pathways and associated genes
library(dplyr)

## Get all the GO terms
library(GO.db)
GO_terms <- select(GO.db, keys= keys(GO.db), columns = c("GOID","TERM","ONTOLOGY","DEFINITION"))


## Filter down to terms that include "glia" or "Gliosis", or "muller"'
library(stringr)
search_terms <- c("gliosis", "glia","müller")
select_GO_terms <- GO_terms[grep(paste(search_terms,collapse = "|"),GO_terms$TERM),1]


## Pull down the genes related to the assigned terms
library(mgsa)
GO_anno <- readGAF("data/external/goa_human.gaf")

## Filter to GO terms found in previous step
hits <- GO_anno@sets[grep(paste(select_GO_terms,collapse="|"), names(GO_anno@sets))] ## Filter to selected lists
hits <- tibble(GO_term_ID = names(hits), gene_index = matrix(hits))          ## Get gene indices for each GO term
hits <- hits %>% tidyr::unnest_longer(gene_index)                     ## Reorganize so each row has a gene/Go term pair
hits$Uniprot_ID <- names(GO_anno@itemName2ItemIndex[hits$gene_index]) ## Add Uniprot gene IDs to tibble

## Add gene info
hits <- hits %>% inner_join(data.frame(Uniprot_ID = rownames(GO_anno@itemAnnotations), 
                                       gene_symbol = GO_anno@itemAnnotations$symbol, 
                                       gene_name = GO_anno@itemAnnotations$name))
## Add GO term info
hits <- hits %>% inner_join(data.frame(GO_term_ID = rownames(GO_anno@setAnnotations), 
                                       GO_term_name = GO_anno@setAnnotations$term, 
                                       GO_term_definition = GO_anno@setAnnotations$definition))



## I previously used biomaRt for this but the results are were not great
## When I manually checked the gene lists associated with particular GO terms
## Many genes were missing for most terms
## Instead, I have downloaded a static copy of the GO annotations from amiGO
## This is more up to date and seems to work better
## Leaving this old code for biomaRt commented out here in case we need in the future

# library("biomaRt")
# ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
# 
# select_GO_genes <- getBM(c("go_id","ensembl_gene_id","hgnc_symbol"), 
#                          filters = c("go","chromosome_name"), 
#                          values = list(select_GO_terms$GOID,c(1:22,"X","Y")), 
#                          ensembl) %>% 
#   filter(go_id %in% select_GO_terms$GOID) %>%
#   arrange(go_id, hgnc_symbol)



## Write results out
write.csv(x = hits, file = "data/processed/Glial_GO_terms.csv",quote = T, row.names = F)

