## Script to pull glia related GO/KEGG/Reactome pathways
library(dplyr)

## Get all the GO terms
library(GO.db)
GO_terms <- select(GO.db, keys= keys(GO.db), columns = c("GOID","TERM","ONTOLOGY","DEFINITION"))

# ## Get all the Reactome terms
# library(reactome.db)
# keys <- keys(reactome.db)
# Reactome_terms <- select(reactome.db, keys= keys, columns=c("PATHID","PATHNAME"))
# Reactome_terms <- Reactome_terms %>% filter(str_detect(PATHNAME, "Homo sapiens"))

## Filter down to terms that include "glia" or "Gliosis", or "muller"'
library(stringr)
search_terms <- c("gliosis", "glia","müller")
select_GO_terms <- GO_terms[grep(paste(search_terms,collapse = "|"),GO_terms$TERM),]
select_Reactome_terms <- Reactome_terms[grep(paste(search_terms,collapse = "|"),Reactome_terms$PATHNAME),]

## Pull down the genes related to the assigned terms
library("biomaRt")
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

select_GO_genes <- getBM(c("go_id","ensembl_gene_id","hgnc_symbol"), 
                         filters = c("go","chromosome_name"), 
                         values = list(select_GO_terms$GOID,c(1:22,"X","Y")), 
                         ensembl) %>% 
  filter(go_id %in% select_GO_terms$GOID) %>%
  arrange(go_id, hgnc_symbol)

## Clean up and export
export <- select_GO_terms %>% 
  inner_join(select_GO_genes, by=c("GOID"="go_id")) %>%
  group_by(GOID, TERM, ONTOLOGY, DEFINITION) %>%
  summarise(gene_ids = paste(ensembl_gene_id, collapse = "; "), gene_names = paste(hgnc_symbol, collapse = "; "))

write.csv(x = export, file = "Glial_GO_terms.csv",quote = T, row.names = F)
