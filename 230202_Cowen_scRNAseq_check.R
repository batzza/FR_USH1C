
library(Seurat)
library(scCustomize)
library(fcoex)
library(org.Hs.eg.db)

## Load and normalize dataset
load("/Users/batzza/Desktop/EGAD00001006350.seurat.Rdata")
scEiaD <- NormalizeData(object = scEiaD)
scEiaD <- ScaleData(object = scEiaD)

## Set genes to check
# USH1C: ENSG00000006611
features <- c("ENSG00000006611")

## Set cell type as identity for object
scEiaD <- SetIdent(scEiaD, value = scEiaD@meta.data$CellType)

## Subset to remove NA cells
not_na_cells <- scEiaD@meta.data %>% filter(!is.na(CellType)) %>% rownames()
scEiaD <- subset(scEiaD, cells = not_na_cells)

## Find percent expressing features of interest by cell type
sort(unlist(Percent_Expressing(seurat_object = scEiaD, features = features)), decreasing = T)

## Check expression across all cell types
VlnPlot(scEiaD, features = features)
DotPlot(scEiaD, features = features) + RotatedAxis()

## ID Muller Glia Cells
mg_cells <- scEiaD@meta.data %>% filter(CellType == "Muller Glia") %>% rownames()

## Subset to MG cells
mg_seurat <- subset(x=scEiaD, cells = mg_cells)
mg_seurat <- SetIdent(mg_seurat, value = mg_seurat@meta.data$SubCellType)

## Check expression of USH1C
FeaturePlot(object = mg_seurat, features = features, split.by = 'SubCellType')
VlnPlot(mg_seurat, features = features)
Percent_Expressing(seurat_object = mg_seurat, features = features)
AverageExpression(mg_seurat, features=features)


## Set identity based on expression of USH1C
Idents(mg_seurat, WhichCells(object = mg_seurat, expression = ENSG00000006611 > 0, slot = 'data')) <- 'USH1C_positive'
Idents(mg_seurat, WhichCells(object = mg_seurat, expression = ENSG00000006611 == 0, slot = 'data')) <- 'USH1C_negative'

# Find differentially expressed features between between MG cells with and without USH1C
ush1c_mg.de.markers <- FindMarkers(mg_seurat, ident.1 = "USH1C_positive", ident.2 = "USH1C_negative", min.diff.pct = 0.001, min.pct = 0.001, logfc.threshold = 0.001)

markersPlot <- data.frame(ENSEMBL=ush1c_mg.de.markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames())

## Change from row IDs to gene names
annots <- select(org.Hs.eg.db, keys=markersPlot$ENSEMBL, 
                 columns="SYMBOL", keytype="ENSEMBL")

annots <- annots %>% mutate(NAMES = coalesce(SYMBOL,ENSEMBL))

## Plot top markers
DotPlot(mg_seurat, features = markersPlot) + RotatedAxis() + scale_x_discrete(labels=annots$NAMES)
