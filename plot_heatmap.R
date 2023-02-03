## Load library for heatmap plotting
library(ComplexHeatmap)
library(tidyverse)
library(circlize)

## If you get an error that there is no library called ComplexHeatmap
## Then you need to install it first
## You can do that by uncommenting the next 3 lines and running this code 
## you should only need to do this once and only if you don't already have the package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")



## If you get an error that there is no library called dplyr
## Uncomment and run the following line once:
#install.packages("tidyverse")



## If you get an error that there is no library called circilize
## Uncomment and run the following line once:
#install.packages("circlize")




## Set the working directory to the script file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load dataset in same folder as script
cpm_data <- read.csv("Emily_cpm.csv")

## Extract sample metadata from CPM file
meta <- data.frame(sample_name = colnames(cpm_data[2:ncol(cpm_data)]))
meta <- cbind(meta,data.frame(do.call('rbind', strsplit(as.character(meta$sample_name),'.',fixed=TRUE))))
colnames(meta) <- c('sample_name','genotype','age','replicate')

## Change meta age to numeric for sorting puprposes
meta$age <- as.numeric(gsub(x=meta$age,pattern="D",replacement = ""))


## Select genes of interest to map
genelist <- c("TSPAN6", "DPM1", "FUCA2", "OTX2", "CRX")


## Filter data
mat <- cpm_data %>%
  filter(X %in% genelist) %>%              ## Filter the CPM data
  column_to_rownames(var = "X") %>%        ## Move gene name to rowname to allow conversion to numeric matrix
  as.matrix()                              ## Conver to matrix

## Convert to log2(CPM+1)
mat <- log2(mat+1)

## Set colors for heatmap
col_fun = colorRamp2(c(0,max(mat)), c("white", "red"))  ## Currently using "max(mat)" so the high end of the color
                                                        ## scale matches the maximum value in the selected genes
                                                        ## If you want to set arbitrarily, change to soemthing like:
                                                        ## colorRamp2(c(0,5), c("white", "red"))
  
## Plot
p <-Heatmap(mat,
    col = col_fun,                         ## Set heatmap colors to function defined above
    cluster_columns=F,                     ## Turn sample clustering on or off -- if clustering is off, samples are in order of input data
    cluster_rows=F,                        ## Turn gene clustering on or off -- if off, genes are in order of input data
    column_split=meta$age,                 ## Split heatmap display by metadata -- can be "meta$genotype" or "meta$age"
    column_gap = unit(5, "mm") ,           ## Set size of gaps between split groups
    name="meta",                           ## Set source for split group names
    bottom_annotation = 
      HeatmapAnnotation(
        Genotype = meta$genotype,               ## Add annotation row for patient/control
        col = list(Genotype = c("Control" = "red", "Patient" = "blue")) ## Set annotation colors
        ),
    heatmap_legend_param = list(title="log2(CPM+1)")  ## Set title of legend
    )

## Show plot
p

## Save plot
pdf(file = paste0(format(Sys.time(), "%Y-%m-%d_%H-%M"),"_Gene_heatmap.pdf"),  ## Automatically name file with date/time
    height = 8, width = 10.5, onefile=FALSE)                            ## Set size of plot
p
dev.off()


