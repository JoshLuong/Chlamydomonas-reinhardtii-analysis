## Download and install packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("dplyr")
install.packages("ggplot2")
BiocManager::install("DESeq2")
BiocManager::install("GO.db")
BiocManager::install("topGO")
BiocManager::install("tidyverse")
install.packages("pheatmap")

# Load libraries
library(tidyverse)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(topGO)

# Load the R Data
# https://drive.google.com/drive/folders/16UkXBTaFY6Qjo28LgcBAK7yvG0VBvLtH
load('./project_essentials-20221124T191810Z-001/project_essentials/R data/micb405_R.RData')

# Drop observations with NA 
res_no_NA <- res %>% drop_na()

# Drop statistically insignificant results (p > 0.05)
res_filtered <- res_no_NA %>% filter(padj <= 0.05)

# Read in the "gene universe" file
# Download chamy_GOIDS from lecture; move it to your working directory if it's not there
geneID2GO <- readMappings('chlamy_GOIDs.tsv')
geneUniverse <- names(geneID2GO)

# Upregulated genes
up_genes <- res_filtered %>% 
  # Filter for the log fold changes for upregulated genes
  filter(log2FoldChange >= 1) %>%
  rownames_to_column("gene_id") # Put the rownames as a column so they can be saved in your CSV file

# Downregulated genes
down_genes <- res_filtered %>% 
  # Filter for the log fold changes for downregulated genes
  filter(log2FoldChange <= -1) %>%
  rownames_to_column("gene_id") # Put the rownames as a column so they can be saved in your CSV file

# Process data
upregulated_genes <- as.character(up_genes$gene_id)
downregulated_genes <- as.character(down_genes$gene_id)

# Factor the names
up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_genes))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_genes))
names(up_gene_list) <- geneUniverse
names(down_gene_list)<- geneUniverse

# Build GOdata object for upregulated
up_GO_data <- new("topGOdata",
                  description = "UV_0h_1h_up",
                  # Biological process
                  ontology = "BP",
                  allGenes = up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)

# Build the GOdata object in topGO for downregulated
down_GO_data <- new("topGOdata",
                    description = "UV_0h_1h_down",
                    # Biological process
                    ontology = "BP",
                    allGenes = down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

# Run Fisher's exact test with weight01 algorithm
up_result <- runTest(up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")
down_result <- runTest(down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")

# Summarize the results
up_summary <- GenTable(up_GO_data,
                       weight01 = up_result,
                       orderBy = "up_result",
                       ranksOf = "up_result",
                       topNodes = 20)
down_summary <- GenTable(down_GO_data,
                       weight01 = down_result,
                       orderBy = "down_result",
                       ranksOf = "down_result",
                       topNodes = 20)
