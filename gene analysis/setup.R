## Download packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install packages
install.packages("dplyr")
install.packages("ggplot2")
BiocManager::install("DESeq2")

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


## Setup
# Data here: 
# https://drive.google.com/drive/folders/1S3p01KEA_BJ4dlH-gSu4K6pd3ZUPi8qd

# Accomondate your local file structure of the downloaded folder (assumes that you've moved the whole project_essentials folder under your project directory AND unzipped it). Also, rename the .tab files to .tsv before running the code below
setwd('./project_essentials-20221124T191810Z-001/project_essentials/star_tabs')

# Treatment rep 1
UV_1h_rep1 <- read_tsv("treat_SRX1023257ReadsPerGene.out.tsv",
                      col_names = c("gene_id", "total",
                                    "antisense", "sense"),
                      skip = 4)
# Treatment rep 2
UV_1h_rep2 <- read_tsv("treat_SRX1023258ReadsPerGene.out.tsv",
                      col_names = c("gene_id", "total",
                                    "antisense", "sense"),
                      skip = 4)
# Control rep 1
UV_0h_rep1 <- read_tsv("ctrl_SRX1023254ReadsPerGene.out.tsv",
                      col_names = c("gene_id", "total",
                                    "antisense", "sense"),
                      skip = 4)

# Control rep 2
UV_0h_rep2 <- read_tsv("ctrl_SRX1023255ReadsPerGene.out.tsv",
                      col_names = c("gene_id", "total",
                                    "antisense", "sense"),
                      skip = 4)

# Compile all of our loaded files into a dataframe
# We have stranded data (single) and use the sense strand
dat <- data.frame(row.names = UV_1h_rep1$gene_id,
                                UV_0h_rep1 = UV_0h_rep1$sense,
                                UV_0h_rep2 = UV_0h_rep2$sense,
                                UV_1h_rep1 = UV_1h_rep1$sense,
                                UV_1h_rep2 = UV_1h_rep2$sense)

# Transform dat into a matrix
dat_matrix <- as.matrix(dat) 

# Metadata file whose rownames match the name and the order of the columns of our matrix above
# Make our Metadata file that contains our column information for our matrix
metadata <- data.frame(row.names = colnames(dat_matrix), 
                       condition = c("UV_0h", "UV_0h", "UV_1h", "UV_1h")
                       )

# Double check that the row names from our metadata file match the column names of our matrix
colnames(dat_matrix) == rownames(metadata)

## DESeq2

# Create our DESeq2 object
dds_matrix <- DESeqDataSetFromMatrix(countData = dat_matrix, #matrix 
                              colData = metadata, #metadata file
                              design = ~condition)

# Set control condition using the relevel function for our 0h control
dds_matrix$condition <- relevel(dds_matrix$condition, ref = "UV_0h")

# Run DESeq2 on the dataset
dds <- DESeq(dds_matrix)

# Because the distribution of RNA Seq data is highly skewed, with a few high abundance genes and many low abundance genes, it is helpful to transform our data to visualize clustering. DESeq2 comes with the function rlog()

rld <- rlog(dds)
# Use PCA to identify which samples are more similar and if they group by one or more of the independent variables (in our case, we have only a single variable that can take “control” or “treated”).
plotPCA(rld, intgroup = "condition")

# Names of the results that DESeq2 calculated for next step
resultsNames(dds)

# Now we will extract the results for our comparison between the 1h timepoint and the 0h timepoint
res <- results(dds, name = "condition_UV_1h_vs_UV_0h") %>% as.data.frame() # 

# Filter out any row with an NA value
res_no_NA <- res %>% drop_na()

# Put the rownames as a column so they can be saved in your CSV file
res_no_NA_final <- res_no_NA %>% rownames_to_column("gene_id")