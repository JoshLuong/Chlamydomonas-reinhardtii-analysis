# Analyzing Chlamydomonas reinhardtii’s stress response: Will Chlamydomonas reinhardtii experience effects similar gene expression to other higher order organisms 

**Authors: Bea Liston, Pardis Shirkan, James Lo, Joshua Luong, Uzak Terzioglu, Nicole Robinson**

This study aims to investigate differences in gene expression in Chlamydomonas reinhardtii in exposure to no UV-B light and with UV-B light. The focus will initially start with photosynthetic genes (e.g. LHC-related) and then further explore other stress-related genes.

In order to investigate this, conditions of no UV-B light and after a 1-hour exposure to UV-B light will be retrieved from the study by Tilbrook et al. (2016) with 2 replicates of each condition.

Bioinformatic tools were used such as STAR for single read RNA alignment, FASTQC to check for quality control, gene ontology (GO) terms to investigate the pathway and related genes, and R (DeSeq2 to analyze RNA-Seq data, plotting the log2FoldChange of top 10 genes that were up or down regulated, and other plotting tools) to visualize the data collected.

This repository houses the scripts we used for:

1. STAR alignment, FASTQC analysis (`/STAR/alignment.sh`)
2. GO Enrichment (`/GO enrichment`)
3. Top 10 upregulated and downregulated gene analysis (`/gene analysis/top10Genes.R`)
4. Heatmap analysis (`/gene analysis/heatmap.R`)

## STAR
The `alignment.sh` script was used to retrieve necessary data, check `md5sums`, generate FASTQCs, and finally, align our data using STAR.

**This should be run before any subsequent analysis.**

## Gene Analysis
The `gene analysis\setup.R` **must be run** before running any dependent analysis under `gene analysis/`.

### Heatmap
Run this code to generate a heatmap.
This heatmap uses the top 10 upregulated and top 10 downregulated genes to see the gene expression of our top 20 genes before and after exposure.

This heatmap shows that our replicates show very similar gene expression. This allows us to know that the replicates we used are reliable.

### Top 10 Upregulated and Downregulated Genes
Run the code in this file to generate both top 10 upregulated genes and top 10 downregulated genes.

## GO Enrichment
The `GO enrichment\setup.R` **must be run** before running any dependent analysis under `GO enrichment/`.

### Down and upregulated
Run `downregulated.R` to find GO term enrichments for downregulated genes. Run `upregulated.R` to find GO term enrichments for upregulated genes.

### Both Down and upregulated GO enrichments
**Run `downregulated.R` and `upregulated.R` before proceeding.** This code (`summary.R`) is intended to generate a GO enrichment that summarizes both the upregulated genes and the downregulated genes in a single plot.
