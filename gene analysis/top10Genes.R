## Visualizing the top 10 highest and lowest up/down-regulated genes

# Top 10 genes positively differentially expressed
top10_genes <- res_no_NA_final %>%
  arrange(desc(log2FoldChange)) %>%
  # NOTE that we use the desc() function to organize the column in descending order
  head(n = 10)

# Top 10 genes negatively differentially expressed
bot10_genes <- res_no_NA_final %>%
  arrange(log2FoldChange) %>% 
  # NOTE since we don't use desc(), the column is organized in ascending order
  head(n = 10)

# Bar plots
top10_plot_horizontal_SE <- ggplot(top10_genes, aes(x = gene_id, y = log2FoldChange)) +
  # Lets add errorbars! We also change the width= of our error bars to make them nicer!
  geom_errorbar(aes(ymin= log2FoldChange - lfcSE, ymax =log2FoldChange + lfcSE), width = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = order) +
  coord_flip()

bot10_plot_horizontal_SE <- ggplot(bot10_genes, aes(x = gene_id, y = log2FoldChange)) +
  # Lets add errorbars! We also change the width= of our error bars to make them nicer!
  # DESeq2 calculates the standard error for the log2FoldChange called lfcSE
  geom_errorbar(aes(ymin= log2FoldChange - lfcSE, ymax =log2FoldChange + lfcSE), width = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = order) +
  coord_flip()