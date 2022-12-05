# Extracting specific gene IDs we want to visualize in the heatmap:

# Top 10 genes that we upregulated
up_reg_vec_c <- c("CHLRE_04g217948v5","CHLRE_06g278267v5", "CHLRE_04g211850v5", "CHLRE_08g367500v5", "CHLRE_01g021800v5", "CHLRE_07g344500v5", "CHLRE_07g320450v5", "CHLRE_07g320400v5", "CHLRE_01g016750v5", "CHLRE_01g016600v5")

# Top 10 genes that we downregulated
down_reg_vec_c <- c("CHLRE_03g146667v5","CHLRE_02g118350v5", "CHLRE_06g296700v5", "CHLRE_03g143967v5", "CHLRE_17g728483v5", "CHLRE_16g680700v5", "CHLRE_06g297049v5", "CHLRE_06g263200v5", "CHLRE_09g391650v5", "CHLRE_01g039450v5")

# Top 10 genes that were upregulated and downregulated
up_down_reg_vec_c <- c(
  "CHLRE_04g217948v5","CHLRE_06g278267v5", "CHLRE_04g211850v5", "CHLRE_08g367500v5", "CHLRE_01g021800v5", "CHLRE_07g344500v5", "CHLRE_07g320450v5", "CHLRE_07g320400v5", "CHLRE_01g016750v5", "CHLRE_01g016600v5",
"CHLRE_03g146667v5","CHLRE_02g118350v5", "CHLRE_06g296700v5", "CHLRE_03g143967v5", "CHLRE_17g728483v5", "CHLRE_16g680700v5", "CHLRE_06g297049v5", "CHLRE_06g263200v5", "CHLRE_09g391650v5", "CHLRE_01g039450v5"
)

# Matrix of top genes that we upregulated and downregulated
up_down_reg_matrix <- dat_matrix[up_down_reg_vec_c ,]
# Create heatmap data from matrix
data <- as.matrix(up_down_reg_matrix)

# Generate heatmap. Remove dendrogram
gene_heatmap <- heatmap(data, scale= "column", margins=c(10,7), Rowv = NA, Colv = NA)+
  geom_tile()