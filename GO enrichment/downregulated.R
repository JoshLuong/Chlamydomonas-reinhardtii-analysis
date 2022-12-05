# Load data for downregulated genes
down_GO <- down_summary %>% mutate(up_down = "downregulated")

# Filter out any non-significant data and calculate the gene ratio
down_GO_filtered <- down_GO %>%

# Only retain significant rows
filter(weight01 <= 0.05) %>%

# Calculate the enrichment ratio
mutate(GeneRatio = Significant/Annotated) %>%

# Plot the downregulated go terms
down_GO_plot <- ggplot(down_GO_filtered, aes(x= Term, y = GeneRatio)) +
  # Lets add very thin bars (this goes really well with geom point!)
  geom_col(width = 0.1) +
  geom_point(size = 3) +
  
  # Flip so the x-axis labels are readable
  coord_flip()

# First, lets arrange the data based on the enrichment ratio
down_GO_filtered_arranged <- down_GO_filtered %>%
arrange(GeneRatio) %>%
# Make our 'Term' column a factor so that we can organize it!
mutate(Term = factor(Term))

# Now lets extract the order of the term column
order_term <- down_GO_filtered_arranged %>% select(Term) %>% pull()

down_GO_plot_final <- ggplot(down_GO_filtered_arranged, aes(x= Term, y = GeneRatio, color = as.numeric(weight01))) +
  # Lets add very thin bars (this goes really well with geom point!)
  geom_col(width = 0.05) +
  
  # Since we just want the size of the point, we will add it just to the geom_point()
  geom_point(aes(size= Significant),) +
  # scale the point size up
  scale_size(range = c(2,12)) +
  coord_flip() +
  scale_x_discrete(limits = order_term) +
  # Scale to 1.1 to see last point more easily
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.25), expand = c(0,0)) +
  # Lets change the color to have a better gradient
  scale_color_gradient(low = "orange", high = "blue") + 
  theme_light() +
  # Title
  ggtitle("Go Term Enrichment for Downregulated Genes") +
  
  # Change the axis tittles
  labs(x = "GO term description", y = "Enrichment Ratio", color = "p-value", size= "# of sig. genes") + 
  
  # Make the plot prettier
  theme(panel.border = element_rect(color = "black"),
      panel.grid = element_line(colour = "grey96")) +
  # Change text size
  theme(text = element_text(size = 20))     
