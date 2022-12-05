## Find go terms for both upregulated and downregulated
joined_GO_filtered_arranged <- bind_rows(up_GO, down_GO) %>%
  # Only retain significant rows
  filter(weight01 <= 0.05) %>%
  # Calculate the enrichment ratio
  mutate(GeneRatio = Significant/Annotated) %>%
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term)) %>%
  head(n=40)

# Now lets extract the order of the term column
order_term_joined <- joined_GO_filtered_arranged %>% select(Term) %>% pull()

down_GO_joined <- ggplot(joined_GO_filtered_arranged, aes(x= Term, y = GeneRatio, color = as.numeric(weight01))) +
  geom_point(aes(size= Significant)) +
  coord_flip() +
  scale_x_discrete(limits = order_term_joined) +
  scale_color_gradient(low = "orange", high = "blue") + 
  scale_y_continuous(limits = c(0,0.7), breaks = seq(0,0.7,0.1), expand = c(0,0)) + 
  theme_light() +
  labs(x = "GO term description", y = "Enrichment Ratio", color = "p-value", size= "# of sig. genes") + 
  ggtitle("Go Term Enrichment for Downregulated and Upregulated Genes") +
  facet_grid(.~ up_down) +
  theme(panel.border = element_rect(color = "black"),
      panel.grid = element_line(colour = "grey96"))
