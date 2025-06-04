library(tidyverse)
library(pheatmap)

fCountsRPKMTOP <- read.xlsx("fCountsRPKMTOP.xlsx")
annotation_column <- read.xlsx("annotation_column.xlsx")

ann_colors = list(
  Sample = c("1" = "forestgreen",
             "3" = "brown",
             "4" = "#999900"),
  Condition = c("LPS"="#FFFF66", "PH"="orange", "PM" = "pink", "SE" = "#66B2FF", "Vehicle" = "darkgreen"))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)

HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         annotation_colors = ann_colors,
                         cluster_rows = T,
                         cluster_cols = T,
                         #main = 'Heatmap: 500 most variable genes - RPKM',
                         show_rownames = F,
                         cutree_rows = 4,
                         #cutree_cols = 5,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                         display_numbers = F,
                         col=colors,
                         filename = 'Heatmap_top500rpkm_sample_all_samp_nosamp2_t.pdf',
                         width = 6, height = 6)
#qui

gene_clusters <- cutree(HP$tree_row, k = 4)  # You choose the number of clusters
gene_clusters <- data.frame(Gene = names(gene_clusters),
                            Cluster = factor(gene_clusters))

expr_long <- fCountsRPKMTOP %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(gene_clusters, by = "Gene") %>%
  left_join(annotation_column %>% rownames_to_column("SampleID"), by = c("Sample" = "SampleID"))

expr_long <- expr_long %>%
  group_by(Gene) %>%
  mutate(Expression_z = scale(Expression)) %>%
  ungroup()

expr_long$Condition <- factor(expr_long$Condition, 
                              levels = c("Vehicle", "SE", "PM", "PH", "LPS"))


pdf("spaghetti.pdf")
ggplot(expr_long, aes(x = Condition, y = Expression_z, group = Gene)) +
  geom_line(alpha = 0.2, color = "grey") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", size = 1) +
  facet_wrap(~ Cluster, scales = "free_y", ncol = 1) +
  scale_x_discrete(limits = c("Vehicle", "SE", "PM", "PH", "LPS")) +
  theme_minimal() +
  labs(y = "Z-scored Expression", x = NULL, title = "Spaghetti Plots by Gene Cluster") +
  theme(strip.text = element_text(size = 12))
dev.off()


