#devtools::install_github("xmc811/Scillus", ref = "development")
library(Scillus)
library(Seurat)
library(ggplot2)

stromaCAFs <- readRDS("stromaCAFs")

cluster.markers <- FindAllMarkers(stromaCAFs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Seurat find all markers heatmap with double annotatiotion coloured by condition and celltype

pdf("AllMarkers_heatmap_CAF_RNA_avgexp_split_by_cond.pdf", 10, 10)
Scillus::plot_heatmap(dataset = stromaCAFs, 
                      markers = top10$gene,
                      sort_var = c("customclassif","cond"),
                      anno_var = c("customclassif","cond"),
                      hm_colors = c("purple","black","yellow"),
                      hm_limit = c(-1,0,1),
                      anno_colors = list(c('apCAF'='#662d91','iCAF' = '#9e1f63', 'myCAF' = '#fbb040'),
                                         c("Ctrl" = "#F8766D","leuk-A" = "#00BA38", "leuk-B" = "#619CFF"))
)