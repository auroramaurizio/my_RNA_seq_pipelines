library(Seurat)
library(ggplot2)

stromaCAFs <- readRDS("stromaCAFs")

stromaCAFs$customclassif.stim <- paste(Idents(stromaCAFs), stromaCAFs$cond, sep = " ")
stromaCAFs$customclassif <- Idents(stromaCAFs)
Idents(stromaCAFs) = stromaCAFs$customclassif.stim

AverageExpression_new <- AverageExpression(stromaCAFs, return.seurat = T, group.by = "customclassif.stim")

AverageExpression_new$customclassif.stim <- names(Idents(AverageExpression_new))


AverageExpression_new$cond <- c("Ctrl","leuk-A","leuk-B","Ctrl","leuk-A","leuk-B","Ctrl","leuk-A","leuk-B" )
AverageExpression_new$customclassif <- c("apCAF","apCAF","apCAF","iCAF","iCAF","iCAF","myCAF","myCAF","myCAF")

group.colors = c("apCAF Ctrl" = "#F8766D", "iCAF Ctrl" = "#F8766D", "myCAF Ctrl" = "#F8766D",
                 "apCAF leuk-A" = "#00BA38", "iCAF leuk-A" = "#00BA38", "myCAF leuk-A" = "#00BA38",
                 "apCAF leuk-B" = "#619CFF", "iCAF leuk-B" = "#619CFF", "myCAF leuk-B" = "#619CFF") 


pdf("AllMarkers_heatmap_cond_test.pdf", 2.8, 12)
DoHeatmap(AverageExpression_new, group.by = "customclassif.stim", size = 2, raster = FALSE, group.colors = c("apCAF Ctrl" = "#F8766D", "iCAF Ctrl" = "#F8766D", "myCAF Ctrl" = "#F8766D",
                                                                                                             "apCAF leuk-A" = "#00BA38", "iCAF leuk-A" = "#00BA38", "myCAF leuk-A" = "#00BA38",
                                                                                                             "apCAF leuk-B" = "#619CFF", "iCAF leuk-B" = "#619CFF", "myCAF leuk-B" = "#619CFF"), 
          features = top10$gene, draw.lines = F, angle = 90) + NoLegend()
dev.off()

pdf("AllMarkers_heatmap_customc_test.pdf", 2.8, 12)
DoHeatmap(AverageExpression_new, group.by = "customclassif", size = 2, raster = FALSE, features = top10$gene, draw.lines = F, angle = 90, group.colors = c('myCAF' = '#fbb040',
                                                                                                                                                           'iCAF' = '#9e1f63',
                                                                                                                                                           'apCAF'='#662d91') ) + NoLegend()
dev.off()
