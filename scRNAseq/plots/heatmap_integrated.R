
#suppressPackageStartupMessages(library(monocle3)) # OK
suppressPackageStartupMessages(library(monocle)) # OK
suppressPackageStartupMessages(library(SeuratWrappers)) # NO
suppressPackageStartupMessages(library(Seurat)) # OK
#suppressPackageStartupMessages(library(SeuratData)) # NO
suppressPackageStartupMessages(library(ggplot2)) # OK
suppressPackageStartupMessages(library(patchwork)) # OK
suppressPackageStartupMessages(library(magrittr)) # OK
suppressPackageStartupMessages(library(future)) # OK
suppressPackageStartupMessages(library(cowplot)) # OK
suppressPackageStartupMessages(library(dplyr)) # OK
library("ggalluvial")
###################################################################

#saveRDS(integrated, file = "integrated_subset_clust_renamed_30pc.Rds")



C_high = readLines("C_high.txt")
C_medium = readLines("C_medium.txt")


C_high_medium = c(C_high,C_medium)




annotation_row = C_high_medium
length(C_medium)
pathways_C_high = rep("C_high", 93)
pathways_C_medium = rep("C_medium", 154)
pathways = c(pathways_C_high, pathways_C_medium)

DefaultAssay(integrated) = "RNA"

# sample 1 - intact
sample1 <- integrated_3TIP_new
avgexp1 = AverageExpression(sample1, return.seurat = T)
counts1 <- GetAssayData(avgexp1, assay="RNA", slot="data")
genes <- annotation_row
counts1 <- as.matrix(counts1[rownames(counts1) %in% genes, ])
cluster1 <- colnames(counts1)
colnames(counts1) = colnames(counts1)
sample1 <- colnames(counts1)


sample <- sample1
cluster <- cluster1
stim1 <- rep("integrated",6)
#stim2 <- rep("2", 9)
stim <- stim1

df <- data.frame(sample, cluster)
head(df)
df = df[order(df$sample),]

head(df)
counts = counts1
head(counts)
annotation_column = as.data.frame(df[,c('cluster')])
colnames(annotation_column)= "cluster"
rownames(annotation_column) <- df$cluster
head(annotation_column)



ann_colors = list(
  pathways = c("C_high" = '#FF007F',
               "C_medium" = '#FFCCFF'),
  cluster = c('A' = '#FF6666',
              'B' = '#990000',
              'C' = '#336666',
              'D' = '#6600CC',
              'E' = '#FF99CC',
              'F' = '#FF00FF',
              'G' = '#399933',
              'H' = '#0066FF',
              'I' = '#99CC33'))



annotation_row = annotation_row
counts_ordered = counts[,row.names(annotation_column)]
counts_ordered = counts_ordered[annotation_row,]

head(counts_ordered)

sub_samp_ordered <- counts_ordered


tail(counts_ordered)
tail(annotation_row)
#sub_samp_ordered <- counts_ordered[annotation_row,,drop=FALSE]
#sub_samp_ordered = sub_samp_ordered[sub_samp_ordered$name != "Gm42418", ]

#my.breaks <- c(seq(-3, 3, by=0.01))

#remove = c("Fabp4", "Gm42418")
C <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = C(255)
minH = -1; maxH=1
myb = seq(minH, maxH, by = 0.01)
#C <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- C(length(myb))

annotation_row
pathways

length(pathways)
df_row <- data.frame(annotation_row, CellTypes)
row_annotation = df_row[2]
row.names(df_row) = df_row$annotation_row


draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)


library("pheatmap")

length(pathways)
df_row <- data.frame(annotation_row, pathways)
row_annotation = df_row[2]
row.names(df_row) = df_row$annotation_row


pheatmap(sub_samp_ordered,
         show_rownames = F,
         show_colnames = F,
         cellheight = 1,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = annotation_column,
         annotation_colors = ann_colors,
         width = 5,
         height = 5,
         scale = 'row',
         color = myc,
         annotation_row = df_row[2],
         filename = 'heatmap_C_small_RNA_new.pdf')



pheatmap::pheatmap(sub_samp_ordered1,
                   show_rownames = T,
                   show_colnames = F,
                   cluster_rows = F,
                   cellheight = 8,
                   cluster_cols = F,
                   annotation_col = annotation_column,
                   annotation_colors = ann_colors,
                   scale = 'row',
                   color = myc,
                   annotation_row = df_row[2],
                   filename = 'heatmap_C_big_RNA_new.pdf')




