
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




DefaultAssay(integrated) = "integrated"
res= 1
sub <- FindClusters(integrated, resolution = res)
DimPlot(sub, split.by = "stim")

subs = subset(sub, idents = c("0","1", "2", "3", "4", "5", "6","7","9","10"))

DimPlot(subs, split.by = "stim")

new.cluster.ids.lit <- c('A', 
                         'B',
                         'C',
                         'D',
                         'E',
                         'F',
                         'G',
                         'H',
                         'I',
                         'L'
)


names(new.cluster.ids.lit) <- levels(subs)
integrated_new <- RenameIdents(subs, new.cluster.ids.lit)

DimPlot(integrated_new)

DefaultAssay(integrated) = "RNA"

C_high = readLines("C_high.txt")
C_medium = readLines("C_medium.txt")



annotation_row = c(C_high,
                   C_medium)

length(C_high)
length(C_medium)
pathways_C_high = rep("C_high", 93)
pathways_C_medium = rep("C_medium", 154)

pathways  = c(pathways_C_high,
              pathways_C_medium)


# sample 1 - intact
sample1 <- subset(x = integrated_new, subset = stim == "1")
avgexp1 = AverageExpression(sample1, return.seurat = T)
counts1 <- GetAssayData(avgexp1, assay="RNA", slot="data")
genes <- annotation_row 
counts1 <- as.matrix(counts1[rownames(counts1) %in% genes, ])
cluster1 <- colnames(counts1)
colnames(counts1) = paste(colnames(counts1), "1", sep = "_")
sample1 <- colnames(counts1)

# sample 2 - injured
sample2 <- subset(x = integrated_new, subset = stim == "2")
avgexp2 = AverageExpression(sample2, return.seurat = T)
counts2 <- GetAssayData(avgexp2, assay="RNA", slot="data")
genes <- annotation_row 
counts2 <- as.matrix(counts2[rownames(counts2) %in% genes, ])
cluster2 <- colnames(counts2)
colnames(counts2) = paste(colnames(counts2), "2", sep = "_")
sample2 <- colnames(counts2)

sample <- append(sample1, sample2)
cluster <- append(cluster1, cluster2)
stim1 <- rep("1", 9)
stim2 <- rep("2", 9)
stim <- append(stim1, stim2)

df <- data.frame(sample, cluster, stim)
head(df)
#df$sample = paste(df$cluster, "_", df$stim, sep = "")
df = df[order(df$sample),]

head(df)
counts = cbind(counts1, counts2)
head(counts)
annotation_column = df[,c('cluster','stim')]

rownames(annotation_column) <- df$sample
head(annotation_column)

#ann_colors = list(
#    stim = c("limegreen","darkviolet"),
#    cluster = c('))

ann_colors = list(
  pathways = c("C_high" = '#FF007F',
               "C_medium" = '#FFCCFF'),
  stim = c("1" = '#4C9900',
           "2" = '#FFB266'),
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
head(counts_ordered)

sub_samp_ordered <- counts_ordered[annotation_row,,drop=FALSE]


C <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = C(255)
minH = -1; maxH=1
myb = seq(minH, maxH, by = 0.01)
#C <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- C(length(myb))


length(pathways)
df_row <- data.frame(annotation_row, pathways)
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



pheatmap(sub_samp_ordered, 
         show_rownames = F, 
         show_colnames = F,
         cellheight = 1,
         cluster_rows = F,  
         cluster_cols = F, 
         #width = 5, 
         #height = 5, 
         annotation_col = annotation_column, 
         annotation_colors = ann_colors, 
         scale = 'row', 
         color = myc, 
         annotation_row = df_row[2], 
         fontsize_row = 5,
         filename = 'heatmap_C_small_RNA_split.pdf')



pheatmap(sub_samp_ordered, 
         show_rownames = T, 
         show_colnames = F,
         cellheight = 5,
         cluster_rows = F,  
         cluster_cols = F, 
         #width = 5, 
         #height = 5, 
         annotation_col = annotation_column, 
         annotation_colors = ann_colors, 
         scale = 'row', 
         color = myc, 
         annotation_row = df_row[2], 
         fontsize_row = 5,
         filename = 'heatmap_C_big_RNA_split.pdf')






toMatch = c("gene1","gene2","gene3")


matches <- grep(paste(toMatch,collapse="|"), 
                df_row$annotation_row)

matched <- grep(paste(toMatch,collapse="|"), 
                df_row$annotation_row, value = T)

ha = rowAnnotation(foo = anno_mark(at = matches, labels = matched))

pdf("heatmap.pdf", 20, 10)
ComplexHeatmap::pheatmap(sub_samp_ordered, 
                         show_rownames = F, 
                         show_colnames = F,
                         cellheight = 3,
                         cellwidth = 50,
                         cluster_rows = F,  
                         cluster_cols = F, 
                         #width = 5, 
                         #height = 5, 
                         annotation_col = annotation_column, 
                         annotation_colors = ann_colors, 
                         scale = 'row', 
                         color = myc, 
                         annotation_row = df_row[2], 
                         fontsize_row = 1,
                         right_annotation = ha)
dev.off()

##############################################################

