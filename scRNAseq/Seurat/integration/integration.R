suppressMessages(library(ggbeeswarm))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
library(pheatmap)
library(reshape2)
library(scales)
library(viridis)
library(SeuratWrappers)
library(slingshot)
require(BiocStyle)
library(SingleCellExperiment)
#print(version) R version 3.6.1
#Bonanomi
packageVersion("Seurat")
library(dittoSeq)

data_only_integrated_stim.Rds

data = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/public_data/data/data_only_integrated_stim.Rds")
DimPlot(data)
DefaultAssay(data)

Metadata = read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/public_data/data/GSE162610_barcode_metadata.tsv", sep = "\t")


Endo_Meta = Metadata[grep("CellsOfInterest", Metadata$celltype), ]
head(Endo_Meta)

cId_Endo_Meta = row.names(Endo_Meta)

DimPlot(data, label=T, cells.highlight=cId_Endo_Meta, cols.highlight = c("cyan"), cols= "grey")

data_EC <- subset(data, cells = cId_Endo_Meta)

DimPlot(data_EC, split.by = "stim", ncol = 2)

DefaultAssay(data_EC) = "RNA"


D0 <- subset(x = data_EC, subset = stim == "1_uninj")
D1 <- subset(x = data_EC, subset = stim == "1dpi")
D3 <- subset(x = data_EC, subset = stim == "3dpi")
D7 <- subset(x = data_EC, subset = stim == "7dpi")

object_clean_new1.list <- lapply(X = c(D0,D1,D3,D7), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# find integration anchors
integrated <- FindIntegrationAnchors(object.list = object_clean_new1.list, dims = 1:30)

features.to.integrate = integrated@anchor.features

integrated <- IntegrateData(anchorset = integrated, dims = 1:30, features.to.integrate = features.to.integrate)

DefaultAssay(integrated) <- "integrated"

res = 0.1

DefaultAssay(integrated) = "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = res)

saveRDS(integrated, "data_from_Metadata.RDS")

#get some info
table(integrated$orig.ident)
table(integrated$stim)


pdf("DimPlot_group_by_dev_split_by_stim.pdf", 12, 5)
DimPlot(integrated, reduction = "umap", group.by = "development", split.by = "stim")
dev.off()

pdf("DimPlot_group_by_dev_split_by_orig.ident.pdf", 12, 5)
DimPlot(integrated, reduction = "umap", group.by = "development", split.by = "orig.ident")
dev.off()


DefaultAssay(integrated) = "RNA"

# make a fatureplot of one or multiple genes
FeaturePlot(integrated, "Prox1", order = T)


some_genes= c("Ccl8", "Ccl24", "Cd209", "Clec10a", "Fcna", "Folr2", "Fxyd2", "Retnla", "Tslp")


Plot_sign <- function(Seraut.object, signature, operator = sum, titlename) {
  x <- Seraut.object
  DefaultAssay(x) <- "RNA"
  
  x[["Sign_exp"]] <- apply(FetchData(object = x, 
                                     vars = signature),
                           1,
                           operator)
  FP <- FeaturePlot(x, reduction = "tsne", 
                    features = 'Sign_exp', 
                    label = T, 
                    pt.size = 2,
                    order = T,
                    split.by="stim",
                    cols = c("lightgrey", "blue")) +
    theme(plot.title = element_text(color="blue", size=15, face="bold.italic"),
          plot.subtitle = element_text(color="dodgerblue2", size=8, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'dodgerblue4', size=10, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 10),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=10),
          axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 10),
          legend.text = element_text(face = "bold", color = "dodgerblue2", size = 10),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    #labs(title= titlename, subtitle = paste('',toString(signature), sep=''), 
    labs(title=titlename)
  #, 
  #     x = "tSNE 1", y = "tSNE 2") 
  return(FP)
}



Plot_sign(integrated,
          signature= some_genes, 
          operator = median, titlename = "some_genes")



#make a violinplot



Plot_sign_vln <- function(Seraut.object, signature, operator = sum,  titlename) {
  x <- Seraut.object
  DefaultAssay(x) <- "RNA"
  
  x[["Sign_exp"]] <- apply(FetchData(object = x, 
                                     vars = signature),
                           1,
                           operator)
  FP <- VlnPlot(x,
                features = 'Sign_exp', 
                split.by="stim",
                cols = c("limegreen","darkviolet"),
                split.plot = TRUE,
                pt.size = 0,
  ) +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=12, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'black', size=10, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 10),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=10),
          axis.title.y = element_text(face = "bold", color = "black", size = 10),
          #legend.text = element_text(face = "bold", color = "black", size = 6),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    #geom_boxplot(width=0.1, fill="white")+
    #labs(title= titlename, subtitle = paste('',toString(signature), sep=''), 
    labs(title=titlename)
  #, 
  #     x = "tSNE 1", y = "tSNE 2") 
  return(FP)
}



pdf("VlnPlot_pathway.pdf", 60, 20)
VlnPlot(integrated, 
        features = some_genes, 
        split.by = "stim", 
        cols = c("limegreen","darkviolet","pink","cyan"))
dev.off()

#make a dotplot

pdf("dotPlot_pathways.pdf", 60, 20)
DotPlot(integrated, 
        features = some_genes, 
        dot.scale = 8, split.by = "stim", 
        cols = c("limegreen","darkviolet","pink","cyan"), 
        assay = "RNA", idents = c("A","B", "C","D")) #+ 
#coord_flip()
dev.off()



#rename clusters
new.cluster.ids.lit <- c('celltype1', 
                         'celltype2',
                         'celltype3',
                         'celltype4')


names(new.cluster.ids.lit) <- levels(integratd)
integrated_new <- RenameIdents(integrated, new.cluster.ids.lit)

DimPlot(integrated_rn, cols = c('celltype1' = '#F8766D', 
                                'celltype2' = 'brown', 
                                'celltype3' = '#7CAE00', 
                                'celltype4' = '#CD9600'))



#make a barplot
pt <- table(Idents(integrated_new), integrated_new$stim)
pt <- as.data.frame(pt)
pt$Cluster <- as.character(pt$Var1)


ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") 








#make a heatmap

thresh.use = 0.25
min.pct = 0.25
min.diff.pct = -Inf
test.use ="wilcox"


DefaultAssay(integrated_new) <- "integrated"

cluster.markers = FindAllMarkers(integrated_new, thresh.use = thresh.use, test.use=test.use, min.pct=min.pct, min.diff.pct=min.diff.pct, only.pos=TRUE)

filename_xls <- 'FindAllMarkers.xlsx'
write.xlsx(cluster.markers,
           file= filename_xls, 
           row.names = F,
           asTable = T)


top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(integrated_new, features = top10$gene, size = 3) + NoLegend() + theme(text = element_text(size = 7))



DoHeatmap(integrated_new, features = y, size = 0, group.colors = c('celltype1' = '#F8766D', 
                                                               'celltype2' = 'brown', 
                                                               'celltype3' = '#7CAE00', 
                                                               'celltype4' = '#CD9600'), 
          angle = 90) + NoLegend() +  theme(text = element_text(size = 6))


reord=integrated_new
reord$orig.ident <- factor(x = reord$orig.ident , levels = c("celltype4", "celltype2", "celltype3","celltype1")
table(reord$orig.ident)



#extract from the integrated object only celltypes of interest

#sample2_bis <- subset(x = integrated_new, subset = stim == "2")

sample2_bis = subset(integrated_new, idents = c("celltype3", "celltype4"))





