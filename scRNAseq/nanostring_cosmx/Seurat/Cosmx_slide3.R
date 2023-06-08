# libraries ---------------------------------------------------------------
library(Seurat)
library(future)
library(tidyverse)
library(ragg)
plan("multisession", workers = 10)




obj = readRDS("nano.obj_filter_tris.RDS")

DimPlot(obj, label = T)


ImageDimPlot(obj, fov = "R5630_Slide3", cols = "polychrome", axes = TRUE)+coord_cartesian(ylim = c(300000, 290000), xlim = c(-30000,-45000 ))


ImageDimPlot(obj, fov = "R5630_Slide3", axes = TRUE, cols = "polychrome")

p2 <- ImageDimPlot(obj, fov = "R5630_Slide3", alpha = 0.3, molecules = "RGS5", nmols = 10000) +
  NoLegend()


p1 <- ImageFeaturePlot(obj, fov = "R5630_Slide3", features = "CD68", max.cutoff = "q95")

ImageDimPlot(nano.obj, fov = "lung5.rep1", axes = TRUE, cols = "glasbey")



# -------------------------------------------------------------------------
# read in the data from Auroro
nano.obj <- LoadNanostring(data.dir = "/beegfs/scratch/ric.cosr/ric.cosr/DellabonaP/CosMX/5_Raw_data/R5630_Slide3", fov = "R5630_Slide3")



# data.dir <- "/beegfs/scratch/ric.cosr/ric.cosr/DellabonaP/CosMX/5_Raw_data/R5630_Slide3"
# test2 <- ReadNanostring(data.dir = data.dir, type = "centroids")
# cents <- CreateCentroids(test2$centroids)
# coords <- CreateFOV(coords = list("centroids" = cents), type = "centroids")
# obj <- CreateSeuratObject(counts = test2$matrix)
# obj[["fov"]]<-subset(coords, cell=Cells(obj))
# 
# str(obj@images$fov)



# filter the barcodes based on the final object
meta_input <- nano.obj@meta.data %>%
  rownames_to_column("barcode")



# read in the original object and filter for the metadata
test <- readRDS("/beegfs/scratch/ric.cosr/ric.cosr/DellabonaP/CosMX/Cosmx_seurat_object.Rds")
meta_ref <- test@meta.data %>% 
  filter(Run_Tissue_name == "R5630_Slide3") %>%
  rownames_to_column("barcode") %>% 
  mutate(barcode_fix=str_remove_all(barcode,pattern = "R5630.Slide3_"))



meta_full <- left_join(meta_input,meta_ref,by=c("barcode"="barcode_fix")) %>% 
  mutate(keep = case_when(is.na(Width)~0,
                          T~1)) %>% 
  column_to_rownames("barcode")



#
nano.obj@meta.data <- meta_full



nano.obj_filter <- subset(nano.obj,subset = keep == 1)




# standard preprocessing
nano.obj_filter <- SCTransform(nano.obj_filter, assay = "Nanostring") %>%
  RunPCA(npcs = 30, features = rownames(nano.obj_filter)) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.3)



# sample plotting
nano.obj_filter@meta.data
DimPlot(nano.obj_filter, reduction = "umap",label = T,group.by = "nb_clus")


DimPlot(nano.obj_filter, reduction = "umap",label = T,split.by = "nb_clus")
ImageDimPlot(nano.obj_filter, fov = "R5630_Slide3", cols = "polychrome", axes = TRUE)
ImageDimPlot(nano.obj_filter, fov = "R5630_Slide3", cols = "polychrome", axes = TRUE,group.by = "nb_clus")



ImageDimPlot(nano.obj_filter, fov = "R5630_Slide3", cols = "polychrome", axes = TRUE,group.by = "fov")
ImageDimPlot(nano.obj_filter, fov = "R5630_Slide3", cols = "polychrome", axes = TRUE)+coord_cartesian(ylim = c(300000, 290000), xlim = c(-50000,0))



dim(nano.obj_filter@meta.data)
# nano.obj@images$R5630_Slide3@boundaries$centroids
nano.obj_filter$responder <- case_when(nano.obj_filter@meta.data$fov %in% c(1,2)~"NR",
                                       T~"CR")



nano.obj_filter$cellType.responder <- paste0(nano.obj_filter$nb_clus,"_",nano.obj_filter$responder)



nano.obj_filter@meta.data$nb_clus %>% table()



#saveRDS(nano.obj_filter, "nano.obj_filter_slide3.RDS")

Idents(nano.obj_filter) <- "cellType.responder"



plasmablast <- FindMarkers(nano.obj_filter,ident.1 = "plasmablast_CR",ident.2 = "plasmablast_NR",logfc.threshold = 0)
plasmablast$gene <- rownames(plasmablast)

fibroblast <- FindMarkers(nano.obj_filter,ident.1 = "fibroblast_CR",ident.2 = "fibroblast_NR",logfc.threshold = 0)
fibroblast$gene <- rownames(fibroblast)


endothelial <- FindMarkers(nano.obj_filter,ident.1 = "endothelial_CR",ident.2 = "endothelial_NR",logfc.threshold = 0)
endothelial$gene <- rownames(endothelial)

a <- FindMarkers(nano.obj_filter,ident.1 = "a_CR",ident.2 = "a_NR",logfc.threshold = 0)
a$gene <- rownames(a)


b <- FindMarkers(nano.obj_filter,ident.1 = "b_CR",ident.2 = "b_NR",logfc.threshold = 0)
b$gene <- rownames(b)

c <- FindMarkers(nano.obj_filter,ident.1 = "c_CR",ident.2 = "c_NR",logfc.threshold = 0)
c$gene <- rownames(c)


d <- FindMarkers(nano.obj_filter,ident.1 = "d_CR",ident.2 = "d_NR",logfc.threshold = 0)
d$gene <- rownames(d)


e <- FindMarkers(nano.obj_filter,ident.1 = "e_CR",ident.2 = "e_NR",logfc.threshold = 0)
e$gene <- rownames(e)

f <- FindMarkers(nano.obj_filter,ident.1 = "f_CR",ident.2 = "f_NR",logfc.threshold = 0)
f$gene <- rownames(f)


B_cell <- FindMarkers(nano.obj_filter,ident.1 = "B-cell_CR",ident.2 = "B-cell_NR",logfc.threshold = 0)
B_cell$gene <- rownames(B_cell)


macrophage <- FindMarkers(nano.obj_filter,ident.1 = "macrophage_CR",ident.2 = "macrophage_NR",logfc.threshold = 0)
macrophage$gene <- rownames(macrophage)


mast <- FindMarkers(nano.obj_filter,ident.1 = "mast_CR",ident.2 = "mast_NR",logfc.threshold = 0)
mast$gene <- rownames(mast)

mDC <- FindMarkers(nano.obj_filter,ident.1 = "mDC_CR",ident.2 = "mDC_NR",logfc.threshold = 0)
mDC$gene <- rownames(mDC)


monocyte <- FindMarkers(nano.obj_filter,ident.1 = "monocyte_CR",ident.2 = "monocyte_NR",logfc.threshold = 0)
monocyte$gene <- rownames(monocyte)


NK <- FindMarkers(nano.obj_filter,ident.1 = "NK_CR",ident.2 = "NK_NR",logfc.threshold = 0)
NK$gene <- rownames(NK)


pDC <- FindMarkers(nano.obj_filter,ident.1 = "pDC_CR",ident.2 = "pDC_NR",logfc.threshold = 0)
pDC$gene <- rownames(pDC)


TCD4_memory <- FindMarkers(nano.obj_filter,ident.1 = "T CD4 memory_CR",ident.2 = "T CD4 memory_NR",logfc.threshold = 0)
TCD4_memory$gene <- rownames(TCD4_memory)


TCD4_naive <- FindMarkers(nano.obj_filter,ident.1 = "T CD4 naive_CR",ident.2 = "T CD4 naive_NR",logfc.threshold = 0)
TCD4_naive$gene <- rownames(TCD4_naive)


TCD8_naive <- FindMarkers(nano.obj_filter,ident.1 = "T CD8 naive_CR",ident.2 = "T CD8 naive_NR",logfc.threshold = 0)
TCD8_naive$gene <- rownames(TCD8_naive)


TCD8_memory <- FindMarkers(nano.obj_filter,ident.1 = "T CD8 memory_CR",ident.2 = "T CD8 memory_NR",logfc.threshold = 0)
TCD8_memory$gene <- rownames(TCD8_memory)


Treg <- FindMarkers(nano.obj_filter,ident.1 = "Treg_CR",ident.2 = "Treg_NR",logfc.threshold = 0)
Treg$gene <- rownames(Treg)


filename_xls <- 'slide3_plasmablast_FindMarkers.xlsx'
write.xlsx(plasmablast,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_fibroblast_FindMarkers.xlsx'
write.xlsx(fibroblast,
           file= filename_xls, 
           rownames = F,
           asTable = T)



filename_xls <- 'slide3_endothelial_FindMarkers.xlsx'
write.xlsx(endothelial,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_a_FindMarkers.xlsx'
write.xlsx(a,
           file= filename_xls, 
           rownames = F,
           asTable = T)

filename_xls <- 'slide3_Bcell_FindMarkers.xlsx'
write.xlsx(B_cell,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_b_FindMarkers.xlsx'
write.xlsx(b,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_c_FindMarkers.xlsx'
write.xlsx(c,
           file= filename_xls, 
           rownames = F,
           asTable = T)



filename_xls <- 'slide3_d_FindMarkers.xlsx'
write.xlsx(d,
           file= filename_xls, 
           rownames = F,
           asTable = T)




filename_xls <- 'slide3_e_FindMarkers.xlsx'
write.xlsx(e,
           file= filename_xls, 
           rownames = F,
           asTable = T)

filename_xls <- 'slide3_f_FindMarkers.xlsx'
write.xlsx(f,
           file= filename_xls, 
           rownames = F,
           asTable = T)




filename_xls <- 'slide3_macrophage_FindMarkers.xlsx'
write.xlsx(macrophage,
           file= filename_xls, 
           rownames = F,
           asTable = T)



filename_xls <- 'slide3_mast_FindMarkers.xlsx'
write.xlsx(mast,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_mDC_FindMarkers.xlsx'
write.xlsx(mDC,
           file= filename_xls, 
           rownames = F,
           asTable = T)



filename_xls <- 'slide3_monocyte_FindMarkers.xlsx'
write.xlsx(monocyte,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_NK_FindMarkers.xlsx'
write.xlsx(NK,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_pDC_FindMarkers.xlsx'
write.xlsx(pDC,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_plasmablast_FindMarkers.xlsx'
write.xlsx(plasmablast,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_TCD4memory_FindMarkers.xlsx'
write.xlsx(TCD4_memory,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_TCD8memory_FindMarkers.xlsx'
write.xlsx(TCD8_memory,
           file= filename_xls, 
           rownames = F,
           asTable = T)



filename_xls <- 'slide3_TCD8naive_FindMarkers.xlsx'
write.xlsx(TCD8_naive,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_TCD4naive_FindMarkers.xlsx'
write.xlsx(TCD4_naive,
           file= filename_xls, 
           rownames = F,
           asTable = T)


filename_xls <- 'slide3_Treg_FindMarkers.xlsx'
write.xlsx(Treg,
           file= filename_xls, 
           rownames = F,
           asTable = T)











FindMarkers(object = nano.obj_filter)



Idents(nano.obj_filter) <- "nb_clus"
test_markers <- FindAllMarkers(nano.obj_filter,verbose = T)


library(openxlsx)

filename_xls <- 'slide2_Findallmarkers.xlsx'
write.xlsx(test_markers,
           file= filename_xls, 
           rownames = F,
           asTable = T)




FindMarkers(object = nano.obj_filter)



test_markers %>% 
  filter(cluster=="endothelial")



FeaturePlot(nano.obj_filter,features = c("RGS5","MYH11","PECAM1","VWF"),order = T)



table(nano.obj_filter@meta.data$responder,nano.obj_filter@meta.data$fov)



test_plasmablast <- subset(nano.obj_filter,subset =  nb_clus=="plasmablast")



DimPlot(test_plasmablast, reduction = "umap",label = T,group.by = "fov")




# sample cropping ---------------------------------------------------------
# create a Crop
cropped.coords <- Crop(nano.obj_filter[["R5630_Slide3"]], y = c(300000, 291000), x = c(-42500, -37500), coords = "plot")
# set a new field of view (fov)
nano.obj_filter[["test"]] <- cropped.coords



# visualize FOV using default settings (no cell boundaries)
p1 <- ImageDimPlot(nano.obj, fov = "test", axes = TRUE, size = 0.7, border.color = "white", cols = "polychrome",
                   coord.fixed = FALSE)



# visualize FOV with full cell segmentations
DefaultBoundary(nano.obj[["test"]]) <- "segmentation"
p2 <- ImageDimPlot(nano.obj, fov = "test", axes = TRUE, border.color = "white", border.size = 0.1,
                   cols = "polychrome", coord.fixed = FALSE)



# simplify cell segmentations
nano.obj[["test"]][["simplified.segmentations"]] <- Simplify(coords = nano.obj[["test"]][["segmentation"]],
                                                             tol = 3)
DefaultBoundary(nano.obj[["test"]]) <- "simplified.segmentations"



# visualize FOV with simplified cell segmentations
DefaultBoundary(nano.obj[["test"]]) <- "simplified.segmentations"
p3 <- ImageDimPlot(nano.obj, fov = "test", axes = TRUE, border.color = "white", border.size = 0.1,
                   cols = "polychrome", coord.fixed = FALSE)



p1 + p2 + p3





# -------------------------------------------------------------------------



DefaultFOV(test)



meta_final %>% 
  filter(Run_Tissue_name=="R5630_Slide3")
