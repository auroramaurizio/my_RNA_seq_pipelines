
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

TU.subset = readRDS("EC_T.Rds")

TU.subset$CellTypes = Idents(TU.subset)

DefaultAssay(TU.subset) = "RNA"

integrated_3TIP <-readRDS("integrated_EC.Rds")


DefaultAssay(integrated_3TIP) = "integrated"
res= 1
sub <- FindClusters(integrated_3TIP, resolution = res)
DimPlot(sub, split.by = "stim")

subs = sub

DimPlot(subs, split.by = "stim")

new.cluster.ids.lit <- c('1',
                         '2',
                         '3'
)

names(new.cluster.ids.lit) <- levels(subs)
integrated_3TIP_new_eight <- RenameIdents(subs, new.cluster.ids.lit)

integrated_3TIP_new_eight$CellTypes = Idents(integrated_3TIP_new_eight)

integrated_3TIP_new_eight <- RunUMAP(integrated_3TIP_new_eight, reduction = "pca", dims = 1:35, verbose = FALSE, return.model=TRUE)

DimPlot(integrated_3TIP_new_eight)


TU.anchors <- FindTransferAnchors(reference = integrated_3TIP_new_eight, query = TU.subset,
                                       dims = 1:30)


predictions <- TransferData(anchorset = TU.anchors, refdata = integrated_3TIP_new_eight$CellTypes,
                            dims = 1:30)



predic_fil = predictions[predictions$predicted.id == "3", ]

summary(predic_fil$prediction.score.max)
TU.subset <- AddMetaData(TU.subset, metadata = predictions)


TU.subset <- MapQuery(anchorset = TU.anchors, reference = integrated_3TIP_new_eight, query = TU.subset,
                                 refdata = list(celltype = "CellTypes"), reference.reduction = "pca", reduction.model = "umap")


p1 <- DimPlot(integrated_3TIP_new_eight, reduction = "umap",
              label = TRUE, label.size = 3,
              repel = TRUE,  cols = c('1' = '#0066FF',
                                      '2' = '#336666',
                                      '3' = '#399933'),
              pt.size = 4) + NoLegend() + ggtitle("Reference annotations")



p2 <- DimPlot(TU.subset, reduction = "ref.umap", 
              group.by = "predicted.celltype", 
              label = TRUE, label.size = 3, 
              repel = TRUE,  cols = c('1' = '#0066FF',
                                      '2' = '#336666',
                                      '3' = '#399933'),
             pt.size = 4) + NoLegend() + ggtitle("Query transferred labels Toma Uninj")


pdf("label_tarnsfer.pdf", 12, 10)
p1+p2
dev.off()

DimPlot(TU.anchors)


