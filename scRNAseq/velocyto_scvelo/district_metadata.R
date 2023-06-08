#open from Rstudio with 64G
#srun -t 08:00:00 --mem=64g --pty bash


#suppressPackageStartupMessages(library(monocle3)) # OK
#suppressPackageStartupMessages(library(SeuratWrappers)) # OK
#suppressPackageStartupMessages(library(Seurat)) # OK
#suppressPackageStartupMessages(library(SeuratData)) # OK
#suppressPackageStartupMessages(library(SeuratWrappers)) # OK
#suppressPackageStartupMessages(library(ggplot2)) # OK
#suppressPackageStartupMessages(library(patchwork)) # OK
#suppressPackageStartupMessages(library(magrittr)) # OK
#suppressPackageStartupMessages(library(future)) # OK
#suppressPackageStartupMessages(library(cowplot)) # OK
#suppressPackageStartupMessages(library(dplyr)) # OK

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))

setwd('/beegfs/scratch/ric.cosr/ric.brendolan/BrendolanA_1280_scRNAseq/dataset/090521/Cd34subset')


all.combined = readRDS(file = "integratedCd34.Rds")

#extract single integrated objects

mLN <- subset(x = all.combined, subset = district == "mLN")

iLN <- subset(x = all.combined, subset = district == "iLN")

brLN <- subset(x = all.combined, subset = district == "brLN")


# generating meta file with all info

meta_data_general <- cbind(rownames(mLN@meta.data),as.data.frame(Idents(mLN)),as.data.frame(mLN@meta.data$development),as.data.frame(mLN[["umap"]]@cell.embeddings) )
colnames(meta_data_general)<- c("Cell","cell_type","dev_state","UMAP_1","UMAP_2")
write.csv(meta_data_general, file = "meta_data_dev_mLN_seurat.csv")


meta_data_general <- cbind(rownames(iLN@meta.data),as.data.frame(Idents(iLN)),as.data.frame(iLN@meta.data$development),as.data.frame(iLN[["umap"]]@cell.embeddings) )
colnames(meta_data_general)<- c("Cell","cell_type","dev_state","UMAP_1","UMAP_2")
write.csv(meta_data_general, file = "meta_data_dev_iLN_seurat.csv")


meta_data_general <- cbind(rownames(brLN@meta.data),as.data.frame(Idents(brLN)),as.data.frame(brLN@meta.data$development),as.data.frame(brLN[["umap"]]@cell.embeddings) )
colnames(meta_data_general)<- c("Cell","cell_type","dev_state","UMAP_1","UMAP_2")
write.csv(meta_data_general, file = "meta_data_dev_brLN_seurat.csv")







#later remove suffix from cells due to integration




Scrivo qui update sulle cellule isolate con i nuovi marcatori specifici per bladder cancer validati da Doglioni tra quelli proposti da me e Sofia.
Come anticipavo lo scorso meeting le uniche cellule CNV+ tra quelle isolate dai pazienti 271022, 081122, e 091122 sono: BLAD271022_A3(PERCP/EPCAM/CLAU/CD24),
BLAD091122_A3(CD24 PICCOLE), BLAD091122_C3 (CD24 PICCOLE), URIN091122_E3 (EPCAM/CLAU/CD318/CD24), URIN091122_A5(DAPI). Quindi 1 cellula nel 271022, nessuna
cellula nel 081122 e 4 cellule nel 091122. Gianna poi diceva che vediamo cellule positive a questi marcatori anche nei controlli sani.. un po' come
Succedeva con EpCAMneg/Clau/Vim ... Peccato.
