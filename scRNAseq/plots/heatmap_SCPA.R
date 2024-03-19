#########################       
library(circlize)
library(ComplexHeatmap)
library(Seurat)
library(magrittr)
library(tidyverse)
library(msigdbr)
library(SCPA)
library(grid)
library(scales)
library(openxlsx)

setwd("/beegfs/scratch/ric.cosr/dataset/160224")

stromaCAFs <- readRDS("stroma2_CAFs_Califano_scTYPE.RDS")

idents <- c('apCAF','iCAF','myCAF')

Idents(stromaCAFs) <- stromaCAFs$customclassif

stromaCAFs <- subset(x = stromaCAFs, idents = c("Unknown"), invert = TRUE)

stromaCAFs$Cell_Type <- stromaCAFs$customclassif
################################################################################
pathways <- msigdbr("mouse", "H") %>%
  format_pathways()

stromaCAFs$Cell_Type <- stromaCAFs$customclassif

atlas <- SplitObject(stromaCAFs, split.by = "stim")


pathways <- msigdbr("mouse", "H") %>%
  format_pathways()

col_hm <- colorRamp2(colors = c("red", "white", "blue"), breaks = c(-6, 0, 6))


scpa_out_allz <- list()
for (i in cell_types) {
  
  CTRL <- seurat_extract(atlas$combined_ctrl, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  LEUB <- seurat_extract(atlas$LEUB, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_allz[[i]] <- compare_pathways(list(CTRL, LEUB), pathways, downsample = 2000, parallel = TRUE, cores = 4)
  # For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
  scpa_out_allz[[i]]$QvalFC <-  ifelse(scpa_out_allz[[i]]$FC<0 , (-scpa_out_allz[[i]]$qval),scpa_out_allz[[i]]$qval)
  scpa_out_allz[[i]] <- scpa_out_allz[[i]][c("Pathway", "QvalFC")]
  colnames(scpa_out_allz[[i]]) <- c("Pathway", paste(i, "QvalFC", sep = "_"))
}


scpa_out_allzbk <- scpa_out_allz

scpa_out_allz <- scpa_out_allz %>% 
  reduce(full_join, by = "Pathway") %>% 
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_QvalFC", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  column_to_rownames("Pathway")


pdf("heatmap_Hallmarks_CTRL_LEUB_celltypes_newheat.pdf", 8, 10)
Heatmap(as.matrix(scpa_out_allz),
        col = col_hm,
        heatmap_legend_param = list(at = c(-6, 0, 6)),
        name = "Qval",
        show_row_names = T,
        #right_annotation = row_an,
        column_names_gp = gpar(fontsize = 6),
        width = ncol(scpa_out_allz)*unit(3, "mm"), 
        height = nrow(scpa_out_allz)*unit(3, "mm"),
        row_names_gp = gpar(fontsize = 4),  
        border = T,
        #column_km = 2,
        row_km = 3,
        column_labels = idents)
dev.off()




scpa_out_all <- list()
for (i in cell_types) {
  
  CTRL <- seurat_extract(atlas$combined_ctrl, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  LEUB <- seurat_extract(atlas$LEUB, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_all[[i]] <- compare_pathways(list(CTRL, LEUB), pathways, downsample = 2000, parallel = TRUE, cores = 4)
  
  # For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
  
}

scpa_out_allbk <- scpa_out_all

write.xlsx(scpa_out_all, "heatmap_Hallmarks_CTRL_LEUB_celltypes_newheat.xlsx")
########################




################################################################################



scpa_out_allz <- list()
for (i in cell_types) {
  
  CTRL <- seurat_extract(atlas$combined_ctrl, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  LEUA <- seurat_extract(atlas$LEUA, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_allz[[i]] <- compare_pathways(list(CTRL, LEUA), pathways, downsample = 2000, parallel = TRUE, cores = 4)
  # For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
  scpa_out_allz[[i]]$QvalFC <-  ifelse(scpa_out_allz[[i]]$FC<0 , (-scpa_out_allz[[i]]$qval),scpa_out_allz[[i]]$qval)
  scpa_out_allz[[i]] <- scpa_out_allz[[i]][c("Pathway", "QvalFC")]
  colnames(scpa_out_allz[[i]]) <- c("Pathway", paste(i, "QvalFC", sep = "_"))
}


scpa_out_allzbk <- scpa_out_allz

scpa_out_allz <- scpa_out_allz %>% 
  reduce(full_join, by = "Pathway") %>% 
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_QvalFC", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  column_to_rownames("Pathway")


pdf("heatmap_Hallmarks_CTRL_LEUA_celltypes_newheat.pdf", 8, 10)
Heatmap(as.matrix(scpa_out_allz),
        col = col_hm,
        heatmap_legend_param = list(at = c(-6, 0, 6)),
        name = "Qval",
        show_row_names = T,
        #right_annotation = row_an,
        column_names_gp = gpar(fontsize = 6),
        width = ncol(scpa_out_allz)*unit(3, "mm"), 
        height = nrow(scpa_out_allz)*unit(3, "mm"),
        row_names_gp = gpar(fontsize = 4),  
        border = T,
        #column_km = 2,
        row_km = 3,
        column_labels = idents)
dev.off()



scpa_out_all <- list()
for (i in cell_types) {
  
  CTRL <- seurat_extract(atlas$combined_ctrl, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  LEUA <- seurat_extract(atlas$LEUA, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_all[[i]] <- compare_pathways(list(CTRL, LEUA), pathways, downsample = 2000, parallel = TRUE, cores = 4)
  
  # For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
  
}

scpa_out_allbk <- scpa_out_all

write.xlsx(scpa_out_all, "heatmap_Hallmarks_CTRL_LEUA_celltypes_newheat.xlsx")
########################




################################################################################


scpa_out_allz <- list()
for (i in cell_types) {
  
  LEUB <- seurat_extract(atlas$LEUB, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  LEUA <- seurat_extract(atlas$LEUA, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_allz[[i]] <- compare_pathways(list(LEUB, LEUA), pathways, downsample = 2000, parallel = TRUE, cores = 4)
  # For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
  scpa_out_allz[[i]]$QvalFC <-  ifelse(scpa_out_allz[[i]]$FC<0 , (-scpa_out_allz[[i]]$qval),scpa_out_allz[[i]]$qval)
  scpa_out_allz[[i]] <- scpa_out_allz[[i]][c("Pathway", "QvalFC")]
  colnames(scpa_out_allz[[i]]) <- c("Pathway", paste(i, "QvalFC", sep = "_"))
}


scpa_out_allzbk <- scpa_out_allz

scpa_out_allz <- scpa_out_allz %>% 
  reduce(full_join, by = "Pathway") %>% 
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_QvalFC", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  column_to_rownames("Pathway")


pdf("heatmap_Hallmarks_LEUB_LEUA_celltypes_newheat.pdf", 8, 10)
Heatmap(as.matrix(scpa_out_allz),
        col = col_hm,
        heatmap_legend_param = list(at = c(-6, 0, 6)),
        name = "Qval",
        show_row_names = T,
        #right_annotation = row_an,
        column_names_gp = gpar(fontsize = 6),
        width = ncol(scpa_out_allz)*unit(3, "mm"), 
        height = nrow(scpa_out_allz)*unit(3, "mm"),
        row_names_gp = gpar(fontsize = 4),  
        border = T,
        #column_km = 2,
        row_km = 3,
        column_labels = idents)
dev.off()



scpa_out_all <- list()
for (i in cell_types) {
  
  LEUB <- seurat_extract(atlas$LEUB, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  LEUA <- seurat_extract(atlas$LEUA, 
                         meta1 = "Cell_Type", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_all[[i]] <- compare_pathways(list(LEUB, LEUA), pathways, downsample = 2000, parallel = TRUE, cores = 4)
  
  # For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
  
}


write.xlsx(scpa_out_all, "heatmap_Hallmarks_LEUB_LEUA_celltypes_newheat.xlsx")
################################################################################
