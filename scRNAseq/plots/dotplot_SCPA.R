##############

library(circlize)
library(ComplexHeatmap)
library(Seurat)
library(magrittr)
library(tidyverse)
library(msigdbr)
library(SCPA)
library(grid)
library(scales)


setwd("/beegfs/scratch/ric.cosr/dataset/160224")

stromacelltypess <- readRDS("stroma2_281023.RDS")

stromacelltypess$customclassif <- Idents(stromacelltypess)

DimPlot(stromacelltypess, group.by = "customclassif")

#stromacelltypess <- subset(x = stromacelltypess, idents = c("Unknown"), invert = TRUE)

stromacelltypess$Cell_Type <- stromacelltypess$customclassif

DimPlot(stromacelltypess, group.by = "Cell_Type")

cell_types <- unique(stromacelltypess$Cell_Type)
atlas <- SplitObject(stromacelltypess, split.by = "stim")


pathways <- msigdbr("mouse", "H") %>%
  format_pathways()


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

table(stromacelltypess$Cell_Type)


meso <- scpa_out_all$`mesothelial cells`
meso$celltypes <- "mesothelial_cells"

Ly6c <- scpa_out_all$`Ly6c+  perivascular adventitial cells`
Ly6c$celltypes <- "Ly6c_peri_adv"

Cxcl9 <- scpa_out_all$`Cxcl9 positive TRC cells`
Cxcl9$celltypes <- "Cxcl9"

Red_Pulp_Reticular <- scpa_out_all$`Red Pulp Reticular Cells 1`
Red_Pulp_Reticular$celltypes <- "Red_Pulp_Reticular"

Pericytes1 <- scpa_out_all$`Pericytes 1`
Pericytes1$celltypes <- "Pericytes1"

T_B_reticular <- scpa_out_all$`T-B reticular cells`
T_B_reticular$celltypes <- "T_B_reticular"

Pericytes2 <- scpa_out_all$`Pericytes 2`
Pericytes2$celltypes <- "Pericytes2"

Proliferating <- scpa_out_all$Proliferating
Proliferating$celltypes <- "Proliferating"

FDCtype1cells <- scpa_out_all$`FDC type 1 cells`
FDCtype1cells$celltypes <- "FDCtype1cells"

Thy1 <- scpa_out_all$`Thy1+ TRCs`
Thy1$celltypes <- "Thy1_TRC"


celltypess <- rbind(meso,Ly6c,Cxcl9,Red_Pulp_Reticular,Pericytes1,
                    T_B_reticular,Pericytes2,Proliferating,FDCtype1cells,Thy1)
#celltypess <- rbind(celltypes_temp,icelltypes)

celltypess_sign <- celltypess[celltypess$adjPval < 0.05,]

celltypess_sign$log2FC <- log2(abs(celltypess_sign$FC))

table(celltypess_sign$celltypes)
library(dplyr)

#dat <- celltypess_sign %>% filter_all(any_vars(!is.na(.)))

#celltypess_sign[grep("HALLMARK_E2F_TARGETS", celltypess_sign$Pathway),]


pdf("dotPlot_CTRL_LEUA_log2fc_celltypess_new.pdf", 10, 15)
# Create dotplot
ggplot(celltypess_sign, aes(x = celltypes, y = reorder(Pathway, -FC), color = -FC, size = qval)) +
  geom_point() +
  scale_color_gradient2(high = muted("red"),
                        mid = "white",
                        low = muted("blue"),
                        midpoint = 0)  +  # Center color gradient at zero
  #  scale_size_continuous(range = c(1, 10)) +  # Adjust size range
  #  labs(x = "LeukB - LeukA ", y = "Pathway") +  # Update axis labels
  theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 60)) 
dev.off()

