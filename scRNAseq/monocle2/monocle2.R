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

###################################################################

#saveRDS(integrated, file = "integrated_subset_clust_renamed_30pc.Rds")

setwd("/beegfs/scratch/ric.cosr/ric.brendolan/BrendolanA_1280_scRNAseq/dataset/090521/Cd34subset/br")

#load seurat object
integrated = readRDS(file ="integratedCd34.Rds")

#select assay
DefaultAssay(integrated) = "integrated"


#########################################
#########################################
#########################################

# select the district
br.combined <- subset (integrated,  subset = `district` == 'brLN')

#Extract data, phenotype data, and feature data from the SeuratObject --> here we decided to work on the scaled data
data <- as(as.matrix(br.combined@assays$integrated@scale.data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = br.combined@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
br_monocle2 <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              #lowerDetectionLimit = 0.5,
                              expressionFamily = uninormal())# since I have already normalized, thresholded and scaled data

#extract variable genes
var_genes <- br.combined[["integrated"]]@var.features
ordering_genes <- var_genes

br_monocle2 <- setOrderingFilter(br_monocle2, ordering_genes)

#we decided to consider 2 PC. 
br_monocle2 <- reduceDimension(br_monocle2, max_components = 2,
                               method = 'DDRTree', norm_method="none", pseudo_expr=0,scaling=TRUE)

## order cells change colors and theta to match your plot
br_monocle2 <- orderCells(br_monocle2)

#plot the trajectories colored by celltype, developmental stage, state, pseudotime etc.

P1 = plot_cell_trajectory(br_monocle2, color_by = "CellType", cell_name_size = 8) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20, face="bold"))


#help(plot_cell_trajectory)

pdf("plot_cell_trajectory_CellType_br_subset_2pc_Cd34_scaled.pdf") 
P1
dev.off()

P2 = plot_cell_trajectory(br_monocle2, color_by = "development") + scale_color_manual(values=c('embryo' = '#00BA38', 'P2' = '#619CFF',
                                                                                           'adult' = '#F8766D'))+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))


pdf("plot_cell_trajectory_development_br_subset_2pc_Cd34_scaled.pdf")
P2
dev.off()

P3 = plot_cell_trajectory(br_monocle2, color_by = "State")+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))


pdf("plot_cell_trajectory_state_br_subset_2pc_Cd34_scaled.pdf")
P3
dev.off()

P3_bis = plot_cell_trajectory(br_monocle2, color_by = "State") + facet_wrap(~development, nrow = 1) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))

pdf("plot_cell_trajectory_state_facet_wrap_br_subset_2pc_Cd34_scaled.pdf", 30, 10)
P3_bis
dev.off()


Ptest= plot_cell_trajectory(test, color_by = "State") + facet_wrap(~development, nrow = 1) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))

pdf("ptest.pdf", 30, 10)
Ptest
dev.off()


P4 = plot_cell_trajectory(br_monocle2, color_by = "CellType") + facet_wrap(~development, nrow = 1)+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))



pdf("plot_cell_trajectory_facetwrap_dev_CellType_br_subset_2pc_Cd34_scaled.pdf", 30, 10) 
P4
dev.off()


P4_bis = plot_cell_trajectory(br_monocle2, color_by = "CellType") + facet_wrap(~CellType, nrow = 1)+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))



pdf("plot_cell_trajectory_facetwrap_CellType_CellType_br_subset_2pc_Cd34_scaled.pdf", 40, 10) 
P4_bis
dev.off()

#use this function to place the origin of the trajectory in cluster (CellType) 4, the one with the highest concentration on prolif cells.

GM_state <- function(cds){
  if (length(unique(cds@phenoData@data$State)) > 1){
    T0_counts <- table(cds@phenoData@data$State, cds@phenoData@data$CellType)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}


br_monocle2 <- orderCells(br_monocle2, root_state = GM_state(br_monocle2))


P6 = plot_cell_trajectory(br_monocle2, color_by = "Pseudotime") + facet_wrap(~development, nrow = 1)+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))

pdf("plot_cell_trajectory_facetwrap_development_Pseudotime_br_2pc_Cd34_scaled.pdf", 30, 10)
P6
dev.off()


P6bis = plot_cell_trajectory(br_monocle2, color_by = "Pseudotime") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))

pdf("plot_cell_trajectory_development_Pseudotime_br_2pc_Cd34_scaled.pdf")
P6bis
dev.off()

#plot_multiple_branches_heatmap accepts only root and terminal branches, no intermediate ones

plot_table <- plot_multiple_branches_heatmap(br_monocle2,
                                             branches = c(1,2,6),
                                             branches_name = c("1","2","6"),
                                             cores = 32,
                                             return_heatmap = TRUE,
                                             use_gene_short_name = T,
                                             show_rownames = T)


#saveRDS(plot_table, "plot_table.RDS")

#saveRDS(br_monocle2, "br_monocle2.RDS")

#write.table(plot_table$annotation_row, 'plot_newheatmap_br_subset_2pc_bp2_Cd34.txt')

#head(plot_table$tree_row)

pdf("plot_genes_branched_heatmap_br_Cd34.pdf", 20, 120)
plot_multiple_branches_heatmap(br_monocle2,
                               branches = c(1,2,6),
                               branches_name = c("1","2","6"),
                               cores = 32,
                               use_gene_short_name = T,
                               show_rownames = T)
dev.off()

saveRDS(i_monocle2, "i_monocle2.RDS")

test = readRDS("i_monocle2.RDS")

help(plot_multiple_branches_heatmap)

pdf("plot_genes_branched_heatmap_br_Cd34_small_1module.pdf")
plot_multiple_branches_heatmap(br_monocle2,
                               branches = c(1,2,6),
                               branches_name = c("1","2","6"),
                               cores = 16,
                               num_clusters = 1,
                               use_gene_short_name = T,
                               show_rownames = F)
dev.off()



pdf("plot_genes_branched_heatmap_br_Cd34_small_6modules.pdf")
plot_multiple_branches_heatmap(br_monocle2,
                               branches = c(1,2,6),
                               branches_name = c("1","2","6"),
                               cores = 16,
                               num_clusters = 6,
                               use_gene_short_name = T,
                               show_rownames = F)
dev.off()



br_monocle2_res <- BEAM(br_monocle2, branch_point = 1, cores = 16)
br_monocle2_res <- br_monocle2_res[order(br_monocle2_res$qval),]
br_monocle2_res <- br_monocle2_res[,c("gene_short_name", "pval", "qval")]



plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp1_Cd34 <- plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                                                                                       qval < 1e-4)),],
                                                                                          branch_point = 1,
                                                                                          cores = 16,
                                                                                          num_clusters = 6,
                                                                                          return_heatmap = TRUE,
                                                                                          use_gene_short_name = T,
                                                                                          show_rownames = F)


write.table(plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp1_Cd34$annotation_row, 'plot_heatmap_br_subset_2pc_bp1_noname_Cd34_scaled.txt')


pdf("plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp1_Cd34_scaled.pdf", 20, 60)
plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                         qval < 1e-4)),],
                            branch_point = 1,
                            cores = 16,
                            num_clusters = 6,
                            show_rownames = T)
dev.off()


pdf("plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp1_noname_clust_Cd34_scaled.pdf")
plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                         qval < 1e-4)),],
                            branch_point = 1,
                            cores = 16,
                            num_clusters = 6,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()




pdf("plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp1_noname_noclust_Cd34_scaled.pdf")
plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                         qval < 1e-4)),],
                            branch_point = 1,
                            cores = 16,
                            num_clusters = 1,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()





br_monocle2_res <- BEAM(br_monocle2, branch_point = 2, cores = 16)
br_monocle2_res <- br_monocle2_res[order(br_monocle2_res$qval),]
br_monocle2_res <- br_monocle2_res[,c("gene_short_name", "pval", "qval")]


plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp2_Cd34 <- plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                         qval < 1e-4)),],
                            branch_point = 2,
                            cores = 16,
                            num_clusters = 6,
                            use_gene_short_name = T,
                            return_heatmap = TRUE,
                            show_rownames = T)



write.table(plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp2_Cd34$annotation_row, 'plot_heatmap_br_subset_2pc_bp2_Cd34_scaled.txt')


pdf("plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp2_Cd34.pdf", 20, 60)
plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                         qval < 1e-4)),],
                            branch_point = 2,
                            cores = 16,
                            num_clusters = 6,
                            show_rownames = T)
dev.off()


pdf("plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp2_noname_clust_Cd34_bis_scaled.pdf")
plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                         qval < 1e-4)),],
                            branch_point = 2,
                            cores = 16,
                            num_clusters = 6,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()



pdf("plot_cell_trajectory_plot_heatmap_br_subset_2pc_bp2_noname_noclust_Cd34_scaled.pdf")
plot_genes_branched_heatmap(br_monocle2[row.names(subset(br_monocle2_res,
                                                         qval < 1e-4)),],
                            branch_point = 2,
                            cores = 16,
                            num_clusters = 1,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()


##########
##########

