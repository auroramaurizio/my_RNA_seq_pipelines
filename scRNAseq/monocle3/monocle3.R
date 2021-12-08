suppressPackageStartupMessages(library(monocle3)) # OK
suppressPackageStartupMessages(library(SeuratWrappers)) # NO
suppressPackageStartupMessages(library(Seurat)) # OK
suppressPackageStartupMessages(library(SeuratData)) # NO
suppressPackageStartupMessages(library(ggplot2)) # OK
suppressPackageStartupMessages(library(patchwork)) # OK
suppressPackageStartupMessages(library(magrittr)) # OK
suppressPackageStartupMessages(library(future)) # OK
suppressPackageStartupMessages(library(cowplot)) # OK
suppressPackageStartupMessages(library(dplyr)) # OK
#maybe you do not need all the above libraries. 

#suppressPackageStartupMessages(library(monocle)) 
#check monocle version 
print(version)
packageVersion("monocle3")
sessionInfo()

#import integrated seurat object 
#m.combined = ReadRDS...
setwd("/beegfs/scratch/ric.cosr/ric.brendolan/BrendolanA_1280_scRNAseq/dataset/090521/Cd34subset/m")

#load seurat object
integrated = readRDS(file ="integratedCd34.Rds")

#select assay
DefaultAssay(integrated) = "integrated"


#########################################
#########################################
#########################################

# select the district
m.combined <- subset (integrated,  subset = `district` == 'mLN')


#set assay to integrated
DefaultAssay(m.combined) = "integrated"
# create celldataset
cds.m <- as.cell_data_set(m.combined)

cds.m <- cluster_cells(cds.m, k = 20)

#p1 <- plot_cells(cds.m, cell_size = 0.9,show_trajectory_graph = FALSE) + ggtitle(label = "Clusters by Monocle3") + theme(
#  plot.title = element_text(color = "black", size = 20, face = "bold"))

#p2 <- plot_cells(cds.m, cell_size = 0.9,color_cells_by = "partition", show_trajectory_graph = FALSE) + ggtitle("Partitions (super clusters) by Monocle3") + theme(
#  plot.title = element_text(color = "black", size = 20, face = "bold"))

#wrap_plots(p1, p2)

cds.m <- learn_graph(cds.m,use_partition = TRUE)

plot_cells(cds.m, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)


# I define the root as the point with the highest density of "embryo" cells

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="embryo"){
  cell_ids <- which(colData(cds)[, "development"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


# I define the root as the point with the highest density of "embryo" cells
cds.m <- order_cells(cds.m, root_pr_nodes=get_earliest_principal_node(cds.m))



plot1 <- plot_cells(cds.m,
           color_cells_by = "development",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_roots = TRUE,
           cell_size = 0.9,
           label_branch_points=TRUE,
           graph_label_size=4) + ggtitle("mesenteric LN by development") + theme(
  plot.title = element_text(color = "black", size = 20, face = "bold"), legend.text = element_text(face = "bold", color = "black", size = 10),legend.title = element_text(face = "bold", color = "black", size = 12))



plot2 <- plot_cells(cds.m,
                    color_cells_by = "pseudotime",
                    label_cell_groups=FALSE,
                    label_leaves=FALSE,
                    cell_size = 0.9,
                    label_roots = TRUE,
                    label_branch_points=FALSE,
                    graph_label_size=4) + ggtitle("mesenteric LN by pseudotime") + theme(
                      plot.title = element_text(color = "black", size = 20, face = "bold"), legend.text = element_text(face = "bold", color = "black", size = 10),legend.title = element_text(face = "bold", color = "black", size = 12))


plot3 <- plot_cells(cds.m,
           color_cells_by = "CellType",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_roots = TRUE,
           cell_size = 0.9,
           label_branch_points=FALSE,
           graph_label_size=4) + ggtitle("mesenteric LN by development") + theme(
  plot.title = element_text(color = "black", size = 20, face = "bold"), legend.text = element_text(face = "bold", color = "black", size = 10),legend.title = element_text(face = "bold", color = "black", size = 12))



#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

#se vuoi valutare il comportamento di particolari geni lungo la traject usa la funzione plot_genes_in_pseudotime()

# se vuoi selezionare alcune cellule/branches e focalizzarti su quelle usa choose_cells()