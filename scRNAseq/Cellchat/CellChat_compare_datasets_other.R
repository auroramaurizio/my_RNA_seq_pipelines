#https://stackoverflow.com/questions/62372181/unavailable-to-install-systemfonts-packages-in-bioconductor-3-12-r-4-0-1
#install.packages("systemfonts", type = "source")
#install.packages("svglite")
#devtools::install_github("sqjin/CellChat")
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
suppressMessages(library(openxlsx))

setwd("/Users/maurizio.aurora")

#load R object
condition1 = readRDS("condition1.Rds")
DimPlot(condition1)

table(Idents(condition1))
#CellChat requires two user inputs: one is the gene expression data of cells, 
#and the other is either user assigned cell labels (i.e., label-based mode) or 
#a low-dimensional representation of the single-cell data (i.e., label-free mode). 
#For the latter, CellChat automatically groups cells by building a shared neighbor 
#graph based on the cell-cell distance in the low-dimensional space or the pseudotemporal trajectory space.

#For the gene expression data matrix, genes should be in rows with rownames and cells in columns with colnames. 
#Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) 
#is required as input for CellChat analysis. 
#If user provides count data, we provide a normalizeData function to account for library size and then do log-transformed. 
#For the cell group information, a dataframe with rownames is required as input for CellChat.
DefaultAssay(condition1) = "RNA"
#DimPlot()
# take raw data and normalise it. In seurat 3 raw.data are called counts
#count_raw_condition1 =as.matrix(GetAssayData(object = object, slot = "counts"))
#?normalizeData
#i can obtain the normalized data in this way with the normalize data function
#normalizeData(count_raw, scale.factor = 10000, do.log = TRUE)[1:5,1:5]

#or taking advantage of seurat's normalization data 
count_norm_condition1 = condition1@assays$RNA@data # normalized data matrix
#data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
#labels <- Idents(condition1)
#meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

meta_data_condition1 <- cbind(rownames(condition1@meta.data),as.data.frame(Idents(condition1)))
colnames(meta_data_condition1)<- c("Cell","labels")
cellchat_condition1 <- createCellChat(object = count_norm_condition1, meta = meta_data_condition1, group.by = "labels")
#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT
cellchat_condition1 <- addMeta(cellchat_condition1, meta = meta_data_condition1)
cellchat_condition1 <- setIdent(cellchat_condition1, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_condition1@idents) # show factor levels
groupSize <- as.numeric(table(cellchat_condition1@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat_condition1@DB <- CellChatDB.use
cellchat_condition1 <- subsetData(cellchat_condition1) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) #

cellchat_condition1 <- identifyOverExpressedGenes(cellchat_condition1)
cellchat_condition1 <- identifyOverExpressedInteractions(cellchat_condition1)

cellchat_condition1 <- projectData(cellchat_condition1, PPI.mouse)
cellchat_condition1 <- computeCommunProb(cellchat_condition1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_condition1 <- filterCommunication(cellchat_condition1, min.cells = 10)
cellchat_condition1 <- computeCommunProbPathway(cellchat_condition1)
cellchat_condition1 <- aggregateNet(cellchat_condition1)
groupSize <- as.numeric(table(cellchat_condition1@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("condition1_only_circleNofI.pdf")
netVisual_circle(cellchat_condition1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("condition1_only_circleSofI.pdf")
netVisual_circle(cellchat_condition1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat_condition1 <- cellchat_condition1@net$weight
par(mfrow = c(3,4), xpd=TRUE)

pdf("condition1_only.pdf")
for (i in 1:nrow(mat_condition1)) {
  mat2 <- matrix(0, nrow = nrow(mat_condition1), ncol = ncol(mat_condition1), dimnames = dimnames(mat_condition1))
  mat2[i, ] <- mat_condition1[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_condition1), title.name = rownames(mat_condition1)[i])
}
dev.off()


getwd()


############
############
############

getwd()
#load R object
condition2 = readRDS("condition2.Rds")
pdf("condition2.pdf")
DimPlot(condition2, split.by = "orig.ident", ncol = 2)
dev.off()
#CellChat requires two user inputs: one is the gene expression data of cells, 
#and the other is either user assigned cell labels (i.e., label-based mode) or 
#a low-dimensional representation of the single-cell data (i.e., label-free mode). 
#For the latter, CellChat automatically groups cells by building a shared neighbor 
#graph based on the cell-cell distance in the low-dimensional space or the pseudotemporal trajectory space.

#For the gene expression data matrix, genes should be in rows with rownames and cells in columns with colnames. 
#Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) 
#is required as input for CellChat analysis. 
#If user provides count data, we provide a normalizeData function to account for library size and then do log-transformed. 
#For the cell group information, a dataframe with rownames is required as input for CellChat.
DefaultAssay(condition2) = "RNA"
#DimPlot()
# take raw data and normalise it. In seurat 3 raw.data are called counts
#count_raw_condition2 =as.matrix(GetAssayData(object = object, slot = "counts"))
#?normalizeData
#i can obtain the normalized data in this way with the normalize data function
#normalizeData(count_raw, scale.factor = 10000, do.log = TRUE)[1:5,1:5]

#or taking advantage of seurat's normalization data 
count_norm_condition2 = condition2@assays$RNA@data # normalized data matrix
#data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
#labels <- Idents(condition2)
#meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

meta_data_condition2 <- cbind(rownames(condition2@meta.data),as.data.frame(Idents(condition2)))
colnames(meta_data_condition2)<- c("Cell","labels")
cellchat_condition2 <- createCellChat(object = count_norm_condition2, meta = meta_data_condition2, group.by = "labels")
#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT
cellchat_condition2 <- addMeta(cellchat_condition2, meta = meta_data_condition2)
cellchat_condition2 <- setIdent(cellchat_condition2, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_condition2@idents) # show factor levels
groupSize <- as.numeric(table(cellchat_condition2@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data
showDatabaseCategory(CellChatDB)


head(CellChatDB$cofactor)
filename_xls <- 'CellChatDB_interactionDB_human.xlsx'
write.xlsx(CellChatDB$interaction,
           file= filename_xls, 
           row.names = T,
           asTable = T)

(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat_condition2@DB <- CellChatDB.use
cellchat_condition2 <- subsetData(cellchat_condition2) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) #

cellchat_condition2 <- identifyOverExpressedGenes(cellchat_condition2)
cellchat_condition2 <- identifyOverExpressedInteractions(cellchat_condition2)

cellchat_condition2 <- projectData(cellchat_condition2, PPI.mouse)
cellchat_condition2 <- computeCommunProb(cellchat_condition2)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_condition2 <- filterCommunication(cellchat_condition2, min.cells = 10)
cellchat_condition2 <- computeCommunProbPathway(cellchat_condition2)
cellchat_condition2 <- aggregateNet(cellchat_condition2)
groupSize <- as.numeric(table(cellchat_condition2@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("condition2_only_circleNofI.pdf")
netVisual_circle(cellchat_condition2@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("condition2_only_circleSofI.pdf")
netVisual_circle(cellchat_condition2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
mat_condition2 <- cellchat_condition2@net$weight
par(mfrow = c(3,4), xpd=TRUE)
pdf("condition1_only.pdf")
for (i in 1:nrow(mat_condition2)) {
  mat2 <- matrix(0, nrow = nrow(mat_condition2), ncol = ncol(mat_condition2), dimnames = dimnames(mat_condition2))
  mat2[i, ] <- mat_condition2[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_condition2), title.name = rownames(mat_condition2)[i])
}
dev.off()
getwd()


getwd()


#########
#########
#########
saveRDS(cellchat_condition2, "cellchat_condition2.Rds")
saveRDS(cellchat_condition1, "cellchat_condition1.Rds")

data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)

object.list <- list(condition2 = cellchat_condition2, condition1 = cellchat_condition1)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

pdf("comparison_barplots.pdf")
gg1 + gg2
dev.off()

pdf("diffInteraction_net.pdf")
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf("diffInteraction_heat.pdf", 10, 5)
gg1 + gg2
dev.off()


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("diffInteraction_net_split.pdf")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()



for (i in 1:length(object.list)) {
  object.list[[i]] = netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP") 
}


for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("inter_2D_space.pdf", 10, 5)
patchwork::wrap_plots(plots = gg)
dev.off()

?netAnalysis_computeCentrality

cellchat_condition2 <- netAnalysis_computeCentrality(cellchat_condition2, slot.name = "netP") 
cellchat_condition1 <- netAnalysis_computeCentrality(cellchat_condition1, slot.name = "netP") 

?netAnalysis_signalingChanges_scatter
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MES_PERINEUR")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MES_PERINEUR", signaling.exclude = c("MES_TIP1"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))


#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("netVisual_embeddingPairwisefunctional.pdf")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()
#> 2D visualization of signaling networks from datasets 1 2

#Identify signaling groups based on structure similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("netVisual_embeddingPairwisestructural.pdf")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
dev.off()

#> 2D visualization of signaling networks from datasets 1 2

#Compute and visualize the pathway distance in the learned joint manifold
#We can identify the signaling networks with larger (or less) 
#difference based on their Euclidean distance in the shared two-dimensions space. 
#Larger distance implies larger difference of the communication networks between 
#two datasets in terms of either functional or structure similarity. 
#NB: We only compute the distance of overlapped signaling pathways between two datasets. 
#Those signaling pathways that are only identified in one dataset are not considered here. 
#If there are more than three datasets, one can do pairwise comparisons by defining 
#comparison in the function rankSimilarity.

pdf("rankSimilarity_functional.pdf")
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2
dev.off()


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

pdf("ranknet.pdf")
gg1 + gg2
dev.off()
#Identify and visualize the conserved and context-specific signaling pathways

#Compare outgoing (or incoming) signaling associated with each cell population

library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.7.1.1010
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, )
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16)

pdf("pathway.union_heat1.pdf",10,10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

pdf("pathway.union_heat2.pdf",10,10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "OrRd")
pdf("pathway.union_heat3.pdf",10,10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()


#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
pdf("netVisual_bubble_mergedobject.pdf",15,20)
netVisual_bubble(cellchat, sources.use = "MES_PERINEUR", targets.use = c(1:15),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object
dev.off()



gg1 <- netVisual_bubble(cellchat, sources.use = "MES_PERINEUR", targets.use = c(1:15),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = "MES_PERINEUR", targets.use = c(1:15),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("netVisual_bubble_mergedobject_inc_dec.pdf",15,20)
gg1 + gg2
dev.off()


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "condition1"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "condition1",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "condition2",ligand.logFC = -0.1, receptor.logFC = -0.1)

#Identify dysfunctional signaling by using differential expression analysis

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)



pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object

pdf("bubbleplot_up_dw_regulated_signaling_LRpairs.pdf",12,12)
gg1 + gg2
dev.off()

# Chord diagram
pdf("chord_up_dw_regulated_signaling_LRpairs.pdf",12,12)
netVisual_chord_gene(object.list[[2]], sources.use = "MES_PERINEUR", targets.use = c(1:15), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = "MES_PERINEUR", targets.use = c(1:15), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
#Part IV: Visually compare cell-cell communication using Hierarchy plot,
#Circle plot or Chord diagram



pathways.show <- c("VEGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("VEGF_graph.pdf", 15, 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()




pathways.show <- c("VEGF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
pdf("VEGF_heatmap.pdf")
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()


# Chord diagram
pathways.show <- c("VEGF") 
par(mfrow = c(2,2), xpd=TRUE)

#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

?netVisual_aggregate


pathways.show <- c("VEGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 


pdf("prova.pdf")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(object.list[[i]], signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy" )
dev.off()

?netVisual_aggregate

pdf("hierarchy.pdf", 15, 10)
for (i in 1:length(object.list)) {
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = c(2,3,4,5,6,7,9,10,11,12,13,14,15) # a numeric vector. 
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "hierarchy")
}
dev.off()

pdf("contribution_of_each_ligand-receptor_pair_to_the_overall_signaling_pathway.pdf", 15, 10)
for (i in 1:length(object.list)) {
  # Hierarchy plot
  netAnalysis_contribution(object.list[[2]], signaling = pathways.show)
}
dev.off()




netAnalysis_contribution(object.list[[i]], signaling = pathways.show)

head(object.list[[i]])

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

?netVisual_aggregate

pdf(file ="VEGF_Chorddiagram.pdf", width = 20, height =16)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

?netVisual_chord_gene


pdf(file ="VEGF_ChorddiagramGene_EPINEUR_to_all.pdf", width = 20, height =16)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = "MES_EPINEUR", targets.use = c(1:15), lab.cex = 0.5, title.name = paste0("Signaling from Perineurial cells to all", names(object.list)[i]))
}
dev.off()

pdf(file ="VEGF_ChorddiagramGene_ARTERIAL_to_all.pdf", width = 20, height =16)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[2]], sources.use = "ARTERIAL", targets.use = c(1:15), lab.cex = 0.5, title.name = paste0("Signaling from Perineurial cells to all", names(object.list)[i]))
}
dev.off()


pdf(file ="VEGF_ChorddiagramGene_ARTERIAL_to_TIP3.pdf", width = 20, height =16)
# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = "ARTERIAL", targets.use = "TIP_3",  title.name = paste0("Signaling from Arterial cells to TIP3", names(object.list)[i]), legend.pos.x = 10)
}
dev.off()


saveRDS(cellchat, file = "cellchat_comparisonAnalysis_mouse_condition2_vs_condition1.rds")


object.list[[1]]@net$weight 

rowSums(colSums(object.list[[1]]@net$weight))

9.607221+9.208371+12.901699+12.491153+3.567118+1.906135+1.550044+1.303640+14.377060+8.292873+1.630612+1.448423+2.983043+11.396089+16.149442

setwd("/Users/maurizio.aurora/comparison")
cellchat = readRDS("cellchat_comparisonAnalysis_mouse_condition2_vs_condition1.rds")


mat_condition2 <- object.list[[1]]@net$weight 

pdf("interactions_condition2.pdf", 15, 15)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat_condition2)) {
  mat2 <- matrix(0, nrow = nrow(mat_condition2), ncol = ncol(mat_condition2), dimnames = dimnames(mat_condition2))
  mat2[i, ] <- mat_condition2[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_condition2), title.name = rownames(mat_condition2)[i])
}
dev.off()

mat_condition1 <- object.list[[2]]@net$weight 

pdf("interactions_condition1.pdf",  15, 15)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat_condition1)) {
  mat2 <- matrix(0, nrow = nrow(mat_condition1), ncol = ncol(mat_condition1), dimnames = dimnames(mat_condition1))
  mat2[i, ] <- mat_condition1[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_condition1), title.name = rownames(mat_condition1)[i])
}
dev.off()

netAnalysis_contribution(cellchat_condition1, signaling = pathways.show)
cellchat_condition1




pdf(file ="VEGF_ChorddiagramGene_ARTERIAL_to_all_condition1.pdf", width = 10, height =8)
netVisual_chord_gene(object.list[[2]], sources.use = "ARTERIAL", targets.use = c(1:15), lab.cex = 0.5, title.name = paste0("Signaling from Perineurial cells to all condition1"))
dev.off()


plotGeneExpression(object.list[[1]], signaling = "VEGF")


for (i in 1:length(object.list)) {
  object.list[[i]] = netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP") 
}

setwd("/Users/maurizio.aurora/comparison")
cellchat_condition2 = readRDS("cellchat_condition2.Rds")
cellchat_condition1 = readRDS("cellchat_condition1.Rds")


cellchat_condition2 <- netAnalysis_computeCentrality(cellchat_condition2, slot.name = "netP") 
cellchat_condition1 <- netAnalysis_computeCentrality(cellchat_condition1, slot.name = "netP") 

object.list <- list(condition2 = cellchat_condition2, condition1 = cellchat_condition1)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

?netAnalysis_signalingChanges_scatter

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "ARTERIAL") + xlim(0, 0.4)+ylim(0, 0.4)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MES_ENDONEUR") + xlim(0, 0.2)+ylim(0, 0.3)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))



library(NMF)
#> Loading required package: pkgmaker
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)
selectK(object.list[[1]], pattern = "outgoing")

object.list[[2]]

nPatterns = 4

pdf("condition2_outgoing_cell_communic_patters.pdf", 12,10)
object.list[[1]] <- identifyCommunicationPatterns(object.list[[1]], pattern = "outgoing", k = nPatterns, width = 10,
                                                  height = 12)
dev.off()

# river plot
pdf("condition2_outgoing_alluvial.pdf")
netAnalysis_river(object.list[[1]], pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
dev.off()


# dot plot
pdf("condition2_outgoing_dotplot.pdf")
netAnalysis_dot(object.list[[1]], pattern = "outgoing")
dev.off()


selectK(object.list[[1]], pattern = "incoming")

nPatterns = 5



pdf("condition2_incoming_cell_communic_patters.pdf", 12,10)
object.list[[1]] <- identifyCommunicationPatterns(object.list[[1]], pattern = "incoming", k = nPatterns, width = 10,
                                                  height = 12)
dev.off()

# river plot
pdf("condition2_incoming_alluvial.pdf")
netAnalysis_river(object.list[[1]], pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function
dev.off()


# dot plot
pdf("condition2_incoming_dotplot.pdf", 10, 10)
netAnalysis_dot(object.list[[1]], pattern = "incoming")
dev.off()



selectK(object.list[[2]], pattern = "outgoing")

?identifyCommunicationPatterns

nPatterns = 3

pdf("condition1_outgoing_cell_communic_patters.pdf", 12,10)
object.list[[2]] <- identifyCommunicationPatterns(object.list[[2]], pattern = "outgoing", k = nPatterns, width = 10,
                                                  height = 12)
dev.off()

# river plot
pdf("condition1_outgoing_alluvial.pdf")
netAnalysis_river(object.list[[2]], pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
dev.off()


# dot plot
pdf("condition1_outgoing_dotplot.pdf",10,10)
netAnalysis_dot(object.list[[2]], pattern = "outgoing")
dev.off()

selectK(object.list[[2]], pattern = "outgoing")

?identifyCommunicationPatterns

nPatterns = 3
#################

selectK(object.list[[2]], pattern = "incoming")
nPatterns = 3



pdf("condition1_incoming_cell_communic_patters.pdf", 12,10)
object.list[[2]] <- identifyCommunicationPatterns(object.list[[2]], pattern = "incoming", k = nPatterns, width = 10,
                                                  height = 12)
dev.off()

# river plot
pdf("condition1_incoming_alluvial.pdf")
netAnalysis_river(object.list[[2]], pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function
dev.off()


# dot plot
pdf("condition1_incoming_dotplot.pdf", 10, 10)
netAnalysis_dot(object.list[[2]], pattern = "incoming")
dev.off()

object.list[[2]]@net$weight
clear
head(object.list[[1]]@LR)

object.list[[2]]@data.signaling


setwd("/Users/maurizio.aurora/")

saveRDS(cellchat, file = "cellchat_comparisonAnalysis_mouse_condition2_vs_condition1_bis.rds")
#########
saveRDS(cellchat_condition2, "cellchat_condition2_bis.Rds")
saveRDS(cellchat_condition1, "cellchat_condition1_bis.Rds")


cellchat = readRDS(file = "cellchat_comparisonAnalysis_mouse_condition2_vs_condition1_bis.rds")

cellchat_condition2 = readRDS("cellchat_condition2_bis.Rds")
cellchat_condition1 = readRDS("cellchat_condition1_bis.Rds")

df.net.cellchat_condition2_LR <- subsetCommunication(cellchat_condition2)
df.net.cellchat_condition1_LR <- subsetCommunication(cellchat_condition1)


filename_xls <- 'df.net.cellchat_condition2_LR.xlsx'

write.xlsx(df.net.cellchat_condition2_LR,
           file= filename_xls, 
           row.names = F,
           asTable = T)



filename_xls <- 'df.net.cellchat_condition1_LR.xlsx'

write.xlsx(df.net.cellchat_condition1_LR,
           file= filename_xls, 
           row.names = F,
           asTable = T)

df.net.cellchat_condition2_sigpat <- subsetCommunication(cellchat_condition2, slot.name = "netP" )
df.net.cellchat_condition1_sigpat <- subsetCommunication(cellchat_condition1, slot.name = "netP" )
df.net.cellchat_sigpat <- subsetCommunication(cellchat, slot.name = "netP" )


head(df.net.cellchat_sigpat)
?subsetCommunication


filename_xls <- 'df.net.cellchat_condition2_sigpat.xlsx'

write.xlsx(df.net.cellchat_condition2_sigpat,
           file= filename_xls, 
           row.names = F,
           asTable = T)



filename_xls <- 'df.net.cellchat_condition1_sigpat.xlsx'

write.xlsx(df.net.cellchat_condition1_sigpat,
           file= filename_xls, 
           row.names = F,
           asTable = T)


slot.name = "netP" 

head(df.net.cellchat_condition2)

head(df.net, 2)




pdf("injury_outgoing_cell_communic_patters_.pdf", 12,10)
object.list[[2]] <- identifyCommunicationPatterns(object.list[[2]], pattern = "outgoing", k = nPatterns, width = 10,
                                                  height = 12)
dev.off()

# river plot
pdf("injury_outgoing_alluvial.pdf")
netAnalysis_river(object.list[[2]], pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
dev.off()


