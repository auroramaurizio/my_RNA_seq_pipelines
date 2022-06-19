library(CellChat)
options(stringsAsFactors = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = './'
)
knitr::opts_chunk$set(eval = FALSE)


CellChatDB <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo
write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
write.csv(complex_input, file = "complex_input_CellChatDB.csv")
write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")



setwd("/Users/maurizio.aurora/comparison")
#######
options(stringsAsFactors = FALSE)
interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")

head(cellchat@DB$geneInfo)

test = cellchat@DB$interaction
CellChatDB$interaction [grep("SEMA4A", rownames(CellChatDB$interaction )), ]
test[grep("PLXND1", rownames(test)), ]

#########


setwd("/Users/maurizio.aurora/miniconda3/lib/R/library/CellChat") # This is the folder of CellChat package downloaded from Github
CellChatDB.mouse <- CellChatDB
usethis::use_data(CellChatDB.mouse, overwrite = TRUE)

interact = CellChatDB.mouse$interaction

interact[grep("Plxnd1", interact$interaction_name_2), ]

library("Seurat")

Injury = readRDS("/Users/maurizio.aurora/Injury_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")
Intact = readRDS("/Users/maurizio.aurora/Intact_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")

DimPlot(Injury)
DimPlot(Intact)

DefaultAssay(Injury) = "RNA"

VlnPlot(Injury, "Plxnd1")

#DimPlot()
# take raw data and normalise it. In seurat 3 raw.data are called counts
#count_raw =as.matrix(GetAssayData(object = object, slot = "counts"))

#?normalizeData
#i can obtain the normalized data in this way with the normalize data function
#normalizeData(count_raw, scale.factor = 10000, do.log = TRUE)[1:5,1:5]

#or taking advantage of seurat's normalization data 
count_norm = Injury@assays$RNA@data # normalized data matrix
meta_data <- cbind(rownames(Injury@meta.data),as.data.frame(Idents(Injury)))

head(meta_data)
colnames(meta_data)<- c("Cell","labels")
head(meta_data)

# Here we load a scRNA-seq data matrix and its associated cell meta data
#load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseases
#load("/Users/suoqinjin/Documents/CellChat/tutorial/data_humanSkin_CellChat.rda")
#data.input = data_humanSkin$data # normalized data matrix
#meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
#cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
# Prepare input data for CelChat analysis
#data.input = data.input[, cell.use]
#head(data.input)
#meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
#unique(meta_data$labels) # check the cell labels
#head(meta)
#>  [1] Inflam. FIB  FBN1+ FIB    APOE+ FIB    COL11A1+ FIB cDC2        
#>  [6] LC           Inflam. DC   cDC1         CD40LG+ TC   Inflam. TC  
#> [11] TC           NKT         
#> 12 Levels: APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 ... NKT

cellchat_injury <- createCellChat(object = count_norm, meta = meta_data, group.by = "labels")
#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT

cellchat_injury <- setIdent(cellchat_injury, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_injury@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_injury@idents)) # number of cells in each cell group



# set the used database in the object
cellchat_injury@DB <- CellChatDB.mouse

?subsetData
# subset the expression data of signaling genes for saving computation cost
cellchat_injury <- subsetData(cellchat_injury) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat_injury <- identifyOverExpressedGenes(cellchat_injury, thresh.p = 0.1)
cellchat_injury <- identifyOverExpressedInteractions(cellchat_injury)
# project gene expression data onto PPI network (optional)
?computeCommunProb
cellchat_injury <- projectData(cellchat_injury, PPI.mouse)

pairLR.use = data.frame(interaction_name='SEMA4A_PLXND1')
#https://rdrr.io/github/sqjin/CellChat/man/computeCommunProbPathway.html
cellchat_injury <- computeCommunProb(cellchat_injury)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_injury <- filterCommunication(cellchat_injury, min.cells = 1)
cellchat_injury <- computeCommunProbPathway(cellchat_injury, thresh = 1)
cellchat_injury <- aggregateNet(cellchat_injury)

pathways.show <- c("SEMA4") 
pdf(file ="SEMA4_Chorddiagram_testtt.pdf", width = 20, height =16)
netVisual_aggregate(cellchat_injury, signaling = pathways.show, layout = "chord")
dev.off()



df.net.cellchat_injury_LR[grep("PLXND1", df.net.cellchat_injury_LR$interaction_name), ]

df.net.cellchat_injury_LR <- subsetCommunication(cellchat_injury)

head(df.net.cellchat_injury_LR)


Intact = readRDS("/Users/maurizio.aurora/Intact_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")

DefaultAssay(Intact) = "RNA"

#DimPlot()
# take raw data and normalise it. In seurat 3 raw.data are called counts
#count_raw =as.matrix(GetAssayData(object = object, slot = "counts"))

#?normalizeData
#i can obtain the normalized data in this way with the normalize data function
#normalizeData(count_raw, scale.factor = 10000, do.log = TRUE)[1:5,1:5]

#or taking advantage of seurat's normalization data 
count_norm = Intact@assays$RNA@data # normalized data matrix
meta_data <- cbind(rownames(Intact@meta.data),as.data.frame(Idents(Intact)))

head(meta_data)
colnames(meta_data)<- c("Cell","labels")
head(meta_data)

# Here we load a scRNA-seq data matrix and its associated cell meta data
#load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseases
#load("/Users/suoqinjin/Documents/CellChat/tutorial/data_humanSkin_CellChat.rda")
#data.input = data_humanSkin$data # normalized data matrix
#meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
#cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
# Prepare input data for CelChat analysis
#data.input = data.input[, cell.use]
#head(data.input)
#meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
#unique(meta_data$labels) # check the cell labels
#head(meta)
#>  [1] Inflam. FIB  FBN1+ FIB    APOE+ FIB    COL11A1+ FIB cDC2        
#>  [6] LC           Inflam. DC   cDC1         CD40LG+ TC   Inflam. TC  
#> [11] TC           NKT         
#> 12 Levels: APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 ... NKT

cellchat_intact <- createCellChat(object = count_norm, meta = meta_data, group.by = "labels")
#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT

cellchat_intact <- setIdent(cellchat_intact, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_intact@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_intact@idents)) # number of cells in each cell grou


cellchat_intact@DB <- CellChatDB.mouse
# subset the expression data of signaling genes for saving computation cost
# subset the expression data of signaling genes for saving computation cost
cellchat_intact <- subsetData(cellchat_intact) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat_intact <- identifyOverExpressedGenes(cellchat_intact, thresh.p = 0.1)
cellchat_intact <- identifyOverExpressedInteractions(cellchat_intact)
# project gene expression data onto PPI network (optional)
cellchat_intact <- projectData(cellchat_intact, PPI.mouse)
cellchat_intact <- computeCommunProb(cellchat_intact)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_intact <- filterCommunication(cellchat_intact, min.cells = 1)
cellchat_intact <- computeCommunProbPathway(cellchat_intact, thresh = 1)
cellchat_intact <- aggregateNet(cellchat_intact)



object.list <- list(intact = cellchat_intact, injury = cellchat_injury)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
 
pathways.show <- c("SEMA3") 
pdf(file ="SEMA4a_Chorddiagram_ndb.pdf", width = 20, height =16)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], # signaling = pathways.show, 
                       pairLR.use = data.frame(interaction_name='SEMA4A_PLXND1'))
}
dev.off()



pathways.show <- c("SEMA4A") 
pdf(file ="SEMA3C_PLXND1_Chorddiagram_ndb.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat_injury, # signaling = pathways.show, 
                       pairLR.use = data.frame(interaction_name='SEMA3C_PLXND1'), slot.name = "net")
dev.off()


pdf(file ="SEMA3E_PLXND1_Chorddiagram_ndb.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat_injury, # signaling = pathways.show, 
                     pairLR.use = data.frame(interaction_name='SEMA3C_PLXND1'), slot.name = "net")
dev.off()

?netVisual_chord_gene

pdf(file ="IL2RA_IL2RB_Chorddiagram_ndb.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat_intact, # signaling = pathways.show, 
                     pairLR.use = data.frame(interaction_name='SEMA4A_PLXNB1'), slot.name = "net")
dev.off()

IL2RA_IL2RB

Sema3E, Sema3C and Sema4A 


?netVisual_chord_cell
?netVisual_chord_gene
interaction_name

pdf(file ="SEMA4a_Chorddiagram_ndb.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, signaling = c("PLXND1"),legend.pos.x = 8)
dev.off()


pathways.show <- c("SEMA4") 
pdf(file ="SEMA4_Chorddiagram_ndb.pdf", width = 20, height =16)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], signaling = pathways.show)
}
dev.off()


pathways.show <- c("PLXND1") 
pdf(file ="PLXN_Chorddiagram_ndb1.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat_intact, signaling = pathways.show)
dev.off()

pathways.show <- c("SEMA4") 
pdf(file ="SEMA4_Chorddiagram_ndb1i.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat_injury, signaling = pathways.show)
dev.off()

cellchat

test = cellchat@DB$interaction
test[grep("SEMA4A", test$interaction_name), ]
str(df.net.cellchat_injury_LR["SEMA3C_PLXND1",])


df.net.cellchat_injury_LR[grep("SEMA3C_PLXND1", df.net.cellchat_injury_LR$interaction_name),]

cellchat_intact = readRDS("cellchat_intact_bis.Rds")
cellchat_injury = readRDS("cellchat_injury_bis.Rds")

df.net.cellchat_intact_LR <- subsetCommunication(cellchat_intact)
df.net.cellchat_injury_LR <- subsetCommunication(cellchat_injury)

cellchat_injury[grep("SEMA4A_PLXND1", cellchat_injury$interaction_name)]




# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "injury"
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


