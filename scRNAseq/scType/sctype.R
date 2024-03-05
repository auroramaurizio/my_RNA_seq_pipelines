#CAF new assignment with sctype
#https://github.com/IanevskiAleksandr/sc-type
#https://www.nature.com/articles/s41467-022-28803-w

library(Seurat)

setwd("/Users/maurizio.aurora/paper_figures/")

seurat_object <- readRDS("seurat_object_281023.RDS")

# start building database

apCAF <- read.delim("/Users/maurizio.aurora/CAFs/CancerDiscovery_2019/apCAF_signature.txt",h=T,row.names = "Associated.Gene.Name")
iCAF <- read.delim("/Users/maurizio.aurora/CAFs/CancerDiscovery_2019/iCAF_signature.txt",h=T,row.names = "Associated.Gene.Name")
myCAF <- read.delim("/Users/maurizio.aurora/CAFs/CancerDiscovery_2019/myCAF_signature.txt",h=T,row.names = "Associated.Gene.Name")


apcafs <- row.names(apCAF)
icafs <- row.names(iCAF)
mycafs <- row.names(myCAF)

# converting vector
apcafs_q <- shQuote(apcafs, type = "cmd")
# combining elements using , 
apcafs_c <- noquote(paste(apcafs_q, collapse = ", "))

icafs_q <- shQuote(icafs, type = "cmd")
# combining elements using , 
icafs_c <- noquote(paste(icafs_q, collapse = ", "))

mycafs_q <- shQuote(mycafs, type = "cmd")
# combining elements using , 
mycafs_c <- noquote(paste(mycafs_q, collapse = ", "))


tissueType <- c("CAF","CAF","CAF")
cellName <- c('apCAF', 'iCAF', 'myCAF')
geneSymbolmore1 <- c(apcafs_c,icafs_c,mycafs_c)
geneSymbolmore2 <- ""
shortName <- c('apCAF', 'iCAF', 'myCAF')

df <- data.frame(tissueType,cellName,geneSymbolmore1,geneSymbolmore2,shortName)

write.xlsx(df,"ScTypeDB_Califano_CAF.xlsx", quotes = F)

# end building database

colors <- c('mesothelial cells' = '#7CAE00', 
            'Ly6c+  perivascular adventicial cells' = '#CD9600',  
            'Cxcl9 positive TRC cells'='#7F00FF',
            'Red Pulp Reticular Cells 1' = 'brown', 
            'Pericytes 1' ='#FF6666', 
            'T-B reticular cells' = '#00BFC4', 
            'Pericytes 2' = '#FF9999',
            'Proliferating' = '#F9D313',
            'Thy1+ TRCs' = '#9999FF',
            'FDC type 1 cells' = '#FF3399')


DimPlot(seurat_object, cols = colors, pt.size = 1)


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

setwd("/Users/maurizio.aurora/Downloads")

#use sc-type to define clusters according to gene expression
library(HGNChelper)
load("TRCissue.RData")
DimPlot(seurat_object, group.by="customclassif")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
db_CAF <- read.xlsx("ScTypeDB_Califano_CAF.xlsx")

setwd("/Users/maurizio.aurora/Downloads")
# DB file
db_ = "ScTypeDB_Califano_CAF.xlsx"

tissue = "CAF"
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = seurat_object[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

seurat_object@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_object@meta.data$customclassif[seurat_object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

pdf("CAFs_Califano_scType.pdf")
DimPlot(seurat_object, group.by = "customclassif", pt.size = 1, cols = c('myCAF' = '#fbb040', 'iCAF' = '#9e1f63', 'apCAF'='#662d91', "Unknown" = "grey"))
dev.off()
