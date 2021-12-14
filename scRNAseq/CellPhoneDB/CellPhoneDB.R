
suppressMessages(library(ggbeeswarm))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
library(scales)
library(viridis)
library(SeuratWrappers)
library(slingshot)
require(BiocStyle)
library(SCENT)
library(SingleCellExperiment)
#print(version) R version 3.6.1
setwd("/Users/maurizio.aurora/Documents/BonanomiD_1287_scRNA_injury_TdT")
packageVersion('Seurat') #3.2.2
getwd()

############
############
# intact
# integrate together stst from Carr
############
############


integrated_rn = readRDS("integrated_rn.rds")

DimPlot(integrated_rn, split.by = "stim")

sample <- subset(x = integrated_rn, subset = stim == c("1", "3"))

DimPlot(sample, split.by = "stim")

DefaultAssay(sample) = "RNA"

#https://www.cellphonedb.org/faq-and-troubleshooting

#DimPlot()
# take raw data and normalise it. In seurat 3 raw.data are called counts
count_raw =as.matrix(GetAssayData(objt = sample, slot = "counts"))

# gene names are stored as rownames
#colnames(count_raw)

# convert mouse gene names in human gene names

require("biomaRt")

convertMouseGeneList <- function(x){
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

converted=convertMouseGeneList(rownames(count_raw))
#converted=unique(converted)

# if same MGI.symbol is associated to multiple human genes keep only the first one
occurrence <- converted[!duplicated(converted$MGI.symbol),]

#occurrence[grep("Alg1", occurrence$MGI.symbol), ]

# remove pseudogenes
#occurrence=occurrence[!grepl("^Gm", occurrence$MGI.symbol),]

# count occurrence of a test gene
#occurrence[grep("Ercc5", occurrence$HGNC.symbol), ]

rownames(occurrence) = occurrence$MGI.symbol

merged=merge(occurrence,count_raw,by="row.names")

head(merged)
merged <- merged[-c(1,2)]


#sum counts of duplicated rows
agg=aggregate(merged[,sapply(merged,is.numeric)],merged["HGNC.symbol"],sum)

rownames(agg) = agg$HGNC.symbol

agg$HGNC.symbol <- NULL

#normalize
count_norm <- apply(agg, 2, function(x) (x/sum(x))*10000)
Gene <- rownames(count_norm)
count_norm <- cbind(Gene,count_norm)

#head(count_norm[1:3,1:3])

#write counts file
write.table(count_norm, 'cellphonedb_count_1_3_MES.txt', sep='\t', quote=F, row.names = F)

# generating meta file
meta_data <- cbind(rownames(sample@meta.data),as.data.frame(Idents(sample)))
colnames(meta_data)<- c("Cell","cell_type")
keeps <- colnames(count_norm)
mySubset <- meta_data[meta_data$Cell %in% keeps, ]

#####  cluster is the user’s spific cluster column
write.table(mySubset, 'cellphonedb_meta_1_3_MES.txt', sep='\t', quote=F, row.names=F)



# cellphonedb command
# python -m venv cpdb-venv
# source cpdb-venv/bin/activate
# pip install cellphonedb
# running the statistical method
# as reported here https://github.com/Teichlab/cellphonedb/issues/266
# pip install --force-reinstall numpy==1.19.5
# mkdir CpDB_S2_D1_Ydens
# cellphonedb method statistical_analysis cellphonedb_meta_1_3_MES.txt cellphonedb_count_1_3_MES.txt --counts-data=gene_name --threads=2 --iterations=100

# cellphonedb plot dot_plot
# cellphonedb plot heatmap_plot cellphonedb_meta_1__3_MES.txt



toplot_matrix = read.table("/Users/maurizio.aurora/Documents/CDB_MES/count_network.txt", header=T, sep = "\t")
toplot_matrix1 <- toplot_matrix[order(toplot_matrix$SOURCE, toplot_matrix$TARGET),]

library(pheatmap)
library(reshape2)
matrix_to_plot= acast(toplot_matrix1, SOURCE~TARGET, value.var="count")

Breaks <- seq(min(0), max(130), length = 100)

# create heatmap using pheatmap
pheatmap(matrix_to_plot,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         breaks= Breaks,
         filename = 'Heatmap_1_3_MESs.pdf',
         main = '_vs_MES',
         legend = TRUE)

