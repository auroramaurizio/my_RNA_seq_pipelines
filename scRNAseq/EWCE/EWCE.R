#############
#EWCE
#############

#install.packages("devtools")
#library(devtools)
#install_github("neurogenomics/EWCE")

library(devtools)
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
theme_set(theme_cowplot())

#the main function we require is 
#help(generate.celltype.data)


#it takes as input 
#exp
#annotLevels

#generate exp from my scdataset
#Numerical matrix with row for each gene and column for each cell. 
#Row names are MGI/HGNC gene symbols. 
#Column names are cell IDs which can be cross referenced against the annot data frame.
integrated_3PIP <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3PIP/integrated_EC_3PIP.Rds")
integrated = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/Subclustering_mm10_M16_TdTomato_cc_2/seurat_objects/all_cell_types_scBonanomi_full.Rds")

#extract cell names for each ident in the EC dataset integrated_3PIP
six <-WhichCells(object=integrated_3PIP, idents="six")
one <-WhichCells(object=integrated_3PIP, idents="one")
two <-WhichCells(object=integrated_3PIP, idents="two")
three <-WhichCells(object=integrated_3PIP, idents="three")
five <-WhichCells(object=integrated_3PIP, idents="five")
four_M <-WhichCells(object=integrated_3PIP, idents="four-")
four_P <-WhichCells(object=integrated_3PIP, idents="four+")
seven_M <-WhichCells(object=integrated_3PIP, idents="seven-")
seven_P <-WhichCells(object=integrated_3PIP, idents="seven+")

#highligh specific cells in the full dataset integrated
DimPlot(integrated, label=T, cells.highlight=one, cols.highlight = c("cyan"), cols= "grey")
DimPlot(integrated)

# since the full dataset object lacks three cells labels. Assign them.
Idents(object = integrated, cells = three) <- 'three'

# check if everything went well
DimPlot(integrated)

DefaultAssay(integrated) = "RNA"

#https://www.cellphonedb.org/faq-and-troubleshooting

# take raw data. In seurat 3 raw.data are called counts
count_raw =as.matrix(GetAssayData(object = integrated, slot = "counts"))

# think about performing some kind of normalization

# gene names are stored as rownames
#rownames(count_raw)
# cell names are stored as colnames
#colnames(count_raw)

#A different number of reads is found across each cell. 
#We suggest using scTransform to normalise for differences due to cell size, 
#then linearly scale. Note that this might be slow.

# devtools::install_github(repo = 'ChristophH/sctransform')

#Note that this is optional (and was not used in the original EWCE publication) 
#so by all means ignore this scTransform step.


#library(sctransform)
#scT = sctransform::vst(count_raw, return_cell_attr = TRUE)

#count_raw$exp_scT = correct_counts(scT, count_raw) # umi_corrected

#count_raw$exp_scT_normed = Matrix::t(Matrix::t(count_raw$exp_scT)*(1/Matrix::colSums(count_raw$exp_scT)))

#generate annotLevels
#List with arrays of strings containing the cell type names associated with each column in exp

# create meta data object
meta_data <- cbind(rownames(integrated@meta.data),as.data.frame(Idents(integrated)))
# create meta data "level2class" definition. It is more precise with EC subset definition
colnames(meta_data)<- c("cell_id","level2class")
# create "level1class". all EC subtypes will be called EC
meta_data$level1class <- gsub("one", "EC", meta_data$level2class)
meta_data$level1class <- gsub("two", "EC", meta_data$level1class)
meta_data$level1class <- gsub("three", "EC", meta_data$level1class)
meta_data$level1class <- gsub("four\\-", "EC", meta_data$level1class)
meta_data$level1class <- gsub("four\\+", "EC", meta_data$level1class)
meta_data$level1class <- gsub("five", "EC", meta_data$level1class)
meta_data$level1class <- gsub("six", "EC", meta_data$level1class)
meta_data$level1class <- gsub("seven\\-", "EC", as.character(meta_data$level1class))
meta_data$level1class <- gsub("seven\\+", "EC", meta_data$level1class)
#check if everything went well
unique(meta_data$level1class)


level1class = meta_data$level1class

level2class = meta_data$level2class
  
head(level2class)

#cell_id                        level1class                          level2class


#exp_merged_DROPPED = drop.uninformative.genes(exp=merged_KI$exp, level2annot = merged_KI$annot$level2class)




gene="Lyve1"
annot = meta_data
exp = count_raw
cellExpDist = data.frame(e=exp[gene,],l1=annot[colnames(exp),]$level1class)
ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(gene)


#1) Drop genes which do not show significant evidence of varying between level 
exp_DROPPED = drop.uninformative.genes(exp=count_raw,level2annot = level2class)
#2) Calculate cell type averages and specificity for each gene 
annotLevels = list(level1class=level1class,level2class=level2class)
fNames = generate.celltype.data(exp=exp_DROPPED,annotLevels=annotLevels,groupName="DB")
print(fNames)
#help(generate.celltype.data)
#3) Drop all genes which do not have 1:1 mouse:human orthologs
fNames = filter.genes.without.1to1.homolog(fNames)
print(fNames)

#load(fNames[2])

# use the List generated using generate.celltype.data and a differential expression 
# table and determines the probability of cell-type enrichment in the up & down regulated genes
help(ewce_expression_data)

#load dge tables
crush_D1_vs_crush_D2 = read.xlsx(
  "/DGE_results.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

crush_D1_vs_crush_D2$MGI.symbol = rownames(crush_D1_vs_crush_D2)
  
crush_D2_vs_intact = read.xlsx(
    "/DGE_results.xlsx",
    sheet = 2,
    startRow = 1,
    colNames = TRUE,
    rowNames = TRUE,
)

crush_D2_vs_intact$MGI.symbol = rownames(crush_D2_vs_intact)

crush_D1_vs_intact = read.xlsx(
    "/DGE_results.xlsx",
    sheet = 3,
    startRow = 1,
    colNames = TRUE,
    rowNames = TRUE,
)

KO_vs_WT = read.xlsx(
  "DGE_results_C.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

KO_vs_WT$MGI.symbol = rownames(KO_vs_WT)

head(KO_vs_WT)

load(file = "~//CellTypeData_DB.rda") #now everything is in ctd

help(ewce_expression_data)

tt_results_KO_vs_WT = ewce_expression_data(sct_data=ctd,tt=KO_vs_WT,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange", 
                                                     reps = 10000)

tt_results_crush_D1_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange", 
                                                     reps = 10000)

tt_results_crush_D2_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D2_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                      reps = 10000)

tt_results_crush_D1_vs_crush_D2 = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_crush_D2,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                        reps = 10000)

head(tt_results_KO_vs_WT$bootstrap_data.down)

rownames(tt_results_KO_vs_WT$joint_results)

ewce.plot(tt_results_KO_vs_WT$joint_results)$plain

ewce.plot(tt_results_crush_D1_vs_crush_D2$joint_results)$plain

ewce.plot(tt_results_crush_D2_vs_intact$joint_results)$plain

ewce.plot(tt_results_crush_D1_vs_intact$joint_results)$plain



tt_results_2_KO_vs_WT = ewce_expression_data(sct_data=ctd,tt=KO_vs_WT,annotLevel=2,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange", 
                                                       reps = 10000)

tt_results_2_crush_D1_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_intact,annotLevel=2,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange", 
                                                     reps = 10000)

tt_results_2_crush_D2_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D2_vs_intact,annotLevel=2,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                      reps = 10000)

tt_results_2_crush_D1_vs_crush_D2 = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_crush_D2,annotLevel=2,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                        reps = 10000)

tt_results_2_KO_vs_WT

ewce.plot(tt_results_2_KO_vs_WT$joint_results)$plain

ewce.plot(tt_results_2_crush_D1_vs_crush_D2$joint_results)$plain

ewce.plot(tt_results_crush_D1_vs_crush_D2$joint_results)$plain


ewce.plot(tt_results_2_crush_D2_vs_intact$joint_results)$plain

ewce.plot(tt_results_crush_D2_vs_intact$joint_results)$plain

ewce.plot(tt_results_2_crush_D1_vs_intact$joint_results)$plain

ewce.plot(tt_results_crush_D1_vs_intact$joint_results)$plain



#Generating bootstrap plots for transcriptomes
#A common request is to explain which differentially expressed genes are associated with a cell type…

help(generate.bootstrap.plots.for.transcriptome)

#generate.bootstrap.plots.for.transcriptome takes a genelist and a single 
#cell type transcriptome dataset and generates plots which show how the expression 
#of the genes in the list compares to those in randomly generated gene lists

full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=crush_D1_vs_intact,annotLevel=2,
                                                          full_results=tt_results_2_crush_D1_vs_intact,
                                                          listFileName="crush_D1_vs_intact",reps=10,ttSpecies="mouse",
                                                          sctSpecies="mouse", sortBy="log2FoldChange")


full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=crush_D2_vs_intact,annotLevel=2,
                                                          full_results=tt_results_2_crush_D2_vs_intact,
                                                          listFileName="crush_D2_vs_intact",reps=10,ttSpecies="mouse",
                                                          sctSpecies="mouse", sortBy="log2FoldChange")


full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=crush_D1_vs_crush_D2,annotLevel=2,
                                                          full_results=tt_results_2_crush_D1_vs_crush_D2,
                                                          listFileName="crush_D1_vs_crush_D2",reps=10,ttSpecies="mouse",
                                                          sctSpecies="mouse", sortBy="log2FoldChange")


full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=KO_vs_WT,annotLevel=2,
                                                          full_results=tt_results_2_KO_vs_WT,
                                                          listFileName="KO_vs_WT",reps=10,ttSpecies="mouse",
                                                          sctSpecies="mouse", sortBy="log2FoldChange")

help(generate.bootstrap.plots.for.transcriptome)

DefaultAssay(integrated) = "RNA"
DefaultAssay(integrated_3PIP) = "RNA"


########################################################
# redo for sample 1 only
########################################################



#install.packages("devtools")
#library(devtools)
#install_github("neurogenomics/EWCE")

library(devtools)
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
theme_set(theme_cowplot())

#the main function we require is 
#help(generate.celltype.data)


#it takes as input 
#exp
#annotLevels

#generate exp from my scdataset
#Numerical matrix with row for each gene and column for each cell. 
#Row names are MGI/HGNC gene symbols. 
#Column names are cell IDs which can be cross referenced against the annot data frame.
integrated_3PIP_all <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3PIP/integrated_EC_3PIP.Rds")
integrated_all = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/Subclustering_mm10_M16_TdTomato_cc_2/seurat_objects/all_cell_types_scBonanomi_full.Rds")


integrated_3PIP <- subset(x = integrated_3PIP_all, subset = stim == "1")

integrated <- subset(x = integrated_all, subset = stim == "1")

#extract cell names for each ident in the EC dataset integrated_3PIP
six <-WhichCells(object=integrated_3PIP, idents="six")
one <-WhichCells(object=integrated_3PIP, idents="one")
two <-WhichCells(object=integrated_3PIP, idents="two")
three <-WhichCells(object=integrated_3PIP, idents="three")
five <-WhichCells(object=integrated_3PIP, idents="five")
four_M <-WhichCells(object=integrated_3PIP, idents="four-")
four_P <-WhichCells(object=integrated_3PIP, idents="four+")
seven_M <-WhichCells(object=integrated_3PIP, idents="seven-")
seven_P <-WhichCells(object=integrated_3PIP, idents="seven+")

#highligh specific cells in the full dataset integrated
DimPlot(integrated, label=T, cells.highlight=one, cols.highlight = c("cyan"), cols= "grey")

DimPlot(integrated)

# since the full dataset object lacks three cells labels. Assign them.
Idents(object = integrated, cells = three) <- 'three'

# check if everything went well
DimPlot(integrated)

DefaultAssay(integrated) = "RNA"

#https://www.cellphonedb.org/faq-and-troubleshooting

# take raw data. In seurat 3 raw.data are called counts
count_raw =as.matrix(GetAssayData(object = integrated, slot = "counts"))




# think about performing some kind of normalization

# gene names are stored as rownames
#rownames(count_raw)
# cell names are stored as colnames
#colnames(count_raw)

#A different number of reads is found across each cell. 
#We suggest using scTransform to normalise for differences due to cell size, 
#then linearly scale. Note that this might be slow.

# devtools::install_github(repo = 'ChristophH/sctransform')

#Note that this is optional (and was not used in the original EWCE publication) 
#so by all means ignore this scTransform step.


library(sctransform)
scT = sctransform::vst(count_raw, return_cell_attr = TRUE)

count_raw$exp_scT = correct_counts(scT, count_raw) # umi_corrected

count_raw$exp_scT_normed = Matrix::t(Matrix::t(count_raw$exp_scT)*(1/Matrix::colSums(count_raw$exp_scT)))

#generate annotLevels
#List with arrays of strings containing the cell type names associated with each column in exp

# create meta data object
meta_data <- cbind(rownames(integrated@meta.data),as.data.frame(Idents(integrated)))
# create meta data "level2class" definition. It is more precise with EC subset definition
colnames(meta_data)<- c("cell_id","level2class")
# create "level1class". all EC subtypes will be called EC
meta_data$level1class <- gsub("one", "EC", meta_data$level2class)
meta_data$level1class <- gsub("two", "EC", meta_data$level1class)
meta_data$level1class <- gsub("three", "EC", meta_data$level1class)
meta_data$level1class <- gsub("four\\-", "EC", meta_data$level1class)
meta_data$level1class <- gsub("four\\+", "EC", meta_data$level1class)
meta_data$level1class <- gsub("five", "EC", meta_data$level1class)
meta_data$level1class <- gsub("six", "EC", meta_data$level1class)
meta_data$level1class <- gsub("seven\\-", "EC", as.character(meta_data$level1class))
meta_data$level1class <- gsub("seven\\+", "EC", meta_data$level1class)
#check if everything went well
unique(meta_data$level1class)


level1class = meta_data$level1class

level2class = meta_data$level2class

head(level2class)

#cell_id                        level1class                          level2class


#exp_merged_DROPPED = drop.uninformative.genes(exp=merged_KI$exp, level2annot = merged_KI$annot$level2class)




gene="Lyve1"
annot = meta_data
exp = count_raw
cellExpDist = data.frame(e=exp[gene,],l1=annot[colnames(exp),]$level1class)
ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(gene)



#1) Drop genes which do not show significant evidence of varying between level 
exp_DROPPED = drop.uninformative.genes(exp=count_raw,level2annot = level2class)
#2) Calculate cell type averages and specificity for each gene 
annotLevels = list(level1class=level1class,level2class=level2class)
fNames = generate.celltype.data(exp=exp_DROPPED,annotLevels=annotLevels,groupName="DB")
print(fNames)
#help(generate.celltype.data)
#3) Drop all genes which do not have 1:1 mouse:human orthologs
fNames = filter.genes.without.1to1.homolog(fNames)
print(fNames)

load(fNames[2])

# use the List generated using generate.celltype.data and a differential expression 
# table and determines the probability of cell-type enrichment in the up & down regulated genes
help(ewce_expression_data)

#load dge tables
crush_D1_vs_crush_D2 = read.xlsx(
  "/DGE_results.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

crush_D1_vs_crush_D2$MGI.symbol = rownames(crush_D1_vs_crush_D2)

crush_D2_vs_intact = read.xlsx(
  "/DGE_results.xlsx",
  sheet = 2,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

crush_D2_vs_intact$MGI.symbol = rownames(crush_D2_vs_intact)

crush_D1_vs_intact = read.xlsx(
  "/DGE_results.xlsx",
  sheet = 3,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

crush_D1_vs_intact$MGI.symbol = rownames(crush_D1_vs_intact)

head(crush_D1_vs_intact)

load(file = "~//CellTypeData_DB.rda") #now everything is in ctd

help(ewce_expression_data)
tt_results_crush_D1_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange", 
                                                     reps = 10000)

tt_results_crush_D2_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D2_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                      reps = 10000)

tt_results_crush_D1_vs_crush_D2 = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_crush_D2,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                        reps = 10000)

ewce.plot(tt_results_crush_D1_vs_crush_D2$joint_results)$plain

ewce.plot(tt_results_crush_D2_vs_intact$joint_results)$plain

ewce.plot(tt_results_crush_D1_vs_intact$joint_results)$plain



tt_results_2_crush_D1_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_intact,annotLevel=2,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange", 
                                                       reps = 10000)

tt_results_2_crush_D2_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D2_vs_intact,annotLevel=2,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                        reps = 10000)

tt_results_2_crush_D1_vs_crush_D2 = ewce_expression_data(sct_data=ctd,tt=crush_D1_vs_crush_D2,annotLevel=2,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                          reps = 10000)

ewce.plot(tt_results_2_crush_D1_vs_crush_D2$joint_results)$plain

ewce.plot(tt_results_2_crush_D2_vs_intact$joint_results)$plain

ewce.plot(tt_results_2_crush_D1_vs_intact$joint_results)$plain


#Generating bootstrap plots for transcriptomes
#A common request is to explain which differentially expressed genes are associated with a cell type…

help(generate.bootstrap.plots.for.transcriptome)

#generate.bootstrap.plots.for.transcriptome takes a genelist and a single 
#cell type transcriptome dataset and generates plots which show how the expression 
#of the genes in the list compares to those in randomly generated gene lists

full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=crush_D1_vs_intact,annotLevel=2,
                                                          full_results=tt_results_2_crush_D1_vs_intact,
                                                          listFileName="crush_D1_vs_intact",reps=10,ttSpecies="mouse",
                                                          sctSpecies="mouse", sortBy="log2FoldChange")


full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=crush_D2_vs_intact,annotLevel=2,
                                                          full_results=tt_results_2_crush_D2_vs_intact,
                                                          listFileName="crush_D2_vs_intact",reps=10,ttSpecies="mouse",
                                                          sctSpecies="mouse", sortBy="log2FoldChange")


full_results = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=crush_D1_vs_crush_D2,annotLevel=2,
                                                          full_results=tt_results_2_crush_D1_vs_crush_D2,
                                                          listFileName="crush_D1_vs_crush_D2",reps=10,ttSpecies="mouse",
                                                          sctSpecies="mouse", sortBy="log2FoldChange")




help(generate.bootstrap.plots.for.transcriptome)

DefaultAssay(integrated) = "RNA"
DefaultAssay(integrated_3PIP) = "RNA"



#HGNC.symbol      ID_REF       logFC  AveExpr         t      P.Value adj.P.Val
###################################################################################
#examples from the vignette https://nathanskene.github.io/EWCE/articles/EWCE.html#application-to-transcriptomic-data-1
###################################################################################


data(cortex_mrna)
gene="Necab1"
cellExpDist = data.frame(e=cortex_mrna$exp[gene,],l1=cortex_mrna$annot[colnames(cortex_mrna$exp),]$level1class)
ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


cortex_mrna
rownames(cortex_mrna$exp[gene,])


download.file("goo.gl/r5Y24y",
              destfile="expression_mRNA_17-Aug-2014.txt") 

path = "expression_mRNA_17-Aug-2014.txt"

cortex_mrna  = load.linnarsson.sct.data(path)

data(cortex_mrna)

gene="Necab1"
cellExpDist = data.frame(e=cortex_mrna$exp[gene,],l1=cortex_mrna$annot[colnames(cortex_mrna$exp),]$level1class)

head(cellExpDist)

e=cortex_mrna$exp[gene,]

l1=cortex_mrna$annot[colnames(cortex_mrna$exp),]$level1class


#Important note: you do NOT have to format your input single cell data like the Linnarsson data.
#Just read it into R such that you have an expression matrix and an annotation data frame. 
#The three columns that you must have in the annotation data frame are “cell_id”, 
#“level1class” and “level2class”.
#i need a df with the genenames as rownames and the cellnames as colnames
# e.g. 
rownames(cortex_mrna$exp)
colnames(cortex_mrna$exp)

#I need another df with cell names as colnames cell_id  level1class level2class
head(cortex_mrna$annot)



ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


if(!file.exists("MRK_List2.rpt")){
  download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", destfile="MRK_List2.rpt")
}
cortex_mrna$exp = fix.bad.mgi.symbols(cortex_mrna$exp,mrk_file_path="MRK_List2.rpt")

head(cortex_mrna$exp)

nKeep = 1000
must_keep = c("Apoe","Gfap","Gapdh")
set.seed(123458)
keep_genes = c(must_keep,sample(rownames(cortex_mrna$exp),997))
cortex_mrna$exp = cortex_mrna$exp[keep_genes,]



# devtools::install_github(repo = 'ChristophH/sctransform')
library(sctransform)
scT = sctransform::vst(cortex_mrna$exp, return_cell_attr = TRUE)

cortex_mrna$exp_scT = correct_counts(scT, cortex_mrna$exp) # umi_corrected


cortex_mrna$exp_scT_normed = Matrix::t(Matrix::t(cortex_mrna$exp_scT)*(1/Matrix::colSums(cortex_mrna$exp_scT)))


# Generate celltype data for just the cortex/hippocampus data
exp_CortexOnly_DROPPED = drop.uninformative.genes(exp=cortex_mrna$exp_scT_normed,level2annot = cortex_mrna$annot$level2class)
annotLevels = list(level1class=cortex_mrna$annot$level1class,level2class=cortex_mrna$annot$level2class)
fNames_CortexOnly = generate.celltype.data(exp=exp_CortexOnly_DROPPED,annotLevels=annotLevels,groupName="kiCortexOnly")
print(fNames_CortexOnly)

fNames_CortexOnly = filter.genes.without.1to1.homolog(fNames_CortexOnly)
print(fNames_CortexOnly)

load(fNames_CortexOnly[1])

# Download the hypothalamus data and unzip
if(!file.exists("GSE74672_expressed_mols_with_classes.xlsx")){
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74672/suppl/GSE74672_expressed_mols_with_classes.xlsx.gz", destfile="GSE74672_expressed_mols_with_classes.xlsx.gz")
  system("gunzip GSE74672_expressed_mols_with_classes.xlsx.gz")
}

# Read in the hypothalamus data
hypo_dat = read_excel("GSE74672_expressed_mols_with_classes.xlsx")

head(hypo_dat)

# Extract the expression data, gene symbols and annotation data
exp = data.matrix(hypo_dat[12:dim(hypo_dat)[1],2:dim(hypo_dat)[2]])
#colnames(exp)
rownames(exp) = data.frame(hypo_dat[12:dim(hypo_dat)[1],1])[,1]

level1class = data.frame(level1class=t(hypo_dat[1,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
level2class = data.frame(leve2class=t(hypo_dat[2,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
cell_id     = colnames(hypo_dat)[2:dim(hypo_dat)[2]]
hypo_annot  = data.frame(cell_id=cell_id,level1class=level1class,level2class=level2class,stringsAsFactors = FALSE)

#annotation file
head(hypo_annot)

# Drop the glia and unclassified cells (which don't have level 2  annotations)
hypo_annot  = hypo_annot[!is.na(hypo_annot$level2class) & !hypo_annot$level2class=="uc",]
hypo_exp    = exp[,hypo_annot$cell_id]

# Make the celltype names more aesthetically pleasing
hypo_annot$level2class=gsub(",",";",hypo_annot$level2class)
hypo_annot$level1class[grep("Oxt;|^Avp",hypo_annot$level2class)] = "Oxytocin / Vasopressin Expressing Neurons"
hypo_annot$level1class[grep("^Th;|^Dopamine",hypo_annot$level2class)] = "Hypothalamic Dopaminergic Neurons"
hypo_annot$level1class[grepl("^Vglut2|^Trh|^Qrfp|^Hcrt|^Pmch|^Adcyap1|^Npvf|^Ghrh|^Hmit|^Nms|^Vip;|^Per2|Tnr$|^Gad-low;Gnrh",hypo_annot$level2class) & grepl("neurons",hypo_annot$level1class)] = "Hypothalamic Glutamatergic Neurons"
hypo_annot$level1class[grepl("GABA|^Sst|^Crh|^Npy|^Pomc|^Galanin|^Otof|Pnoc$|^Calcr-high",hypo_annot$level2class) & grepl("^neurons$",hypo_annot$level1class)] = "Hypothalamic GABAergic Neurons"
hypo_annot$level2class[hypo_annot$level2class!=""] = sprintf("Hypothalamic %s Neuron",hypo_annot$level2class[hypo_annot$level2class!=""])

# Fix bad MGI symbols
hypo_exp_CORRECTED = fix.bad.mgi.symbols(hypo_exp)


# Merge the datasets
merged_KI = merge_two_expfiles(exp1=hypo_exp_CORRECTED,  exp2=cortex_mrna$exp,
                               annot1=hypo_annot,        annot2=cortex_mrna$annot,
                               name1="Hypothalamus (KI)", name2="Cortex/Hippo (KI)")

# Drop genes which don't vary significantly between cell types
exp_merged_DROPPED = drop.uninformative.genes(exp=merged_KI$exp, level2annot = merged_KI$annot$level2class)

# Calculate specificity data
annotLevels = list(level1class=merged_KI$annot$level1class,level2class=merged_KI$annot$level2class)
fNames_MergedKI = generate.celltype.data(exp=exp_merged_DROPPED,annotLevels,"MergedKI")
fNames_MergedKI = filter.genes.without.1to1.homolog(fNames_MergedKI)
load(fNames_MergedKI[2])
fNames_MergedKI

help(generate.celltype.data)



rownames(exp_merged_DROPPED)
colnames(exp_merged_DROPPED)

rownames(annotLevels)
colnames(annotLevels)
######################

data(tt_alzh)
tt_results = ewce_expression_data(sct_data=ctd,tt=tt_alzh,annotLevel=1,ttSpecies="human",sctSpecies="mouse")

ewce.plot(tt_results$joint_results)$plain


head(tt_alzh_BA36)
head(tt_alzh_BA44)

head(ctd[])


head(tt_alzh)

data(celltype_data)

cortex_mrna




integrated_rn = readRDS("ntegration_assigned.Rds")

DimPlot(integrated_rn)

DefaultAssay(integrated_rn) = "RNA"

#https://www.cellphonedb.org/faq-and-troubleshooting

#DimPlot()
# take raw data and normalise it. In seurat 3 raw.data are called counts
count_raw =as.matrix(GetAssayData(object = integrated_rn, slot = "counts"))

cortex_mrna
annot
$exp


ctd


tt_results = ewce_expression_data(sct_data=cortex_mrna,tt=tt_alzh,annotLevel=1,ttSpecies="human",sctSpecies="mouse")

head(cortex_mrna)

help(ewce_expression_data)

tail(cortex_mrna$exp)

hypo_dat

############################

level1class
level2class
cell_id
hypo_annot

exp

###########################

data("ctd")
set.seed(1234)
library(reshape2)
genes = c("Apoe","Gfap","Gapdh")
exp = melt(cbind(ctd[[1]]$mean_exp[genes,],genes),id.vars="genes")
colnames(exp) = c("Gene","Cell","AvgExp")
ggplot(exp)+geom_bar(aes(x=Cell,y=AvgExp),stat="identity")+facet_grid(Gene~.)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

data(tt_alzh)


reps=100 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
level=1 
tt_results = ewce_expression_data(sct_data=ctd,tt=tt_alzh,annotLevel=level,ttSpecies="human",sctSpecies="mouse",reps=reps)

ewce.plot(tt_results$joint_results)$plain

# tt_alzh
# HGNC.symbol      ID_REF       logFC   AveExpr         t      P.Value adj.P.Val

# ctd