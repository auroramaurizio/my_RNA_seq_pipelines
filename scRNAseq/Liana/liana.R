## Load required packages
getwd()
library(tidyverse)
library(magrittr)
library(liana)
library(Seurat)
#install.packages("entropy")
#devtools::install_github("tanlabcode/CytoTalk")
#install.packages("openxlsx", dependencies=TRUE)
library(openxlsx)

#packageVersion("liana")
#0.1.12


#https://saezlab.github.io/liana/ 
#https://www.nature.com/articles/s41467-022-30755-0
#https://saezlab.github.io/liana/articles/liana_tutorial.html # R tutorial
#https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html # Python tutorial


## CCC Resources
#`liana` provides CCC resources obtained and formatted via [`OmnipathR`](https://github.com/saezlab/OmnipathR) which are then converted to the appropriate format to each method.
#```{r liana_resource, warning = FALSE}
show_resources()

## CCC Methods
#Each of the resources can then be run with any of the following methods:
#  ```{r liana_method, warning = FALSE}
# Resource currently included in OmniPathR (and hence `liana`) include:
show_methods()

#Note that the different algorithms (or scoring measures) used in `sca`, `natmi`,
#`connectome`, `cellphonedb`, `cytotalk`'s crosstalk scores, and `logfc` were re-implemented in LIANA.
#Yet, the original method pipelines can be called via the `call_*` functions.

## `liana` wrapper function

#To run `liana`, we will use a down-sampled toy *HUMAN* PBMCs scRNA-Seq data set, obtained from [SeuratData](https://github.com/satijalab/seurat-data).

#`liana` takes `Seurat` and `SingleCellExperiment` objects as input, containing processed counts and clustered cells.
#```{r load_data}
#liana_path <- system.file(package = "liana")
#testdata <-
#  readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
#testdata %>% dplyr::glimpse()
#```

# Load and preprocess my seurat object ("nerve intact" dataset)
test = readRDS(file ="integrated_intact_sub.rds")

a = test
b = a$SantoCellType
c = a$aurora_cellType_int

tmp = rbind(b, c)
vettore_finale = c()
for (i in 1:dim(tmp)[2]){
  if(is.na(tmp[1,i])){
    vettore_finale = c(vettore_finale, tmp[2,i])
  } else {
    vettore_finale = c(vettore_finale, tmp[1,i])
  }
  
}
names(vettore_finale) = names(b)
Idents(test) = vettore_finale

testdata <- subset(test, idents = "UNDET", invert = TRUE)

levels(testdata)
#[1] "MES_ENDONEUR"     "MES_EPINEUR"      "MES_PERINEUR"     "MES_DIFF"        
#[5] "MES_DIVIDING"     "MES_DIFF*"        "BARR_END_CAP"     "ARTERIAL"        
#[9] "CAPILLARY_PLVAP-" "TIP_1"            "CAPILLARY_PLVAP+" "VENOUS_PLVAP+"   
#[13] "TIP_2"            "VENOUS_PLVAP-"    "TIP_3"            "MACROPHAGES"     
#[17] "SCHWANN"  
new.cluster.ids.lit <- c('MES_ENDONEUR', 
                         'MES_EPINEUR',
                         'MES_PERINEUR',
                         'MES_DIFF',
                         'MES_DIVIDING',
                         'MES_DIFF*',
                         'BARR_END_CAP',
                         'ARTERIAL',
                         'CAPILLARY_PLVAP-',
                         'IMMAT',
                         'CAPILLARY_PLVAP+',
                         'VENOUS_PLVAP+',
                         'PROLIF',
                         'VENOUS_PLVAP-',
                         'TIP',
                         'MACROPHAGES',
                         'SCHWANN'
)

names(new.cluster.ids.lit) <- levels(testdata)
object_new <- RenameIdents(testdata, new.cluster.ids.lit)

testdata = object_new

#testdata <- subset(test, idents = "MES_DIFF","MES_DIFF*", "MES_DIVIDING", "MES_ENDONEUR","MES_EPINEUR","MES_PERINEUR")
#testdata <- subset(test, idents = "UNDET", invert = TRUE)

#liana reads the label associated CelltType info
DefaultAssay(testdata) = "RNA"
testdata$label = Idents(testdata)
testdata$seurat_annotations = testdata$label

#`liana_wrap` calls a number of methods and and each method is run with the provided resource(s).

#We will now call all methods that are currently available in liana.

#Here we use only the `Consensus` (Default) CCC resource, but any of the 
#aforementioned ones (available via `show_resources()`) can be added to the `resource` parameter
#```{r liana_run, message = FALSE, print = FALSE, results='hide'}
# Run liana (for mouse data use MouseConsensus)


liana_test <- liana_wrap(testdata,
                         method = c("natmi", "connectome", "logfc", "sca", "cellphonedb","cytotalk", "cytotalk"),
                         resource = c("MouseConsensus"))


#Cell identities with less than 5 cells: MES_DIFF* were removed!
#saveRDS(liana_test, "liana_test_intact_A.Rds")
saveRDS(liana_test, "liana_test_intact_B.Rds")
#liana_test_i= readRDS("liana_test_intact.Rds")


liana_test = readRDS("liana_test_intact_B.Rds")

#liana_test= readRDS("liana_test_intact.Rds")
#liana_test

# save output of liana_test 
natmi = liana_test$natmi
connectome = liana_test$connectome
logfc =liana_test$logfc
sca = liana_test$sca
cellphonedb= liana_test$cellphonedb
cytotalk = liana_test$cytotalk

write.table(natmi, "natmi_intact.txt", sep = "\t", quote = FALSE)
write.table(connectome, "connectome_intact.txt", sep = "\t", quote = FALSE)
write.table(logfc, "logfc_intact.txt", sep = "\t", quote = FALSE)
write.table(sca, "sca_intact.txt", sep = "\t", quote = FALSE)
write.table(cellphonedb, "cellphonedb_intact.txt", sep = "\t", quote = FALSE)
write.table(cytotalk, "cytotalk_intact.txt", sep = "\t", quote = FALSE)

liana_test$natmi$
# Liana returns a list of results, each element of which corresponds to a method
liana_test %>% dplyr::glimpse()
#```
#LIANA currently provides a mixture of re-implemented methods and pipelines which externally call specific LR methods.
#By default, LIANA will call *only the internal scoring function*, i.e. those that are re-implemented in LIANA.

#One can use LIANA to also run the original methods.
#For more about the original methods see [LIANA++](https://saezlab.github.io/liana/articles/liana_devel.html).


## Aggregate and Obiain Consensus Ranks
#`liana` also provides consensus ranks for the results obtained using different
#methods. By default, `liana` will provide mean, median, and aggregate* consensus
#ranks
?liana_aggregate
#```{r liana_ranks, warning=FALSE}
# We can aggregate these results into a tibble with consensus ranks
liana_test <- liana_test %>%
  liana_aggregate()
dplyr::glimpse(liana_test)

saveRDS(liana_test, "liana_test_1_aggregate_intact.Rds")

#Voila! That's it. A very brief intro to LIANA and how to obtain the scoring functions†
#for each method implemented in it, as well as an aggregate_rank* which serves as a consensus across methods.

#(†) Note that here we focus on the scores recommended to be used to prioritize interaction in a 
#single sample system. Most of these, with the exception of SingleCellSignalR's LRscore, take
#the specificity of the cluster pair into account.

#(*) The aggregate consensus rank (`aggregate_rank`) is obtained using a re-implementation of the [`RRA`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278763/) method from the `RobustRankAggreg` package.

#`RRA` scores can be interpreted as p-values and interactions which are
#ranked consistently higher than random are assigned low scores/p-values.

## Simple DotPlot

#We will now plot the results. By default, we use the `LRscore` from SingleCellSignalR
#to represent the magnitude of expression of the ligand and receptor, and NATMI's `specificity weights`
#to show how specific a given interaction is to the `source`(L) and `target`(R) cell types.

#Note that the top 20 interactions (`ntop`), are defined as unique ligand-receptors 
#ordered consequentially, regardless of the cell types. Done subsequent to filtering to the `source_groups` and `target_groups` of interest. In this case, we plot interactions in which B cells are the `source` (express the ~ligands) and the other 3 cell types are the `target` cell types (those which express the ~receptors).

#By default, `liana_aggregate` would order LIANA's 
#output according to `rank_aggregate` in ascending order. 

#Note that using `liana_dotplot` /w `ntop` assumes that liana's results are 
#ranked accordingly! Refer to `rank_method`, if you wish to easily rank
#interactions for a single method (see example below).

#```{r liana_dotplot, warning=FALSE, fig.width=11, fig.height=8}
#liana_test = readRDS("liana_test_1_aggregate_intact.Rds")

# rename cell IDs
liana_test <- liana_test %>% 
  mutate(source = recode(source, TIP_3 = 'TIP', TIP_1 = 'IMMAT', TIP_2 = 'PROLIF',))


pdf("liana_dotplot_intact.pdf", 20, 10)
liana_test %>%
  liana_dotplot(source_groups = c("MES_DIFF", 'MES_DIFF*','MES_EPINEUR','MES_DIVIDING','MES_PERINEUR','MES_ENDONEUR'),
                target_groups = c('ARTERIAL', 'BARR_END_CAP','CAPILLARY_PLVAP-','CAPILLARY_PLVAP+', 'IMMAT', 'PROLIF','TIP','VENOUS_PLVAP+', 'VENOUS_PLVAP-'),
                ntop = 20)
dev.off()


pdf("liana_dotplot_intact_TIP_vs_MES.pdf", 20, 10)
liana_test %>%
  liana_dotplot(source_groups = c("TIP"),
                target_groups = c("MES_DIFF", 'MES_DIFF*','MES_EPINEUR','MES_DIVIDING','MES_PERINEUR','MES_ENDONEUR'),
                ntop = 20)
dev.off()

pdf("liana_dotplot_intact_PROLIF_vs_MES.pdf", 20, 10)
liana_test %>%
  liana_dotplot(source_groups = c("PROLIF"),
                target_groups = c("MES_DIFF", 'MES_DIFF*','MES_EPINEUR','MES_DIVIDING','MES_PERINEUR','MES_ENDONEUR'),
                ntop = 20)
dev.off()


pdf("liana_dotplot_intact_IMMAT_vs_MES.pdf", 20, 10)
liana_test %>%
  liana_dotplot(source_groups = c("PROLIF"),
                target_groups = c("MES_DIFF", 'MES_DIFF*','MES_EPINEUR','MES_DIVIDING','MES_PERINEUR','MES_ENDONEUR'),
                ntop = 20)
dev.off()



# on a separate notebook I performed the same analysis on the "nerve injury" dataset


liana_test_inj = readRDS("liana_test_1.Rds")

#save injury output

natmi = liana_test_inj$natmi
connectome = liana_test_inj$connectome
logfc =liana_test_inj$logfc
sca = liana_test_inj$sca
cellphonedb= liana_test_inj$cellphonedb
cytotalk = liana_test_inj$cytotalk

write.table(natmi, "natmi_inj.txt", sep = "\t", quote = FALSE)
write.table(connectome, "connectome_inj.txt", sep = "\t", quote = FALSE)
write.table(logfc, "logfc_inj.txt", sep = "\t", quote = FALSE)
write.table(sca, "sca_inj.txt", sep = "\t", quote = FALSE)
write.table(cellphonedb, "cellphonedb_inj.txt", sep = "\t", quote = FALSE)
write.table(cytotalk, "cytotalk_inj.txt", sep = "\t", quote = FALSE)


list_of_datasets <- list("natmi" = natmi, "connectome" = connectome, "logfc" = logfc, "sca" = sca, "cellphonedb" = cellphonedb, "cytotalk" = cytotalk )
write.xlsx(list_of_datasets, file = "liana_test_injury.xlsx")

liana_test_inj = readRDS("liana_test_1_aggregate.Rds")


pdf("liana_dotplot_inj_TIP_vs_MES.pdf", 20, 10)
liana_test_inj %>%
  liana_dotplot(source_groups = c("TIP"),
                target_groups = c("MES_DIFF", 'MES_DIFF*','MES_EPINEUR','MES_DIVIDING','MES_PERINEUR','MES_ENDONEUR'),
                ntop = 20)
dev.off()

pdf("liana_dotplot_inj_PROLIF_vs_MES.pdf", 20, 10)
liana_test_inj %>%
  liana_dotplot(source_groups = c("PROLIF"),
                target_groups = c("MES_DIFF", 'MES_DIFF*','MES_EPINEUR','MES_DIVIDING','MES_PERINEUR','MES_ENDONEUR'),
                ntop = 20)
dev.off()


pdf("liana_dotplot_inj_IMMAT_vs_MES.pdf", 20, 10)
liana_test_inj %>%
  liana_dotplot(source_groups = c("PROLIF"),
                target_groups = c("MES_DIFF", 'MES_DIFF*','MES_EPINEUR','MES_DIVIDING','MES_PERINEUR','MES_ENDONEUR'),
                ntop = 20)
dev.off()



#```

#Note that missing dots are interactions which are not expressed in at least 10% of the cells (by default)
#in both cell clusters.

#In this case, we consider the specificity of interactions as defined by 
#[NATMI's](https://www.nature.com/articles/s41467-020-18873-z) edge specificity weights. 
#NATMI's specificity edges range from 0 to 1, where `1` means both the ligand
#  and receptor are uniquely expressed in a given pair of cells. 
#  Expression magnitude, represented by [SingleCellExperiment's LRScore](https://academic.oup.com/nar/article/48/10/e55/5810485), is on the other hand, 
#meant to represent a non-negative regularized score, comparable between datasets.

## Frequency Heatmap
#We will now plot the frequencies of interactions for each pair of *potentially* communicating cell types.

#This heatmap was inspired by CellChat's and CellPhoneDB's heatmap designs.

#First, we filter interactions by `aggregate_rank`, which can itself
#be treated as a p-value of for the robustly, highly ranked interactions.
#Nevertheless, one can also filter according to [CPDB's p-value](https://www.nature.com/articles/s41596-020-0292-x) (<=0.05) or
#[SingleCellExperiments LRScore](https://academic.oup.com/nar/article/48/10/e55/5810485), etc.

liana_trunc <- liana_test %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected


liana_trunc$
write.table(liana_trunc, "liana_filter_aggregate_rank_intact.txt", sep = "\t", quote = FALSE)

pdf("liana_heat_freq_intact.pdf", 10, 10)
heat_freq(liana_trunc)
dev.off()


liana_trunc_injury = read.table("liana_filter_aggregate_rank_D7.txt")


liana_trunc_inj = readRDS("liana_test_1_aggregate.Rds")
#saveRDS(liana_test, "liana_test_1_aggregate_intact.Rds")


# in the intact we do not have MES_DIFF* so, in order to make a comparable heatmap,
# I'll remove them also from the Injury

liana_trunc_inj_tem = liana_trunc_inj[!liana_trunc_inj$source == "MES_DIFF*", ]
liana_trunc_inj_temp = liana_trunc_inj_tem[!liana_trunc_inj_tem$target == "MES_DIFF*", ]


liana_trunc <- liana_trunc %>% 
  mutate(source = recode(source, TIP_3 = 'TIP', TIP_1 = 'IMMAT', TIP_2 = 'PROLIF',))

liana_trunc <- liana_trunc %>% 
  mutate(target = recode(target, TIP_3 = 'TIP', TIP_1 = 'IMMAT', TIP_2 = 'PROLIF',))



# plot a nice heatmap comparing Intact and Injury LR interactions


gg1 <- heat_freq(liana_trunc) # intact
gg2 <- heat_freq(liana_trunc_inj_temp) #injury


gg1_orig = gg1
gg2_orig = gg2

#----------------------------------------
gg1@matrix[is.na(gg1@matrix)] <- 0
gg2@matrix[is.na(gg2@matrix)] <- 0
#---------------------------------------
matrix_color_mapping <- gg2@matrix_color_mapping
matrix_color_mapping@colors <- colorRampPalette(c('dodgerblue4','white','darkred'))(5)
gg1@matrix_color_mapping <- gg2@matrix_color_mapping <- matrix_color_mapping
#---------------------------------------
# gg1@right_annotation@anno_list$Strength@fun@data_scale <- gg2@right_annotation@anno_list$Strength@fun@data_scale
# gg1@top_annotation@anno_list$Strength@fun@data_scale <- gg2@top_annotation@anno_list$Strength@fun@data_scale
#---------------------------------------
gg1 + gg2




mat <- gg2_orig@matrix
mat[is.na(mat)] <- 0
legend.name <- 'Number of interactions'
color.heatmap.use <- colorRampPalette(c('dodgerblue4','snow','darkred'))(5)
font.size = 8
font.size.title = 10
title.name = NULL
cluster.rows = TRUE; cluster.cols = TRUE
color.use <- c(
  'VENOUS_PLVAP+' = '#F8766D', 
  'VENOUS_PLVAP-' = 'brown',
  'IMMAT' = '#CD9600',
  'PROLIF' = 'pink',
  'BARR_END_CAP' = '#7CAE00',
  'CAPILLARY_PLVAP-' = '#00BE67', 
  'ARTERIAL' = '#00BFC4',
  'MACROPHAGES' = '#00A9FF', 
  'MES_DIFF'='#C97398',
  'MES_DIVIDING'='gold',
  'MES_ENDONEUR'='#7DF3FF',
  'MES_EPINEUR'='#3993D0',
  'MES_PERINEUR'='#4339D0',
  'LEC'='#FF61CC', 
  'CAPILLARY_PLVAP+' = 'grey', 
  'SCHWANN' = '#DCE961',
  'TIP' = "purple"
)[rownames(mat)]

library(ComplexHeatmap)

df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))
row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))

ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)



ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
              show_column_dend = F, show_row_dend = F,
              bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
              cluster_rows = cluster.rows,cluster_columns = cluster.rows,
              row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
              # width = unit(width, "cm"), height = unit(height, "cm"),
              column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
              row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                          border = NA, #at = colorbar.break,
                                          legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
)
(injury_heatmap <- ht1)


mat <- gg1_orig@matrix
mat[is.na(mat)] <- 0
legend.name <- 'Number of interactions'
color.heatmap.use <- colorRampPalette(c('dodgerblue4','white','darkred'))(5)
font.size = 8
font.size.title = 10
title.name = NULL
cluster.rows = FALSE; cluster.cols = FALSE
color.use <- c(
  'VENOUS_PLVAP+' = '#F8766D', 
  'VENOUS_PLVAP-' = 'brown',
  'IMMAT' = '#CD9600',
  'PROLIF' = 'pink',
  'MES_DIVIDING'='gold',
  'BARR_END_CAP' = '#7CAE00',
  'CAPILLARY_PLVAP-' = '#00BE67', 
  'ARTERIAL' = '#00BFC4',
  'MACROPHAGES' = '#00A9FF', 
  'MES_DIFF'='#C97398',
  'MES_ENDONEUR'='#7DF3FF',
  'MES_EPINEUR'='#3993D0',
  'MES_PERINEUR'='#4339D0',
  'CAPILLARY_PLVAP+' = 'grey', 
  'SCHWANN' = '#DCE961',
  'TIP' = "purple"
)[rownames(mat)]



df<- data.frame(group = colnames(mat)); 
rownames(df) <- colnames(mat)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))
row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))

ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                            ylim = injury_heatmap@right_annotation@anno_list$Strength@fun@data_scale, 
                                            border = FALSE, gp = gpar(fill = color.use, col=color.use)), 
                    show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                ylim = injury_heatmap@top_annotation@anno_list$Strength@fun@data_scale, 
                                                border = FALSE,gp = gpar(fill = color.use, col=color.use)), 
                        show_annotation_name = FALSE)


ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
              bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
              cluster_rows = cluster.rows,cluster_columns = cluster.rows,
              row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
              # width = unit(width, "cm"), height = unit(height, "cm"),
              column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
              row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                          border = NA, #at = colorbar.break,
                                          legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
)

ht1@column_order <- c(7,9,6,10,8,14,1,4,12,14,15,16,11,2,3,5)
ht1@row_order <- c(7,9,6,10,8,14,1,4,12,14,15,16,11,2,3,5)
intact_heatmap <- ht1
intact_heatmap@matrix_color_mapping <- injury_heatmap@matrix_color_mapping

ht2

injury_heatmap@column_order
intact_heatmap@column_order

pdf('Liana_interactions_heatmap_clustered_new.pdf', height = 4)
intact_heatmap + injury_heatmap
dev.off()


##########





#Here, we see that `NK`-`CD8 T` share a relatively large number of the inferred
#interactions, with many of those being send from `NK` - note large sum of 
#interactions (gray barplot) in which NK is the `Sender`.

#**NB!** Here, an assumption is implied that the number of interactions inferred between
#cell types is informative of the communication events occurring in the system.
#This is a rather strong assumption that one should consider prior to making any 
#conclusions. [Our suggestion](https://www.biorxiv.org/content/10.1101/2021.05.21.445160v1.abstract) is that any conclusions should be complimented with further information, such as biological prior knowledge, spatial information, etc.


## Frequency Chord diagram

#Here, we will generate a chord diagram equivalent to the frequencies heatmap.

#First, make sure you have the `circlize` package installed.
#```{r cricl, eval=TRUE, message=TRUE}
if(!require("circlize")){
  install.packages("circlize", quiet = TRUE,
                   repos = "http://cran.us.r-project.org")
}
#```

#In this case, one could choose the source and target cell type groups that they wish to plot.
#```{r chordfreq, warning=FALSE, fig.width=6, fig.height=6}
p <- chord_freq(liana_trunc,
                source_groups = c("VENOUS_PLVAP+", "VENOUS_PLVAP-"),
                target_groups = c("MES_DIFF", "MES_DIFF*", "MES_PERINEUR","MES_ENDONEUR", 'MES_EPINEUR'))
#```

?chord_freq


#For more advanced visualization options, we kindly refer the user to [SCPubr](https://enblacar.github.io/SCpubr-book/21-LigandReceptorPlot.html).


## Run any method of choice.

#We will now run only the CellPhoneDB's permutation-based algorithm.
#We will also lower the number of permutations that we wish to perform for the sake of computational time.
#Note that one can also parallelize the CPDB algorithm implemented in LIANA (in this case we don't,
#as this would only make sense when working with large datasets).

#Also, here we will use a `SingleCellExperiment` object as input. In reality, LIANA converts `Seurat` objects
#to SingleCellExperiment and is to a large extend based on the BioConductor single-cell infrastructure.

#```{r cpdb_test, message = FALSE, print = FALSE, warning=FALSE, fig.width=11, fig.height=9}
# Load Sce testdata
sce <- readRDS(file.path(liana_path , "testdata", "input", "testsce.rds"))
# RUN CPDB alone
cpdb_test <- liana_wrap(sce,
                        method = 'cellphonedb',
                        resource = c('CellPhoneDB'),
                        permutation.params = list(nperms=100,
                                                  parallelize=FALSE,
                                                  workers=4),
                        expr_prop=0.05)
# identify interactions of interest
cpdb_int <- cpdb_test %>%
  # only keep interactions with p-val <= 0.05
  filter(pvalue <= 0.05) %>% # this reflects interactions `specificity`
  # then rank according to `magnitude` (lr_mean in this case)
  rank_method(method_name = "cellphonedb",
              mode = "magnitude") %>%
  # keep top 20 interactions (regardless of cell type)
  distinct_at(c("ligand.complex", "receptor.complex")) %>%
  head(20)
# Plot toy results
cpdb_test %>%
  # keep only the interactions of interest
  inner_join(cpdb_int, 
             by = c("ligand.complex", "receptor.complex")) %>%
  # invert size (low p-value/high specificity = larger dot size)
  # + add a small value to avoid Infinity for 0s
  mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
  liana_dotplot(source_groups = c("c"),
                target_groups = c("c", "a", "b"),
                specificity = "pvalue",
                magnitude = "lr.mean",
                show_complex = TRUE,
                size.label = "-log10(p-value)")
#```


## SingleCellSignalR, CytoTalk, NATMI, and Connectome scores /w Complexes
#The re-implementation of the aforementioned methods in LIANA enables us to make use of 
#multimeric complex information as provided by the e.g. CellPhoneDB, CellChatDB, and ICELLNET resources

#```{r complex_test, message = FALSE, print = FALSE}
# Run liana re-implementations with the CellPhoneDB resource
complex_test <- liana_wrap(testdata,
                           method = c('natmi', 'sca', 'logfc'),
                           resource = c('CellPhoneDB'))
complex_test %>% liana_aggregate()
#```


## Call `liana` with overwritten default settings
#By default `liana` is run with the default for each method which can be obtained via `liana_default()`

#Alternatively, one can also overwrite the default settings by simply passing them to the liana wrapper function
#```{r liana_params, message = FALSE, print = FALSE}
# define geometric mean
geometric_mean <- function(vec){exp(mean(log(vec)))}
# Overwrite default parameters by providing a list of parameters
liana_test <- liana_wrap(testdata,
                         method = c('cellphonedb', 'sca'),
                         resource = 'Consensus',
                         permutation.params = 
                           list(
                             nperms = 10 # here we run cpdb it only with 10 permutations
                           ),
                         complex_policy="geometric_mean"
)
# This returns a list of results for each method
#liana_test %>% dplyr::glimpse()
#```
#Note that with `liana_call.params` one can change the way that we account for complexes.
#By default this is set to the mean of the subunits in the complex (as done for expression in CellPhoneDB/Squidpy), 
#and it would return the subunit closest to the mean. 

#Here we change it to the geometric_mean, but any alternative function can be passed and other 
#perfectly viable approaches could also include e.g. the Trimean of CellChat.

#Please refer to the `?liana_wrap` documentation for more information on all the parameters that can be tuned in liana. Also, one can obtain a list with all default parameters by calling the `liana_defaults()` function.


## Citation
#```{r liana citation, warning=FALSE, message = FALSE, echo=FALSE}
#citation("liana")
#```


## Session information
#```{r session_info, echo=FALSE}
#options(width = 120)
#sessioninfo::session_info()
#```


liana_defaults()


library(openxlsx)

setwd("/Users/maurizio.aurora")

inj = readRDS("integrated_inj_sub.rds")
unique(inj$orig.ident)


intact = read.table("liana_filter_aggregate_rank_intact.txt", sep = "\t")

write.xlsx(intact,
           file= "liana_filter_aggregate_rank_intact.xlsx", 
           row.names = F,
           asTable = T)


inj = read.table("liana_filter_aggregate_rank_D7.txt", sep = "\t")

write.xlsx(inj,
           file= "liana_filter_aggregate_rank_D7.xlsx", 
           row.names = F,
           asTable = T)


intact_source = intact[grepl(c("TIP_2"),intact$source), ]
intact_target = intact_source[grepl(c("MES"),intact_source$target), ]



colnames(intact_new)
intact_neww = intact_target[order(intact_target$aggregate_rank),]


unique(intact_neww$target)
intact_neww$LR = paste(intact_neww$ligand.complex,intact_neww$receptor.complex, sep ="->")



intact_neww$LR  <- factor(x = STROMA$cond, levels = c("Ctrl","leuk-A" ,"leuk-B"))


lev =  c("Apln->Grm7","Col18a1->Gpc1", 
        "Psen1->Notch3", "Col4a2->Cd44",   
        "Bsg->Slc16a1", "Col4a2->Itgb5",     
        "Lama5->Sdc1", "Col4a2->Sdc1",   
        "Fam3c->Ffar2","Gnai2->Ednra",   
        "Col4a1->Cd44","Ado->Adora1",      
        "Lamc1->Itga6_Itgb4", "Gnai2->Cnr1",      
        "Lama4->Itga6_Itgb4", "Col4a1->Sdc1",   
        "Col4a2->Itgb5", "Lama5->Itga6_Itgb4",
        "Calr->Lrp1")
hin = head(intact_neww, 20)
hin$LR  <- factor(x = hin$LR, levels = c("Apln->Grm7","Col18a1->Gpc1", 
                                         "Psen1->Notch3", "Col4a2->Cd44",   
                                         "Bsg->Slc16a1", "Col4a2->Itgb5",     
                                         "Lama5->Sdc1", "Col4a2->Sdc1",   
                                         "Fam3c->Ffar2","Gnai2->Ednra",   
                                         "Col4a1->Cd44","Ado->Adora1",      
                                         "Lamc1->Itga6_Itgb4", "Gnai2->Cnr1",      
                                         "Col4a1->Sdc1",
                                         "Lama5->Itga6_Itgb4",
                                         "Calr->Lrp1"))

col.range = c(0,11)
p <- ggplot(hin, aes(target,LR)) +
  # + geom_point() 
  geom_point(mapping = aes_string(size = "natmi.edge_specificity", color = "sca.LRscore")) + 
  scale_color_gradientn(colours = c("blue", "red"))+
  scale_size_continuous(range=c(1,7), name="specificity")
  #scale_size_binned(limits = c(0.1, 0.5), n.breaks = 5)+
  theme(axis.text.x = element_text( color='black', face="bold.italic", size=8)) +
  labs(title = "UP COND_B") 


pdf("pathway_COND_B_vs_COND_A_stroma_sc_fdr_up_dotplot_CAMs_test_6.pdf",9, 7)
p
dev.off()
