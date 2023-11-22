library(openxlsx)
library(Seurat)
# Set tibble data frame print options
#options(tibble.max_extra_cols = 10)
library(biomaRt)
library(dplyr)
library(fgsea)
library(enrichR)
library(clusterProfiler)


#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)
organism = "org.Mm.eg.db"



setwd("DEG_logfc.threshold_Inf")

# reading in DEG data 
mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A=read.xlsx("mesothelial cells leuk-B vs mesothelial cells leuk-A DGE.xlsx", rowNames = F)  
# reading in data from deseq2

# we want the log2 fold change 
original_gene_list <- mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A$avg_log2FC

# name the vector
names(original_gene_list) <- mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)






library("biomaRt")
library(dplyr)
human_mouse_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

head(mouse_human_genes$Common.Organism.Name)
convert_human_to_mouse <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (human_mouse_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (human_mouse_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    }
  }
  
  return (output)
}


library(msigdbr)
library(fgsea)



msigdbr_df <- msigdbr(species = "Mus musculus", category = "C2")
MARTENS_TRETINOIN_RESPONSE_DN = msigdbr_df[grepl("MARTENS_TRETINOIN_RESPONSE_DN", msigdbr_df$gs_name),]
geneset_MARTENS_TRETINOIN_RESPONSE_DN = MARTENS_TRETINOIN_RESPONSE_DN[,c("gs_name","gene_symbol")]

#MARTENS_TRETINOIN_RESPONSE_DN_cols = MARTENS_TRETINOIN_RESPONSE_DN[,c("gs_name","human_gene_symbol")]
#genes = MARTENS_TRETINOIN_RESPONSE_DN_cols$human_gene_symbol
#converted = convert_human_to_mouse(genes)
#converted_df = as.data.frame(converted)
#colnames(converted_df) = "mouse_gene_symbol"
#converted_df$gs_name = "MARTENS_TRETINOIN_RESPONSE_DN"
#geneset <- converted_df[, c("gs_name","mouse_gene_symbol")]
#geneset_MARTENS_TRETINOIN_RESPONSE_DN <- geneset





setwd("DEG_logfc.threshold_Inf")
getwd()
mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A=read.xlsx("mesothelial cells leuk-B vs mesothelial cells leuk-A DGE.xlsx", rowNames = F)                      
mesothelial_cells_leuk_B_vs_mesothelial_cells_Ctrl=read.xlsx("mesothelial cells leuk-B vs mesothelial cells Ctrl DGE.xlsx", rowNames = F) 
mesothelial_cells_leuk_A_vs_mesothelial_cells_Ctrl=read.xlsx("mesothelial cells leuk-A vs mesothelial cells Ctrl DGE.xlsx", rowNames = F)



################

rownames(mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A) = mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A$gene
rownames(mesothelial_cells_leuk_B_vs_mesothelial_cells_Ctrl) = mesothelial_cells_leuk_B_vs_mesothelial_cells_Ctrl$gene
rownames(mesothelial_cells_leuk_A_vs_mesothelial_cells_Ctrl) = mesothelial_cells_leuk_A_vs_mesothelial_cells_Ctrl$gene



dfList<-list(mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A,
             mesothelial_cells_leuk_B_vs_mesothelial_cells_Ctrl,
             mesothelial_cells_leuk_A_vs_mesothelial_cells_Ctrl)


names(dfList) = c("mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A",
                  "mesothelial_cells_leuk_B_vs_mesothelial_cells_Ctrl",
                  "mesothelial_cells_leuk_A_vs_mesothelial_cells_Ctrl")

#######
dir.create('GSEA/', showWarnings=FALSE, recursive=TRUE)

# -------------------------
# enrichment Parameters
# -------------------------
# databases to make the enrichment of
GSEA.databases <- c("geneset_MARTENS_TRETINOIN_RESPONSE_DN")

for (i in names(dfList)) {
  df_obj <- dfList[[i]]
  # we want the log2 fold change 
  original_gene_list <- df_obj$avg_log2FC
  # name the vector
  names(original_gene_list) <- df_obj$gene
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gsea_results <- GSEA(gene_list, TERM2GENE=geneset_MARTENS_TRETINOIN_RESPONSE_DN, verbose=FALSE, pvalueCutoff = 0.01, maxGSSize = 2000, minGSSize = 5, pAdjustMethod = "BH", eps=0)
  #gseaplot(gsea_results, by = "all", title = "MARTENS_TRETINOIN_RESPONSE_DN", geneSetID = 1)
  pdf(paste('MARTENS_TRETINOIN_RESPONSE_DN_',i,'.pdf',sep=''), width = 6, height = 5)
  gseaplot2(gsea_results,  title = "MARTENS_TRETINOIN_RESPONSE_DN", geneSetID = 1)
  dev.off()
}



#############

i = mesothelial_cells_leuk_B_vs_mesothelial_cells_leuk_A
#for (i in names(dfList)) {
df_obj <- i
# we want the log2 fold change 
original_gene_list <- df_obj$avg_log2FC
# name the vector
names(original_gene_list) <- df_obj$gene
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gsea_results <- GSEA(gene_list, TERM2GENE=geneset_MARTENS_TRETINOIN_RESPONSE_DN, verbose=FALSE, pvalueCutoff = 0.01, maxGSSize = 2000, minGSSize = 5, pAdjustMethod = "BH", eps=0)
#gseaplot(gsea_results, by = "all", title = "MARTENS_TRETINOIN_RESPONSE_DN", geneSetID = 1)
pdf('MARTENS_TRETINOIN_RESPONSE_DN_meso_B_vs_A.pdf', width = 6, height = 5)
gseaplot2(gsea_results,  title = "MARTENS_TRETINOIN_RESPONSE_DN", geneSetID = 1)
dev.off()
