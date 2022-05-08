suppressMessages(library(ggbeeswarm))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
library(pheatmap)
library(reshape2)
library(scales)
library(viridis)
library(SeuratWrappers)
library(slingshot)
require(BiocStyle)
library(SingleCellExperiment)
#print(version) R version 3.6.1
#Bonanomi
packageVersion("Seurat")
library(dittoSeq)



integrated_new = readRDS("integrated_new.RDS")

table(Idents(integrated_new))

DimPlot(integrated_new)

testM = integrated_new

DefaultAssay(testM) = "RNA"

testM$celltype.stim <- paste(Idents(testM), testM$stim, sep = "_")
testM$celltype <- Idents(testM)
Idents(testM) <- "celltype.stim"

table(Idents(testM))

getwd()
enrichr.list <- list()
setwd("/Users/maurizio.aurora/")
dir.create('newdirD/', showWarnings=FALSE, recursive=TRUE)
dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
setwd("/Users/maurizio.aurora/newdirD")
getwd()
for (i in c("A", "B","C")) {
  for (j in c("1", "2")){
    for (k in c("1","2")){
      if (j !=k ){
        m=paste(i,j,"vs",i,k, sep="_")
        temp=FindMarkers(testM, ident.1 = paste(i,j,sep="_"), ident.2 = paste(i,k,sep="_"), verbose = FALSE) 
        temp$gene <- rownames(temp)
        filename_xls = paste(i,j,"vs",i,k,"DGE.xlsx", sep="_")
        write.xlsx(temp,
                   file= filename_xls, 
                   row.names = T,
                   asTable = T)}}}}}
        




file_list <- list.files(path="/Users/maurizio.aurora/newdirD")

A_1_vs_A_2_DGE=read.xlsx("A_1_vs_A_2_DGE.xlsx", rowNames = T)                      
B_1_vs_B_2_DGE=read.xlsx("B_1_vs_B_2_DGE.xlsx", rowNames = T) 
C_1_vs_C_2_DGE=read.xlsx("C_1_vs_C_2_DGE.xlsx", rowNames = T)

                  

dfList=list(A_1_vs_A_2_DGE,
            B_1_vs_B_2_DGE,
            C_1_vs_C_2_DGE)



names(dfList) = c("A_1_vs_A_2_DGE",
                  "B_1_vs_B_2_DGE",
                  "C_1_vs_C_2_DGE")

getwd()

dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
enrichr.list <- list()
# -------------------------
# enrichment Parameters
# -------------------------
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016",
                      "Jensen_TISSUES",
                      "Jensen_COMPARTMENTS",
                      "Jensen_DISEASES")

for (i in names(dfList)) {
  df_obj <- dfList[[i]]
  signif <- (df_obj[df_obj$p_val_adj <= 0.05, ])
  number_of_sig_genes  <- nrow(signif)
  
  cat(i, number_of_sig_genes, "significant genes\n")
  neg <- nrow(signif[signif$avg_log2FC < 0, ])
  
  neg_list <-  rownames(signif[signif$avg_log2FC < 0, ])
  
  write.table(neg_list, paste('./enrichR/FDRdown_',i,
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  pos  <- nrow(signif[signif$avg_log2FC > 0, ])
  pos_list  <- row.names(signif[signif$avg_log2FC > 0, ])
  write.table(pos_list, paste('./enrichR/FDRup_',i,
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  cat(i, pos, "positive fold change\n")
  print(pos_list)
  cat(i, neg, "negative fold change\n")
  print(neg_list)
  
  #enrichr.list <- list()
  enrichr.list[[i]] <- lapply(list(pos_list,neg_list),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
    
  })
  print("NAAAAAAMMEESSS")
  #print((enrichr.list[[i]]))
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down")
  print(enrichr.list[[i]])
  
  
}


for (i in names(dfList)) {
  for (j in c("fdr_up","fdr_down")){
    filename = paste("./enrichR/",i,j,".xlsx", sep="")
    
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)}}

A_1_vs_A_2_DGE=read.xlsx("A_1_vs_A_2_DGE.xlsx", rowNames = T)                      
B_1_vs_B_2_DGE=read.xlsx("B_1_vs_B_2_DGE.xlsx", rowNames = T) 
C_1_vs_C_2_DGE=read.xlsx("C_1_vs_C_2_DGE.xlsx", rowNames = T)


DefaultAssay(integrated_new) = "RNA"

library(cowplot)
A <- subset(integrated_new, idents = "A")
Idents(A) <- "stim"
A <- log1p(AverageExpression(A, verbose = FALSE)$RNA)
A = as.data.frame(A)
A$gene <- rownames(A)
head(A)
colnames(A) = c("intact", "injury","gene")
pval = subset(A_1_vs_A_2_DGE, p_val_adj < 0.05)
pval1 <- pval[order(pval$avg_log2FC),]

up10 <- pval1 %>% top_n(n = 10, wt = avg_log2FC)
dw10 <- pval1 %>% top_n(n = -10, wt = avg_log2FC)

up = up10$gene
down = dw10$gene

cc = c(up,down)

pdf("bisett_1_vs_2_A.pdf", 5,5)
plot=ggplot(A, aes(intact,injury), label=genes) + geom_point() + ggtitle("A") + theme(
  plot.title = element_text(color="black", size=14, face="bold.italic"))
LabelPoints(plot = plot, points = cc, color="red", repel = TRUE)  
dev.off()



B <- subset(integrated_new, idents = "B")
Idents(B) <- "stim"
B <- log1p(AverageExpression(B, verbose = FALSE)$RNA)
B = as.data.frame(B)
B$gene <- rownames(B)
head(B)
colnames(B) = c("intact", "injury","gene")
pval = subset(B_1_vs_B_2_DGE, p_val_adj < 0.05)

pval1 <- pval[order(pval$avg_log2FC),]
up10 <- pval1 %>% top_n(n = 10, wt = avg_log2FC)
dw10 <- pval1 %>% top_n(n = -10, wt = avg_log2FC)

up = up10$gene
down = dw10$gene

cc = c(up,down)

pdf("bisett_1_vs_2_B.pdf", 5,5)
plot=ggplot(B, aes(intact,injury), label=genes) + geom_point() + ggtitle("B") + theme(
  plot.title = element_text(color="black", size=14, face="bold.italic"))
LabelPoints(plot = plot, points = cc, color="red", repel = TRUE)  
dev.off()



C <- subset(integrated_new, idents = "C")
Idents(C) <- "stim"
C <- log1p(AverageExpression(C, verbose = FALSE)$RNA)
C = as.data.frame(C)
C$gene <- rownames(C)
head(C)
colnames(C) = c("intact", "injury","gene")
pval = subset(C_1_vs_C_2_DGE, p_val_adj < 0.05)
pval1 <- pval[order(pval$avg_log2FC),] 
up10 <- pval1 %>% top_n(n = 10, wt = avg_log2FC)
dw10 <- pval1 %>% top_n(n = -10, wt = avg_log2FC)

up = up10$gene
down = dw10$gene

cc = c(up,down)

pdf("bisett_1_vs_2_C.pdf", 5,5)
plot=ggplot(C, aes(intact,injury), label=genes) + geom_point() + ggtitle("C") + theme(
  plot.title = element_text(color="black", size=14, face="bold.italic"))
LabelPoints(plot = plot, points = cc, color="red", repel = TRUE, max.overlaps = Inf, box.padding = 0.5)
dev.off()




