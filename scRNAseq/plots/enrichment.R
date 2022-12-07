
library(ggplot2)
library(GEOquery)
library(edgeR)
library(DESeq2)
library(enrichR)
library(grid)
library(gridExtra)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(openxlsx)
library(viridis)
library(stringr)



myCAFs_LEUB_vs_LEUA=FindMarkers(STROMA_cafs_DEG, ident.1 = "myCAFs_LEUB", ident.2 = "myCAFs_LEUA", verbose = FALSE, min.pct = 0.25)
myCAFs_LEUB_vs_LEUA$gene <- rownames(myCAFs_LEUB_vs_LEUA)


filename_xls <- 'myCAFs_LEUB_vs_LEUA.xlsx'
write.xlsx(myCAFs_LEUB_vs_LEUA,
           file= filename_xls, 
           row.names = T,
           asTable = T)


iCAFs_LEUB_vs_CTRL=FindMarkers(STROMA_cafs_DEG, ident.1 = "iCAFs_LEUB", ident.2 = "iCAFs_combined_ctrl", verbose = FALSE, min.pct = 0.25)
iCAFs_LEUB_vs_CTRL$gene <- rownames(iCAFs_LEUB_vs_CTRL)


filename_xls <- 'iCAFs_LEUB_vs_CTRL.xlsx'
write.xlsx(iCAFs_LEUB_vs_CTRL,
           file= filename_xls, 
           row.names = T,
           asTable = T)


apCAFs_LEUB_vs_LEUA=FindMarkers(STROMA_cafs_DEG, ident.1 = "apCAFs_LEUB", ident.2 = "apCAFs_LEUA", verbose = FALSE, min.pct = 0.25)
apCAFs_LEUB_vs_LEUA$gene <- rownames(apCAFs_LEUB_vs_LEUA)


filename_xls <- 'apCAFs_LEUB_vs_LEUA.xlsx'
write.xlsx(apCAFs_LEUB_vs_LEUA,
           file= filename_xls, 
           row.names = T,
           asTable = T)

setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Brendolan/BrendolanA_1190_scRNAseq/9_bioinfo/LEUA_LEUB_merged_ctrls_EYFP+Marker_genes/figures/STROMA/CAFs/CAFs_final/DEGs_logfc0_min.pct0.1")


suppressMessages(library(openxlsx))

apCAFs_LEUA_vs_CTRL = read.xlsx("apCAFs_LEUA_vs_CTRL_logfc0_min.pct0.1.xlsx", rowNames = T)
apCAFs_LEUB_vs_CTRL= read.xlsx("apCAFs_LEUB_vs_CTRL_logfc0_min.pct0.1.xlsx", rowNames = T)
apCAFs_LEUB_vs_LEUA= read.xlsx("apCAFs_LEUB_vs_LEUA_logfc0_min.pct0.1.xlsx", rowNames = T)
iCAFs_LEUA_vs_CTRL = read.xlsx("iCAFs_LEUA_vs_CTRL_logfc0_min.pct0.1.xlsx", rowNames = T)
iCAFs_LEUB_vs_CTRL = read.xlsx("iCAFs_LEUB_vs_CTRL_logfc0_min.pct0.1.xlsx", rowNames = T)
iCAFs_LEUB_vs_LEUA = read.xlsx("iCAFs_LEUB_vs_LEUA_logfc0_min.pct0.1.xlsx", rowNames = T)
myCAFs_LEUA_vs_CTRL = read.xlsx("myCAFs_LEUA_vs_CTRL_logfc0_min.pct0.1.xlsx", rowNames = T)
myCAFs_LEUB_vs_CTRL = read.xlsx("myCAFs_LEUB_vs_CTRL_logfc0_min.pct0.1.xlsx", rowNames = T)
myCAFs_LEUB_vs_LEUA = read.xlsx("myCAFs_LEUB_vs_LEUA_logfc0_min.pct0.1.xlsx", rowNames = T)



dfList=list(apCAFs_LEUA_vs_CTRL,
            apCAFs_LEUB_vs_CTRL,
            apCAFs_LEUB_vs_LEUA,
            iCAFs_LEUA_vs_CTRL,
            iCAFs_LEUB_vs_CTRL,
            iCAFs_LEUB_vs_LEUA,
            myCAFs_LEUA_vs_CTRL,
            myCAFs_LEUB_vs_CTRL,
            myCAFs_LEUB_vs_LEUA)

names(dfList)<-c("apCAFs_LEUA_vs_CTRL",
            "apCAFs_LEUB_vs_CTRL",
            "apCAFs_LEUB_vs_LEUA",
            "iCAFs_LEUA_vs_CTRL",
            "iCAFs_LEUB_vs_CTRL",
            "iCAFs_LEUB_vs_LEUA",
            "myCAFs_LEUA_vs_CTRL",
            "myCAFs_LEUB_vs_CTRL",
            "myCAFs_LEUB_vs_LEUA")




getwd()
dir.create('enrichR_bis/', showWarnings=FALSE, recursive=TRUE)
enrichr.list <- list()
# -------------------------
# enrichment Parameters
# -------------------------
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2021",
                      "GO_Cellular_Component_2021",
                      "GO_Molecular_Function_2021",
                      "KEGG_2021_Human",
                      "MSigDB_Hallmark_2020",
                      "WikiPathways_2016",
                      "BioCarta_2016",
                      "Jensen_TISSUES",
                      "Jensen_COMPARTMENTS",
                      "Jensen_DISEASES")

for (i in names(dfList)) {
  df_obj <- dfList[[i]]
  signif <- (df_obj[df_obj$p_val_adj <= 0.05, ])
  number_of_sig_genes  <- nrow(signif)
  #print(head(signif))
  cat(i, number_of_sig_genes, "significant genes\n")
  neg <- nrow(signif[signif$avg_logFC < 0, ])
  
  neg_list <-  rownames(signif[signif$avg_log2FC < 0, ])
  
  write.table(neg_list, paste('./enrichR_bis/FDRdown_',i,
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  pos  <- nrow(signif[signif$avg_logFC > 0, ])
  pos_list  <- rownames(signif[signif$avg_log2FC > 0, ])
  write.table(pos_list, paste('./enrichR_bis/FDRup_',i,
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
    filename = paste("./enrichR_bis/",i,j,".xlsx", sep="")
    
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)}}


#############################################



### Plot the top N most significant pathways

# Setup parameters

fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))

N=15

enrichR.file = "/Users/maurizio.aurora/Documents/FIGURES_BRENDOLAN/enrichR_bis/myCAFs_LEUA_vs_CTRLfdr_up.xlsx"


s = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "MSigDB_Hallmark_2020")
enrichR.table = data.frame()

for (dat in s) {
  
  Table <- openxlsx::read.xlsx(xlsxFile=enrichR.file, sheet=dat, startRow=1, colNames=T, rowNames=T,
                               
                               detectDates=F, skipEmptyRows=T, skipEmptyCols=T, na.strings="NA", fillMergedCells=F)
  
  enrichR.table = rbind(enrichR.table, Table)
  
}

# Considering only significant pathways

p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < 0.05,])



if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
  
  pathways.dataframe = data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value)
  
  # Formatting the dataframe for the plot
  
  pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
  
  pathways.dataframe$Pathway.num = as.factor(dim(pathways.dataframe)[1]:1) # as.factor is somehow necessary to write the pathway names
  
}

up = head(pathways.dataframe[1:N,], 15)
up$Condition = "UP"



enrichR.file = "/Users/maurizio.aurora/Documents/FIGURES_BRENDOLAN/enrichR_bis/myCAFs_LEUA_vs_CTRLfdr_down.xlsx"

getwd()

s = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "MSigDB_Hallmark_2020")
enrichR.table = data.frame()

for (dat in s) {
  
  Table <- openxlsx::read.xlsx(xlsxFile=enrichR.file, sheet=dat, startRow=1, colNames=T, rowNames=T,
                               
                               detectDates=F, skipEmptyRows=T, skipEmptyCols=T, na.strings="NA", fillMergedCells=F)
  
  enrichR.table = rbind(enrichR.table, Table)
  
}

# Considering only significant pathways

p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < 0.05,])



if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
  
  pathways.dataframe = data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value)
  
  # Formatting the dataframe for the plot
  
  pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
  
  pathways.dataframe$Pathway.num = as.factor(dim(pathways.dataframe)[1]:1) # as.factor is somehow necessary to write the pathway names
  
}


dw = head(pathways.dataframe[1:N,], 15)
dw$Condition = "DW"


merged = rbind(up,dw)
dati = merged

head(dati)
dati$Log10Adj.P.value = -log10(dati$p.value.adj)

#x <- parse(text=dati$Overlap)
#for(i in seq_along(x))
#  dati$gene.ratio[i] <- eval(x[[i]])

colnames(dati)
library(stringr)
dati$Pathway = str_replace(dati$Pathway, "\\s*\\([^\\)]+\\)", " ")

p <- ggplot(dati, aes(Condition,Pathway)) +
  # + geom_point() 
  geom_point(mapping = aes_string(size = "gene.ratio", color = "Log10Adj.P.value")) + 
  scale_color_gradientn(colours = c("blue", "red"), name = "-Log10(p.adj.value)")+
  #scale_size_continuous(range=c(1,7), name="Gene ratio")
  scale_size_binned(limits = c(0.1, 0.7), n.breaks = 5)+
  theme(axis.text.x = element_text( color='black', face="bold.italic", size=14)) +
  labs(title = "myCAFs LEUAvsCTRL") 

pdf("pathway_myCAFs_LEUA_vs_CTRL_stroma_sc.pdf",7, 5)
p
dev.off()

write.xlsx(dati, "pathway_myCAFs_LEUA_vs_CTRL_stroma_sc.xlsx")



