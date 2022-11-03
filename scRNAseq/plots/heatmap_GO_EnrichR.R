
library(openxlsx)

setwd("/Users/maurizio.aurora/Downloads")

pathways = read.xlsx("summary_ontology.xlsx")
pathways = unique(pathways) 



setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Brendolan/BrendolanA_1190_scRNAseq/9_bioinfo/DEG_1190_merged_ctrls_EYFP+/DEG_and_MarkerGenes_TRC_divided_in_2_clusters/EnrichR_bis")



LAvsLB_FDCfdr_up <- read.xlsx("LAvsLB_FDCfdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_FDCfdr_up_sub = subset(LAvsLB_FDCfdr_up, Term %in% pathways$Term)
LAvsLB_FDCfdr_up_non = subset(pathways, !(Term %in% LAvsLB_FDCfdr_up$Term))
LAvsLB_FDCfdr_up_NA = rbindlist(list(LAvsLB_FDCfdr_up_sub, LAvsLB_FDCfdr_up_non), fill = TRUE)
LAvsLB_FDCfdr_up_NA$Log10Adj.P.value = -log10(LAvsLB_FDCfdr_up_NA$Adjusted.P.value)
LAvsLB_FDCfdr_up_0 <- LAvsLB_FDCfdr_up_NA %>% replace(is.na(.), 0)
LAvsLB_FDCfdr_up_CS = LAvsLB_FDCfdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_FDCfdr_up_CS$LEU  = "LA"
LAvsLB_FDCfdr_up_CS$CT  = "FDC"



LAvsLB_FDCfdr_down <- read.xlsx("LAvsLB_FDCfdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_FDCfdr_down_sub = subset(LAvsLB_FDCfdr_down, Term %in% pathways$Term)
LAvsLB_FDCfdr_down_non = subset(pathways, !(Term %in% LAvsLB_FDCfdr_down$Term))
LAvsLB_FDCfdr_down_NA = rbindlist(list(LAvsLB_FDCfdr_down_sub, LAvsLB_FDCfdr_down_non), fill = TRUE)
LAvsLB_FDCfdr_down_NA$Log10Adj.P.value = -log10(LAvsLB_FDCfdr_down_NA$Adjusted.P.value)
LAvsLB_FDCfdr_down_0 <- LAvsLB_FDCfdr_down_NA %>% replace(is.na(.), 0)
LAvsLB_FDCfdr_down_CS = LAvsLB_FDCfdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_FDCfdr_down_CS$LEU  = "LB"
LAvsLB_FDCfdr_down_CS$CT  = "FDC"



################


LAvsLB_PRCfdr_up <- read.xlsx("LAvsLB_PRCfdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_PRCfdr_up_sub = subset(LAvsLB_PRCfdr_up, Term %in% pathways$Term)
LAvsLB_PRCfdr_up_non = subset(pathways, !(Term %in% LAvsLB_PRCfdr_up$Term))
LAvsLB_PRCfdr_up_NA = rbindlist(list(LAvsLB_PRCfdr_up_sub, LAvsLB_PRCfdr_up_non), fill = TRUE)
LAvsLB_PRCfdr_up_NA$Log10Adj.P.value = -log10(LAvsLB_PRCfdr_up_NA$Adjusted.P.value)
LAvsLB_PRCfdr_up_0 <- LAvsLB_PRCfdr_up_NA %>% replace(is.na(.), 0)
LAvsLB_PRCfdr_up_CS = LAvsLB_PRCfdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_PRCfdr_up_CS$LEU  = "LA"
LAvsLB_PRCfdr_up_CS$CT  = "PRC"



LAvsLB_PRCfdr_down <- read.xlsx("LAvsLB_PRCfdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_PRCfdr_down_sub = subset(LAvsLB_PRCfdr_down, Term %in% pathways$Term)
LAvsLB_PRCfdr_down_non = subset(pathways, !(Term %in% LAvsLB_PRCfdr_down$Term))
LAvsLB_PRCfdr_down_NA = rbindlist(list(LAvsLB_PRCfdr_down_sub, LAvsLB_PRCfdr_down_non), fill = TRUE)
LAvsLB_PRCfdr_down_NA$Log10Adj.P.value = -log10(LAvsLB_PRCfdr_down_NA$Adjusted.P.value)
LAvsLB_PRCfdr_down_0 <- LAvsLB_PRCfdr_down_NA %>% replace(is.na(.), 0)
LAvsLB_PRCfdr_down_CS = LAvsLB_PRCfdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_PRCfdr_down_CS$LEU  = "LB"
LAvsLB_PRCfdr_down_CS$CT  = "PRC"


LAvsLB_MRCfdr_up <- read.xlsx("LAvsLB_MRCfdr_up.xlsx", sheet = "GO_Biological_Process_2021") 
LAvsLB_MRCfdr_up_sub = subset(LAvsLB_MRCfdr_up, Term %in% pathways$Term)
LAvsLB_MRCfdr_up_non = subset(pathways, !(Term %in% LAvsLB_MRCfdr_up$Term))
LAvsLB_MRCfdr_up_NA = rbindlist(list(LAvsLB_MRCfdr_up_sub, LAvsLB_MRCfdr_up_non), fill = TRUE)
LAvsLB_MRCfdr_up_NA$Log10Adj.P.value = -log10(LAvsLB_MRCfdr_up_NA$Adjusted.P.value)
LAvsLB_MRCfdr_up_0 <- LAvsLB_MRCfdr_up_NA %>% replace(is.na(.), 0)
LAvsLB_MRCfdr_up_CS = LAvsLB_MRCfdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_MRCfdr_up_CS$LEU  = "LA"
LAvsLB_MRCfdr_up_CS$CT  = "MRC"



LAvsLB_MRCfdr_down <- read.xlsx("LAvsLB_MRCfdr_down.xlsx", sheet = "GO_Biological_Process_2021") 
LAvsLB_MRCfdr_down_sub = subset(LAvsLB_MRCfdr_down, Term %in% pathways$Term)
LAvsLB_MRCfdr_down_non = subset(pathways, !(Term %in% LAvsLB_MRCfdr_down$Term))
LAvsLB_MRCfdr_down_NA = rbindlist(list(LAvsLB_MRCfdr_down_sub, LAvsLB_MRCfdr_down_non), fill = TRUE)
LAvsLB_MRCfdr_down_NA$Log10Adj.P.value = -log10(LAvsLB_MRCfdr_down_NA$Adjusted.P.value)
LAvsLB_MRCfdr_down_0 <- LAvsLB_MRCfdr_down_NA %>% replace(is.na(.), 0)
LAvsLB_MRCfdr_down_CS = LAvsLB_MRCfdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_MRCfdr_down_CS$LEU  = "LB"
LAvsLB_MRCfdr_down_CS$CT  = "MRC"




LAvsLB_TNC1fdr_up <- read.xlsx("LAvsLB_TNC1fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TNC1fdr_up_sub = subset(LAvsLB_TNC1fdr_up, Term %in% pathways$Term)
LAvsLB_TNC1fdr_up_non = subset(pathways, !(Term %in% LAvsLB_TNC1fdr_up$Term))
LAvsLB_TNC1fdr_up_NA = rbindlist(list(LAvsLB_TNC1fdr_up_sub, LAvsLB_TNC1fdr_up_non), fill = TRUE)
LAvsLB_TNC1fdr_up_NA$Log10Adj.P.value = -log10(LAvsLB_TNC1fdr_up_NA$Adjusted.P.value)
LAvsLB_TNC1fdr_up_0 <- LAvsLB_TNC1fdr_up_NA %>% replace(is.na(.), 0)
LAvsLB_TNC1fdr_up_CS = LAvsLB_TNC1fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TNC1fdr_up_CS$LEU  = "LA"
LAvsLB_TNC1fdr_up_CS$CT  = "TNC1"


LAvsLB_TNC1fdr_down <- read.xlsx("LAvsLB_TNC1fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TNC1fdr_down_sub = subset(LAvsLB_TNC1fdr_down, Term %in% pathways$Term)
LAvsLB_TNC1fdr_down_non = subset(pathways, !(Term %in% LAvsLB_TNC1fdr_down$Term))
LAvsLB_TNC1fdr_down_NA = rbindlist(list(LAvsLB_TNC1fdr_down_sub, LAvsLB_TNC1fdr_down_non), fill = TRUE)
LAvsLB_TNC1fdr_down_NA$Log10Adj.P.value = -log10(LAvsLB_TNC1fdr_down_NA$Adjusted.P.value)
LAvsLB_TNC1fdr_down_0 <- LAvsLB_TNC1fdr_down_NA %>% replace(is.na(.), 0)
LAvsLB_TNC1fdr_down_CS = LAvsLB_TNC1fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TNC1fdr_down_CS$LEU  = "LB"
LAvsLB_TNC1fdr_down_CS$CT  = "TNC1"



LAvsLB_TNC2fdr_up <- read.xlsx("LAvsLB_TNC2fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TNC2fdr_up_sub = subset(LAvsLB_TNC2fdr_up, Term %in% pathways$Term)
LAvsLB_TNC2fdr_up_non = subset(pathways, !(Term %in% LAvsLB_TNC2fdr_up$Term))
LAvsLB_TNC2fdr_up_NA = rbindlist(list(LAvsLB_TNC2fdr_up_sub, LAvsLB_TNC2fdr_up_non), fill = TRUE)
LAvsLB_TNC2fdr_up_NA$Log10Adj.P.value = -log10(LAvsLB_TNC2fdr_up_NA$Adjusted.P.value)
LAvsLB_TNC2fdr_up_0 <- LAvsLB_TNC2fdr_up_NA %>% replace(is.na(.), 0)
LAvsLB_TNC2fdr_up_CS = LAvsLB_TNC2fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TNC2fdr_up_CS$LEU  = "LA"
LAvsLB_TNC2fdr_up_CS$CT  = "TNC2"


LAvsLB_TNC2fdr_down <- read.xlsx("LAvsLB_TNC2fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TNC2fdr_down_sub = subset(LAvsLB_TNC2fdr_down, Term %in% pathways$Term)
LAvsLB_TNC2fdr_down_non = subset(pathways, !(Term %in% LAvsLB_TNC2fdr_down$Term))
LAvsLB_TNC2fdr_down_NA = rbindlist(list(LAvsLB_TNC2fdr_down_sub, LAvsLB_TNC2fdr_down_non), fill = TRUE)
LAvsLB_TNC2fdr_down_NA$Log10Adj.P.value = -log10(LAvsLB_TNC2fdr_down_NA$Adjusted.P.value)
LAvsLB_TNC2fdr_down_0 <- LAvsLB_TNC2fdr_down_NA %>% replace(is.na(.), 0)
LAvsLB_TNC2fdr_down_CS = LAvsLB_TNC2fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TNC2fdr_down_CS$LEU  = "LB"
LAvsLB_TNC2fdr_down_CS$CT  = "TNC2"

LAvsLB_TRC1fdr_up <- read.xlsx("LAvsLB_TRC1fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TRC1fdr_up_sub = subset(LAvsLB_TRC1fdr_up, Term %in% pathways$Term)
LAvsLB_TRC1fdr_up_non = subset(pathways, !(Term %in% LAvsLB_TRC1fdr_up$Term))
LAvsLB_TRC1fdr_up_NA = rbindlist(list(LAvsLB_TRC1fdr_up_sub, LAvsLB_TRC1fdr_up_non), fill = TRUE)
LAvsLB_TRC1fdr_up_NA$Log10Adj.P.value = -log10(LAvsLB_TRC1fdr_up_NA$Adjusted.P.value)
LAvsLB_TRC1fdr_up_0 <- LAvsLB_TRC1fdr_up_NA %>% replace(is.na(.), 0)
LAvsLB_TRC1fdr_up_CS = LAvsLB_TRC1fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TRC1fdr_up_CS$LEU  = "LA"
LAvsLB_TRC1fdr_up_CS$CT  = "TRC1"


LAvsLB_TRC1fdr_down <- read.xlsx("LAvsLB_TRC1fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TRC1fdr_down_sub = subset(LAvsLB_TRC1fdr_down, Term %in% pathways$Term)
LAvsLB_TRC1fdr_down_non = subset(pathways, !(Term %in% LAvsLB_TRC1fdr_down$Term))
LAvsLB_TRC1fdr_down_NA = rbindlist(list(LAvsLB_TRC1fdr_down_sub, LAvsLB_TRC1fdr_down_non), fill = TRUE)
LAvsLB_TRC1fdr_down_NA$Log10Adj.P.value = -log10(LAvsLB_TRC1fdr_down_NA$Adjusted.P.value)
LAvsLB_TRC1fdr_down_0 <- LAvsLB_TRC1fdr_down_NA %>% replace(is.na(.), 0)
LAvsLB_TRC1fdr_down_CS = LAvsLB_TRC1fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TRC1fdr_down_CS$LEU  = "LB"
LAvsLB_TRC1fdr_down_CS$CT  = "TRC1"

LAvsLB_TRC2fdr_up <- read.xlsx("LAvsLB_TRC2fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TRC2fdr_up_sub = subset(LAvsLB_TRC2fdr_up, Term %in% pathways$Term)
LAvsLB_TRC2fdr_up_non = subset(pathways, !(Term %in% LAvsLB_TRC2fdr_up$Term))
LAvsLB_TRC2fdr_up_NA = rbindlist(list(LAvsLB_TRC2fdr_up_sub, LAvsLB_TRC2fdr_up_non), fill = TRUE)
LAvsLB_TRC2fdr_up_NA$Log10Adj.P.value = -log10(LAvsLB_TRC2fdr_up_NA$Adjusted.P.value)
LAvsLB_TRC2fdr_up_0 <- LAvsLB_TRC2fdr_up_NA %>% replace(is.na(.), 0)
LAvsLB_TRC2fdr_up_CS = LAvsLB_TRC2fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TRC2fdr_up_CS$LEU  = "LA"
LAvsLB_TRC2fdr_up_CS$CT  = "TRC2"


LAvsLB_TRC2fdr_down <- read.xlsx("LAvsLB_TRC2fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LAvsLB_TRC2fdr_down_sub = subset(LAvsLB_TRC2fdr_down, Term %in% pathways$Term)
LAvsLB_TRC2fdr_down_non = subset(pathways, !(Term %in% LAvsLB_TRC2fdr_down$Term))
LAvsLB_TRC2fdr_down_NA = rbindlist(list(LAvsLB_TRC2fdr_down_sub, LAvsLB_TRC2fdr_down_non), fill = TRUE)
LAvsLB_TRC2fdr_down_NA$Log10Adj.P.value = -log10(LAvsLB_TRC2fdr_down_NA$Adjusted.P.value)
LAvsLB_TRC2fdr_down_0 <- LAvsLB_TRC2fdr_down_NA %>% replace(is.na(.), 0)
LAvsLB_TRC2fdr_down_CS = LAvsLB_TRC2fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsLB_TRC2fdr_down_CS$LEU  = "LB"
LAvsLB_TRC2fdr_down_CS$CT  = "TRC2"


bound = rbind(LAvsLB_FDCfdr_up_CS, LAvsLB_FDCfdr_down_CS,
              LAvsLB_MRCfdr_up_CS, LAvsLB_MRCfdr_down_CS,
              LAvsLB_PRCfdr_up_CS, LAvsLB_PRCfdr_down_CS,
              LAvsLB_TRC1fdr_up_CS, LAvsLB_TRC1fdr_down_CS, 
              LAvsLB_TRC2fdr_up_CS, LAvsLB_TRC2fdr_down_CS,
              LAvsLB_TNC1fdr_up_CS, LAvsLB_TNC1fdr_down_CS,
              LAvsLB_TNC2fdr_up_CS, LAvsLB_TNC2fdr_down_CS)



bound$Cond = paste(bound$LEU, bound$CT)
bound = unique(as.data.frame(bound))
nrow(bound)

final = acast(unique(bound), Term ~ Cond, value.var = 'Log10Adj.P.value')
ciaooo <- final %>% replace(is.na(.), 0)
ciaooo = as.data.frame(ciaooo)
ciaooone <- mutate_all(ciaooo, function(x) as.numeric(as.character(x)))
df = ciaooone[rowSums(ciaooone[])>0,]
nrow(df)

head(bound)
annotation_column = as.data.frame(unique(bound[,c('Cond','LEU', 'CT')]))
rownames(annotation_column) <- annotation_column$Cond
head(annotation_column)


category_order = c("ECM", "WNT", "Growth factor beta", "Proliferation", "Chemotaxis", "TNF", "IFN", "Hypoxia", "Cytokine", "Migration", "Lipopolysaccharide", "T cell")



ann_colors = list(
  PATHWAY = c("ECM" = "#FF007F", 
              "WNT" = "#6666FF",
              "Proliferation" = "#CC0066",
              "Growth factor beta" ="#FFFF99", 
              "Chemotaxis" ="pink",
              "T cell" = "#CCFFFF",
              "TNF" = "#660066",
              "IFN" = "#FF0000", 
              "Hypoxia" = "#FF9933",
              "Cytokine" = "darkviolet",
              "Migration" = "#3399FF",
              "Lipopolysaccharide" = "#E5CCFF"),
  LEU = c("LA" = '#4C9900', 
          "LB" = '#3333FF'),
  CT = c('FDC' = '#F564E3', 
         'MRC' = '#F8766D', 
         'PRC' = '#00BFC4', 
         'TRC1' = '#619CFF', 
         'TRC2' = '#BF86BF',
         'TNC1' = '#B79F00', 
         'TNC2' = '#00BA38'))





annotation_column = annotation_column[,c("LEU","CT")]


df_row <- data.frame(pathways)
rownames(df_row) = df_row$Term

#reorder df columns as in annotation_column
df <- df[,match(rownames(annotation_column), colnames(df))]


#reorder by PATHWAY the rownames
df_row = df_row[order(df_row$PATHWAY),]
#or by custom order


df_row <- df_row %>% 
  arrange(factor(PATHWAY, levels = category_order))



df_row = subset(df_row, Term %in% rownames(df))


df_1 <- df[match(df_row$Term, rownames(df)), ]  

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
minH = -1; maxH=1
myb = seq(minH, maxH, by = 0.01)
#crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))

#remove the description in brackets (GO:0071346)
rownames(df_1) = str_replace(rownames(df_1), "\\s*\\([^\\)]+\\)", " ")


pdf("LA_vs_LB_heatmapGO.pdf")
ComplexHeatmap::pheatmap(as.matrix(df_1), 
                         show_rownames = T, 
                         show_colnames = F,
                         cellheight = 10,
                         cellwidth = 10,
                         cluster_rows = F,  
                         cluster_cols = F, 
                         #height = 5, 
                         annotation_col = annotation_column, 
                         annotation_colors = ann_colors, 
                         scale = 'row', 
                         color = myc, 
                         annotation_row = df_row[2], 
                         fontsize_row = 4)
dev.off()