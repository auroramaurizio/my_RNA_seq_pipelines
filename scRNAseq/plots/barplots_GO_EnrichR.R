
barplot = c("cellular response to type I interferon",
            "interferon-gamma-mediated signaling pathway",
            "regulation of response to interferon.gamma",
            "response to interferon.beta",
            "regulation of insulin-like growth factor receptor singling pathway",
            "T cell chemotaxis",
            "T cell migration",
            "Wnt signaling pathway, planar cell polarity pathway",
            "extracellular matrix disassembly",
            "extracellular matrix organization",
            "cell.matrix adhesion",
            "cellular response to cytokine stimulus",
            "cellular response to hypoxia",
            "chemokine-mediated signaling pathway",
            "cytokine-mediated signaling pathway",
            "cotranslational protein targeting to membrane")


barplot = c("cellular response to type I interferon",
            "interferon-gamma-mediated signaling pathway",
            "T cell chemotaxis",
            "T cell migration",
            "Wnt signaling pathway, planar cell polarity pathway",
            "extracellular matrix disassembly",
            "extracellular matrix organization",
            "cellular response to hypoxia",
            "chemokine-mediated signaling pathway",
            "cellular response to cytokine stimulus",
            "cytokine-mediated signaling pathway")

barplott = as.data.frame(barplot)
colnames(barplott) = c("Pathways")

setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Brendolan/BrendolanA_1190_scRNAseq/9_bioinfo/DEG_1190_merged_ctrls_EYFP+/DEG_and_MarkerGenes_TRC_divided_in_2_clusters/EnrichR_bis")

LAvsCTRLs_FDCfdr_up <- read.xlsx("LAvsCTRLs_FDCfdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsCTRLs_FDCfdr_up$Pathways = str_replace(LAvsCTRLs_FDCfdr_up$Term, "\\s*\\([^\\)]+\\)", "")
LAvsCTRLs_FDCfdr_up_sub = subset(LAvsCTRLs_FDCfdr_up, Pathways %in% barplot)
LAvsCTRLs_FDCfdr_up_non = subset(barplott, !(Pathways %in% LAvsCTRLs_FDCfdr_up$Pathways))
LAvsCTRLs_FDCfdr_up_NA = rbindlist(list(LAvsCTRLs_FDCfdr_up_sub, LAvsCTRLs_FDCfdr_up_non), fill = TRUE)
LAvsCTRLs_FDCfdr_up_NA$Log10Adj.P.value = -log10(LAvsCTRLs_FDCfdr_up_NA$Adjusted.P.value)
LAvsCTRLs_FDCfdr_up_0 <- LAvsCTRLs_FDCfdr_up_NA %>% replace(is.na(.), 0)
LAvsCTRLs_FDCfdr_up_CS = LAvsCTRLs_FDCfdr_up_0[,c("Pathways","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsCTRLs_FDCfdr_up_CS$LEU  = "LA"
LAvsCTRLs_FDCfdr_up_CS$CT  = "FDC"

nrow(unique(LAvsCTRLs_FDCfdr_up))






LAvsCTRLs_PRCfdr_up <- read.xlsx("LAvsCTRLs_PRCfdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsCTRLs_PRCfdr_up$Pathways = str_replace(LAvsCTRLs_PRCfdr_up$Term, "\\s*\\([^\\)]+\\)", "")
LAvsCTRLs_PRCfdr_up_sub = subset(LAvsCTRLs_PRCfdr_up, Pathways %in% barplot)
LAvsCTRLs_PRCfdr_up_non = subset(barplott, !(Pathways %in% LAvsCTRLs_PRCfdr_up$Pathways))
LAvsCTRLs_PRCfdr_up_NA = rbindlist(list(LAvsCTRLs_PRCfdr_up_sub, LAvsCTRLs_PRCfdr_up_non), fill = TRUE)
LAvsCTRLs_PRCfdr_up_NA$Log10Adj.P.value = -log10(LAvsCTRLs_PRCfdr_up_NA$Adjusted.P.value)
LAvsCTRLs_PRCfdr_up_0 <- LAvsCTRLs_PRCfdr_up_NA %>% replace(is.na(.), 0)
LAvsCTRLs_PRCfdr_up_CS = LAvsCTRLs_PRCfdr_up_0[,c("Pathways","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsCTRLs_PRCfdr_up_CS$LEU  = "LA"
LAvsCTRLs_PRCfdr_up_CS$CT  = "PRC"



#LAvsCTRLs_MRCfdr_up <- read.xlsx("LAvsCTRLs_MRCfdr_up.xlsx", sheet = "GO_Biological_Process_2021") # empty DF
#LAvsCTRLs_MRCfdr_up_sub = subset(LAvsCTRLs_MRCfdr_up, Term %in% barplott)
#LAvsCTRLs_MRCfdr_up_sub = barplott
#LAvsCTRLs_MRCfdr_up_non = barplott
####LAvsCTRLs_MRCfdr_up_non[setdiff(LAvsCTRLs_MRCfdr_up_non$Term, LAvsCTRLs_MRCfdr_up_sub$Term)] <- 0
LAvsCTRLs_MRCfdr_up_CS = barplott
LAvsCTRLs_MRCfdr_up_CS$Combined.Score = 0
LAvsCTRLs_MRCfdr_up_CS$Log10Adj.P.value = 0
LAvsCTRLs_MRCfdr_up_CS$Overlap = 0
LAvsCTRLs_MRCfdr_up_CS$LEU  = "LA"
LAvsCTRLs_MRCfdr_up_CS$CT  = "MRC"
#LAvsCTRLs_MRCfdr_up_CS = subset(LAvsCTRLs_MRCfdr_up_CS, select = -c(2) )

head(LAvsCTRLs_MRCfdr_up_CS)


LAvsCTRLs_TNC1fdr_up <- read.xlsx("LAvsCTRLs_TNC1fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsCTRLs_TNC1fdr_up$Pathways = str_replace(LAvsCTRLs_TNC1fdr_up$Term, "\\s*\\([^\\)]+\\)", "")
LAvsCTRLs_TNC1fdr_up_sub = subset(LAvsCTRLs_TNC1fdr_up, Pathways %in% barplot)
LAvsCTRLs_TNC1fdr_up_non = subset(barplott, !(Pathways %in% LAvsCTRLs_TNC1fdr_up$Pathways))
LAvsCTRLs_TNC1fdr_up_NA = rbindlist(list(LAvsCTRLs_TNC1fdr_up_sub, LAvsCTRLs_TNC1fdr_up_non), fill = TRUE)
LAvsCTRLs_TNC1fdr_up_NA$Log10Adj.P.value = -log10(LAvsCTRLs_TNC1fdr_up_NA$Adjusted.P.value)
LAvsCTRLs_TNC1fdr_up_0 <- LAvsCTRLs_TNC1fdr_up_NA %>% replace(is.na(.), 0)
LAvsCTRLs_TNC1fdr_up_CS = LAvsCTRLs_TNC1fdr_up_0[,c("Pathways","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsCTRLs_TNC1fdr_up_CS$LEU  = "LA"
LAvsCTRLs_TNC1fdr_up_CS$CT  = "TNC1"



LAvsCTRLs_TNC2fdr_up <- read.xlsx("LAvsCTRLs_TNC2fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsCTRLs_TNC2fdr_up$Pathways = str_replace(LAvsCTRLs_TNC2fdr_up$Term, "\\s*\\([^\\)]+\\)", "")
LAvsCTRLs_TNC2fdr_up_sub = subset(LAvsCTRLs_TNC2fdr_up, Pathways %in% barplot)
LAvsCTRLs_TNC2fdr_up_non = subset(barplott, !(Pathways %in% LAvsCTRLs_TNC2fdr_up$Pathways))
LAvsCTRLs_TNC2fdr_up_NA = rbindlist(list(LAvsCTRLs_TNC2fdr_up_sub, LAvsCTRLs_TNC2fdr_up_non), fill = TRUE)
LAvsCTRLs_TNC2fdr_up_NA$Log10Adj.P.value = -log10(LAvsCTRLs_TNC2fdr_up_NA$Adjusted.P.value)
LAvsCTRLs_TNC2fdr_up_0 <- LAvsCTRLs_TNC2fdr_up_NA %>% replace(is.na(.), 0)
LAvsCTRLs_TNC2fdr_up_CS = LAvsCTRLs_TNC2fdr_up_0[,c("Pathways","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsCTRLs_TNC2fdr_up_CS$LEU  = "LA"
LAvsCTRLs_TNC2fdr_up_CS$CT  = "TNC2"



LAvsCTRLs_TRC1fdr_up <- read.xlsx("LAvsCTRLs_TRC1fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsCTRLs_TRC1fdr_up$Pathways = str_replace(LAvsCTRLs_TRC1fdr_up$Term, "\\s*\\([^\\)]+\\)", "")
LAvsCTRLs_TRC1fdr_up_sub = subset(LAvsCTRLs_TRC1fdr_up, Pathways %in% barplot)
LAvsCTRLs_TRC1fdr_up_non = subset(barplott, !(Pathways %in% LAvsCTRLs_TRC1fdr_up$Pathways))
LAvsCTRLs_TRC1fdr_up_NA = rbindlist(list(LAvsCTRLs_TRC1fdr_up_sub, LAvsCTRLs_TRC1fdr_up_non), fill = TRUE)
LAvsCTRLs_TRC1fdr_up_NA$Log10Adj.P.value = -log10(LAvsCTRLs_TRC1fdr_up_NA$Adjusted.P.value)
LAvsCTRLs_TRC1fdr_up_0 <- LAvsCTRLs_TRC1fdr_up_NA %>% replace(is.na(.), 0)
LAvsCTRLs_TRC1fdr_up_CS = LAvsCTRLs_TRC1fdr_up_0[,c("Pathways","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsCTRLs_TRC1fdr_up_CS$LEU  = "LA"
LAvsCTRLs_TRC1fdr_up_CS$CT  = "TRC1"



LAvsCTRLs_TRC2fdr_up <- read.xlsx("LAvsCTRLs_TRC2fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LAvsCTRLs_TRC2fdr_up$Pathways = str_replace(LAvsCTRLs_TRC2fdr_up$Term, "\\s*\\([^\\)]+\\)", "")
LAvsCTRLs_TRC2fdr_up_sub = subset(LAvsCTRLs_TRC2fdr_up, Pathways %in% barplot)
LAvsCTRLs_TRC2fdr_up_non = subset(barplott, !(Pathways %in% LAvsCTRLs_TRC2fdr_up$Pathways))
LAvsCTRLs_TRC2fdr_up_NA = rbindlist(list(LAvsCTRLs_TRC2fdr_up_sub, LAvsCTRLs_TRC2fdr_up_non), fill = TRUE)
LAvsCTRLs_TRC2fdr_up_NA$Log10Adj.P.value = -log10(LAvsCTRLs_TRC2fdr_up_NA$Adjusted.P.value)
LAvsCTRLs_TRC2fdr_up_0 <- LAvsCTRLs_TRC2fdr_up_NA %>% replace(is.na(.), 0)
LAvsCTRLs_TRC2fdr_up_CS = LAvsCTRLs_TRC2fdr_up_0[,c("Pathways","Combined.Score","Log10Adj.P.value", "Overlap")]
LAvsCTRLs_TRC2fdr_up_CS$LEU  = "LA"
LAvsCTRLs_TRC2fdr_up_CS$CT  = "TRC2"




bound_LA_vs_CTRL = rbind(LAvsCTRLs_FDCfdr_up_CS,
                         LAvsCTRLs_MRCfdr_up_CS, 
                         LAvsCTRLs_PRCfdr_up_CS, 
                         LAvsCTRLs_TRC1fdr_up_CS, 
                         LAvsCTRLs_TRC2fdr_up_CS, 
                         LAvsCTRLs_TNC1fdr_up_CS, 
                         LAvsCTRLs_TNC2fdr_up_CS)



nrow(bound_LA_vs_CTRL)


p = row.names(bound_LA_vs_CTRL)

nrow(unique(bound_LA_vs_CTRL))

bound_LA_vs_CTRL = as.data.frame(bound_LA_vs_CTRL)

bound_LA_vs_CTRL$gene.ratio = 0

for (i in 1:105) {
  print(i) 
  bound_LA_vs_CTRL[i,]$gene.ratio = eval(parse(text=bound_LA_vs_CTRL[i,]$Overlap)) }

bound_LA_vs_CTRL$Cond = paste(bound_LA_vs_CTRL$LEU, bound_LA_vs_CTRL$CT, sep = "_")

head(bound_LA_vs_CTRL)

head(bound_LA_vs_CTRL)

col.range = c(0,50)


#bound_LA_vs_CTRL[grep("cell.matrix adhesion", bound_LA_vs_CTRL$Pathways), ]
#cell.matrix adhesion



for (i in unique(bound_LA_vs_CTRL$Cond)) {
    i = "LA_FDC"
    print(i)  
    dati = bound_LA_vs_CTRL[grep(i, bound_LA_vs_CTRL$Cond), ]
    dati$Pathways <- factor(dati$Pathway, level = barplot)
    pup <-ggplot(dati, aes(x=Pathways, y=gene.ratio, fill=Log10Adj.P.value)) +
        geom_bar(stat="identity")+
        scale_fill_gradient(low = "blue", high = "red", limits=col.range) + 
        #scale_x_discrete(limits = dati$Pathway)+
        scale_y_continuous(limits=0:1)+
        labs(title = i) +
        theme_classic()+
        coord_flip()

  
  pdf(paste(i,"_vs_CTRL_UP_CellTypes.pdf", sep = "_"), 10, 5)
  pup
  dev.off() 
  }

getwd()


getwd()
LA_TRC1_vs_CTRL_UP_CellTypes.pdf

head(bound_LA_vs_CTRL)

bound_LA_vs_CTRL = unique(bound_LA_vs_CTRL)

head(bound_LA_vs_CTRL)
final = acast(unique(bound_LA_vs_CTRL), Term ~ Cond, value.var = 'Log10Adj.P.value')
ciaooo <- final %>% replace(is.na(.), 0)
ciaooo = as.data.frame(ciaooo)
ciaooone <- mutate_all(ciaooo, function(x) as.numeric(as.character(x)))
df = ciaooone[rowSums(ciaooone[])>0,]
head(df)


annotation_column = as.data.frame(unique(bound_LA_vs_CTRL[,c('Cond','LEU', 'CT')]))
rownames(annotation_column) <- annotation_column$Cond
head(annotation_column)


#unique(df_row[2]$PATHWAY)

category_order = c("ECM", "WNT", "Growth factor beta", "Chemotaxis", "T cell", "Proliferation", "TNF", "IFN", "Hypoxia", "Cytokine", "Migration", "Lipopolysaccharide")


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
          "CTRL" = '#FFB266'),
  CT = c('FDC' = '#F564E3', 
         'MRC' = '#F8766D', 
         'PRC' = '#00BFC4', 
         'TRC1' = '#619CFF', 
         'TRC2' = '#BF86BF',
         'TNC1' = '#B79F00', 
         'TNC2' = '#00BA38'))





annotation_column = annotation_column[,c("LEU","CT")]
rownames(annotation_column)
df_row <- data.frame(barplott)
rownames(df_row) = df_row$Term

#reorder df columns as in annotation_column
df <- df[,match(rownames(annotation_column), colnames(df))]


#reorder by PATHWAY the rownames
df_row = df_row[order(df_row$PATHWAY),]
#or by custom order
category_order = c("ECM", "WNT", "Growth factor beta", "Proliferation","Chemotaxis",  "TNF", "IFN", "Hypoxia", "Cytokine", "Migration", "Lipopolysaccharide")

head(df_row)
df_row <- df_row %>% 
  arrange(factor(PATHWAY, levels = category_order))

nrow(df_row[2])

nrow(df)


head(ann_colors)
df_row = subset(df_row, Term %in% rownames(df))
nrow(df_row)

df_1 <- df[match(df_row$Term, rownames(df)), ]  

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
minH = -1; maxH=1
myb = seq(minH, maxH, by = 0.01)
#crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))

#remove the description in brackets (GO:0071346)
rownames(df_1) = str_replace(rownames(df_1), "\\s*\\([^\\)]+\\)", " ")

rownames(df_row[2]) = str_replace(rownames(df_row[2]), "\\s*\\([^\\)]+\\)", " ")

pdf("LA_vs_CTRL.pdf")
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


##############################################
#############################################



setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Brendolan/BrendolanA_1190_scRNAseq/9_bioinfo/DEG_1190_merged_ctrls_EYFP+/DEG_and_MarkerGenes_TRC_divided_in_2_clusters/EnrichR_bis")

LBvsCTRLs_FDCfdr_up <- read.xlsx("LBvsCTRLs_FDCfdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_FDCfdr_up_sub = subset(LBvsCTRLs_FDCfdr_up, Term %in% barplott)
LBvsCTRLs_FDCfdr_up_non = subset(barplott, !(Term %in% LBvsCTRLs_FDCfdr_up$Term))
LBvsCTRLs_FDCfdr_up_NA = rbindlist(list(LBvsCTRLs_FDCfdr_up_sub, LBvsCTRLs_FDCfdr_up_non), fill = TRUE)
LBvsCTRLs_FDCfdr_up_NA$Log10Adj.P.value = -log10(LBvsCTRLs_FDCfdr_up_NA$Adjusted.P.value)
LBvsCTRLs_FDCfdr_up_0 <- LBvsCTRLs_FDCfdr_up_NA %>% replace(is.na(.), 0)
LBvsCTRLs_FDCfdr_up_CS = LBvsCTRLs_FDCfdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_FDCfdr_up_CS$LEU  = "LB"
LBvsCTRLs_FDCfdr_up_CS$CT  = "FDC"
length(unique(barplott))
nrow(unique(LBvsCTRLs_FDCfdr_up_sub))






LBvsCTRLs_PRCfdr_up <- read.xlsx("LBvsCTRLs_PRCfdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_PRCfdr_up_sub = subset(LBvsCTRLs_PRCfdr_up, Term %in% barplott)
LBvsCTRLs_PRCfdr_up_non = subset(barplott, !(Term %in% LBvsCTRLs_PRCfdr_up$Term))
LBvsCTRLs_PRCfdr_up_NA = rbindlist(list(LBvsCTRLs_PRCfdr_up_sub, LBvsCTRLs_PRCfdr_up_non), fill = TRUE)
LBvsCTRLs_PRCfdr_up_NA$Log10Adj.P.value = -log10(LBvsCTRLs_PRCfdr_up_NA$Adjusted.P.value)
LBvsCTRLs_PRCfdr_up_0 <- LBvsCTRLs_PRCfdr_up_NA %>% replace(is.na(.), 0)
LBvsCTRLs_PRCfdr_up_CS = LBvsCTRLs_PRCfdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_PRCfdr_up_CS$LEU  = "LB"
LBvsCTRLs_PRCfdr_up_CS$CT  = "PRC"



LBvsCTRLs_MRCfdr_up <- read.xlsx("LBvsCTRLs_MRCfdr_up.xlsx", sheet = "GO_Biological_Process_2021") 
#LBvsCTRLs_MRCfdr_up_sub = subset(LBvsCTRLs_MRCfdr_up, Term %in% barplott)
#LBvsCTRLs_MRCfdr_up_sub = barplott
#LBvsCTRLs_MRCfdr_up_non = barplott
####LBvsCTRLs_MRCfdr_up_non[setdiff(LBvsCTRLs_MRCfdr_up_non$Term, LBvsCTRLs_MRCfdr_up_sub$Term)] <- 0
LBvsCTRLs_MRCfdr_up_CS = barplott
LBvsCTRLs_MRCfdr_up_CS$Combined.Score = 0
LBvsCTRLs_MRCfdr_up_CS$Log10Adj.P.value = 0
LBvsCTRLs_MRCfdr_up_CS$Overlap = 0
LBvsCTRLs_MRCfdr_up_CS$LEU  = "LB"
LBvsCTRLs_MRCfdr_up_CS$CT  = "MRC"
LBvsCTRLs_MRCfdr_up_CS = subset(LBvsCTRLs_MRCfdr_up_CS, select = -c(2) )




LBvsCTRLs_TNC1fdr_up <- read.xlsx("LBvsCTRLs_TNC1fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TNC1fdr_up_sub = subset(LBvsCTRLs_TNC1fdr_up, Term %in% barplott)
LBvsCTRLs_TNC1fdr_up_non = subset(barplott, !(Term %in% LBvsCTRLs_TNC1fdr_up$Term))
LBvsCTRLs_TNC1fdr_up_NA = rbindlist(list(LBvsCTRLs_TNC1fdr_up_sub, LBvsCTRLs_TNC1fdr_up_non), fill = TRUE)
LBvsCTRLs_TNC1fdr_up_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TNC1fdr_up_NA$Adjusted.P.value)
LBvsCTRLs_TNC1fdr_up_0 <- LBvsCTRLs_TNC1fdr_up_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TNC1fdr_up_CS = LBvsCTRLs_TNC1fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TNC1fdr_up_CS$LEU  = "LB"
LBvsCTRLs_TNC1fdr_up_CS$CT  = "TNC1"



LBvsCTRLs_TNC2fdr_up <- read.xlsx("LBvsCTRLs_TNC2fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TNC2fdr_up_sub = subset(LBvsCTRLs_TNC2fdr_up, Term %in% barplott)
LBvsCTRLs_TNC2fdr_up_non = subset(barplott, !(Term %in% LBvsCTRLs_TNC2fdr_up$Term))
LBvsCTRLs_TNC2fdr_up_NA = rbindlist(list(LBvsCTRLs_TNC2fdr_up_sub, LBvsCTRLs_TNC2fdr_up_non), fill = TRUE)
LBvsCTRLs_TNC2fdr_up_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TNC2fdr_up_NA$Adjusted.P.value)
LBvsCTRLs_TNC2fdr_up_0 <- LBvsCTRLs_TNC2fdr_up_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TNC2fdr_up_CS = LBvsCTRLs_TNC2fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TNC2fdr_up_CS$LEU  = "LB"
LBvsCTRLs_TNC2fdr_up_CS$CT  = "TNC2"



LBvsCTRLs_TRC1fdr_up <- read.xlsx("LBvsCTRLs_TRC1fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TRC1fdr_up_sub = subset(LBvsCTRLs_TRC1fdr_up, Term %in% barplott)
LBvsCTRLs_TRC1fdr_up_non = subset(barplott, !(Term %in% LBvsCTRLs_TRC1fdr_up$Term))
LBvsCTRLs_TRC1fdr_up_NA = rbindlist(list(LBvsCTRLs_TRC1fdr_up_sub, LBvsCTRLs_TRC1fdr_up_non), fill = TRUE)
LBvsCTRLs_TRC1fdr_up_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TRC1fdr_up_NA$Adjusted.P.value)
LBvsCTRLs_TRC1fdr_up_0 <- LBvsCTRLs_TRC1fdr_up_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TRC1fdr_up_CS = LBvsCTRLs_TRC1fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TRC1fdr_up_CS$LEU  = "LB"
LBvsCTRLs_TRC1fdr_up_CS$CT  = "TRC1"



LBvsCTRLs_TRC2fdr_up <- read.xlsx("LBvsCTRLs_TRC2fdr_up.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TRC2fdr_up_sub = subset(LBvsCTRLs_TRC2fdr_up, Term %in% barplott)
LBvsCTRLs_TRC2fdr_up_non = subset(barplott, !(Term %in% LBvsCTRLs_TRC2fdr_up$Term))
LBvsCTRLs_TRC2fdr_up_NA = rbindlist(list(LBvsCTRLs_TRC2fdr_up_sub, LBvsCTRLs_TRC2fdr_up_non), fill = TRUE)
LBvsCTRLs_TRC2fdr_up_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TRC2fdr_up_NA$Adjusted.P.value)
LBvsCTRLs_TRC2fdr_up_0 <- LBvsCTRLs_TRC2fdr_up_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TRC2fdr_up_CS = LBvsCTRLs_TRC2fdr_up_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TRC2fdr_up_CS$LEU  = "LB"
LBvsCTRLs_TRC2fdr_up_CS$CT  = "TRC2"



LBvsCTRLs_FDCfdr_down <- read.xlsx("LBvsCTRLs_FDCfdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_FDCfdr_down_CS = barplott
LBvsCTRLs_FDCfdr_down_CS$Combined.Score = 0
LBvsCTRLs_FDCfdr_down_CS$Log10Adj.P.value = 0
LBvsCTRLs_FDCfdr_down_CS$Overlap = 0
LBvsCTRLs_FDCfdr_down_CS$LEU  = "CTRL"
LBvsCTRLs_FDCfdr_down_CS$CT  = "FDC"
LBvsCTRLs_FDCfdr_down_CS = subset(LBvsCTRLs_FDCfdr_down_CS, select = -c(2) )


LBvsCTRLs_PRCfdr_down <- read.xlsx("LBvsCTRLs_PRCfdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_PRCfdr_down_sub = subset(LBvsCTRLs_PRCfdr_down, Term %in% barplott)
LBvsCTRLs_PRCfdr_down_non = subset(barplott, !(Term %in% LBvsCTRLs_PRCfdr_down$Term))
LBvsCTRLs_PRCfdr_down_NA = rbindlist(list(LBvsCTRLs_PRCfdr_down_sub, LBvsCTRLs_PRCfdr_down_non), fill = TRUE)
LBvsCTRLs_PRCfdr_down_NA$Log10Adj.P.value = -log10(LBvsCTRLs_PRCfdr_down_NA$Adjusted.P.value)
LBvsCTRLs_PRCfdr_down_0 <- LBvsCTRLs_PRCfdr_down_NA %>% replace(is.na(.), 0)
LBvsCTRLs_PRCfdr_down_CS = LBvsCTRLs_PRCfdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_PRCfdr_down_CS$LEU  = "CTRL"
LBvsCTRLs_PRCfdr_down_CS$CT  = "PRC"



LBvsCTRLs_MRCfdr_down <- read.xlsx("LBvsCTRLs_MRCfdr_down.xlsx", sheet = "GO_Biological_Process_2021") # empty DF
LBvsCTRLs_MRCfdr_down_sub = subset(LBvsCTRLs_MRCfdr_down, Term %in% barplott)
LBvsCTRLs_MRCfdr_down_non = subset(barplott, !(Term %in% LBvsCTRLs_MRCfdr_down$Term))
LBvsCTRLs_MRCfdr_down_NA = rbindlist(list(LBvsCTRLs_MRCfdr_down_sub, LBvsCTRLs_MRCfdr_down_non), fill = TRUE)
LBvsCTRLs_MRCfdr_down_NA$Log10Adj.P.value = -log10(LBvsCTRLs_MRCfdr_down_NA$Adjusted.P.value)
LBvsCTRLs_MRCfdr_down_0 <- LBvsCTRLs_MRCfdr_down_NA %>% replace(is.na(.), 0)
LBvsCTRLs_MRCfdr_down_CS = LBvsCTRLs_MRCfdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_MRCfdr_down_CS$LEU  = "CTRL"
LBvsCTRLs_MRCfdr_down_CS$CT  = "MRC"




LBvsCTRLs_TNC1fdr_down <- read.xlsx("LBvsCTRLs_TNC1fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TNC1fdr_down_sub = subset(LBvsCTRLs_TNC1fdr_down, Term %in% barplott)
LBvsCTRLs_TNC1fdr_down_non = subset(barplott, !(Term %in% LBvsCTRLs_TNC1fdr_down$Term))
LBvsCTRLs_TNC1fdr_down_NA = rbindlist(list(LBvsCTRLs_TNC1fdr_down_sub, LBvsCTRLs_TNC1fdr_down_non), fill = TRUE)
LBvsCTRLs_TNC1fdr_down_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TNC1fdr_down_NA$Adjusted.P.value)
LBvsCTRLs_TNC1fdr_down_0 <- LBvsCTRLs_TNC1fdr_down_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TNC1fdr_down_CS = LBvsCTRLs_TNC1fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TNC1fdr_down_CS$LEU  = "CTRL"
LBvsCTRLs_TNC1fdr_down_CS$CT  = "TNC1"



LBvsCTRLs_TNC2fdr_down <- read.xlsx("LBvsCTRLs_TNC2fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TNC2fdr_down_sub = subset(LBvsCTRLs_TNC2fdr_down, Term %in% barplott)
LBvsCTRLs_TNC2fdr_down_non = subset(barplott, !(Term %in% LBvsCTRLs_TNC2fdr_down$Term))
LBvsCTRLs_TNC2fdr_down_NA = rbindlist(list(LBvsCTRLs_TNC2fdr_down_sub, LBvsCTRLs_TNC2fdr_down_non), fill = TRUE)
LBvsCTRLs_TNC2fdr_down_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TNC2fdr_down_NA$Adjusted.P.value)
LBvsCTRLs_TNC2fdr_down_0 <- LBvsCTRLs_TNC2fdr_down_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TNC2fdr_down_CS = LBvsCTRLs_TNC2fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TNC2fdr_down_CS$LEU  = "CTRL"
LBvsCTRLs_TNC2fdr_down_CS$CT  = "TNC2"



LBvsCTRLs_TRC1fdr_down <- read.xlsx("LBvsCTRLs_TRC1fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TRC1fdr_down_sub = subset(LBvsCTRLs_TRC1fdr_down, Term %in% barplott)
LBvsCTRLs_TRC1fdr_down_non = subset(barplott, !(Term %in% LBvsCTRLs_TRC1fdr_down$Term))
LBvsCTRLs_TRC1fdr_down_NA = rbindlist(list(LBvsCTRLs_TRC1fdr_down_sub, LBvsCTRLs_TRC1fdr_down_non), fill = TRUE)
LBvsCTRLs_TRC1fdr_down_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TRC1fdr_down_NA$Adjusted.P.value)
LBvsCTRLs_TRC1fdr_down_0 <- LBvsCTRLs_TRC1fdr_down_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TRC1fdr_down_CS = LBvsCTRLs_TRC1fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TRC1fdr_down_CS$LEU  = "CTRL"
LBvsCTRLs_TRC1fdr_down_CS$CT  = "TRC1"



LBvsCTRLs_TRC2fdr_down <- read.xlsx("LBvsCTRLs_TRC2fdr_down.xlsx", sheet = "GO_Biological_Process_2021")
LBvsCTRLs_TRC2fdr_down_sub = subset(LBvsCTRLs_TRC2fdr_down, Term %in% barplott)
LBvsCTRLs_TRC2fdr_down_non = subset(barplott, !(Term %in% LBvsCTRLs_TRC2fdr_down$Term))
LBvsCTRLs_TRC2fdr_down_NA = rbindlist(list(LBvsCTRLs_TRC2fdr_down_sub, LBvsCTRLs_TRC2fdr_down_non), fill = TRUE)
LBvsCTRLs_TRC2fdr_down_NA$Log10Adj.P.value = -log10(LBvsCTRLs_TRC2fdr_down_NA$Adjusted.P.value)
LBvsCTRLs_TRC2fdr_down_0 <- LBvsCTRLs_TRC2fdr_down_NA %>% replace(is.na(.), 0)
LBvsCTRLs_TRC2fdr_down_CS = LBvsCTRLs_TRC2fdr_down_0[,c("Term","Combined.Score","Log10Adj.P.value", "Overlap")]
LBvsCTRLs_TRC2fdr_down_CS$LEU  = "CTRL"
LBvsCTRLs_TRC2fdr_down_CS$CT  = "TRC2"










bound_LB_vs_CTRL = rbind(LBvsCTRLs_FDCfdr_up_CS, LBvsCTRLs_FDCfdr_down_CS,
                         LBvsCTRLs_MRCfdr_up_CS, LBvsCTRLs_MRCfdr_down_CS,
                         LBvsCTRLs_PRCfdr_up_CS, LBvsCTRLs_PRCfdr_down_CS,
                         LBvsCTRLs_TRC1fdr_up_CS, LBvsCTRLs_TRC1fdr_down_CS, 
                         LBvsCTRLs_TRC2fdr_up_CS, LBvsCTRLs_TRC2fdr_down_CS,
                         LBvsCTRLs_TNC1fdr_up_CS, LBvsCTRLs_TNC1fdr_down_CS,
                         LBvsCTRLs_TNC2fdr_up_CS, LBvsCTRLs_TNC2fdr_down_CS)


bound_LB_vs_CTRL$Cond = paste(bound_LB_vs_CTRL$LEU, bound_LB_vs_CTRL$CT, sep = "_")


LEUB = bound_LB_vs_CTRL[grep("LB", bound_LB_vs_CTRL$LEU), ]
CTRL = bound_LB_vs_CTRL[grep("CTRL", bound_LB_vs_CTRL$LEU), ]

bound_LB_vs_CTRL = unique(bound_LB_vs_CTRL)

head(bound_LB_vs_CTRL)
final = acast(unique(bound_LB_vs_CTRL), Term ~ Cond, value.var = 'Log10Adj.P.value')
ciaooo <- final %>% replace(is.na(.), 0)
ciaooo = as.data.frame(ciaooo)
ciaooone <- mutate_all(ciaooo, function(x) as.numeric(as.character(x)))
df = ciaooone[rowSums(ciaooone[])>0,]
head(df)


annotation_column = as.data.frame(unique(bound_LA_vs_CTRL[,c('Cond','LEU', 'CT')]))
rownames(annotation_column) <- annotation_column$Cond
head(annotation_column)


#unique(df_row[2]$PATHWAY)

category_order = c("ECM", "WNT", "Growth factor beta", "Chemotaxis", "T cell", "Proliferation", "TNF", "IFN", "Hypoxia", "Cytokine", "Migration", "Lipopolysaccharide")


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
  LEU = c("LB" = '#3333FF', 
          "CTRL" = '#FFB266'),
  CT = c('FDC' = '#F564E3', 
         'MRC' = '#F8766D', 
         'PRC' = '#00BFC4', 
         'TRC1' = '#619CFF', 
         'TRC2' = '#BF86BF',
         'TNC1' = '#B79F00', 
         'TNC2' = '#00BA38'))





annotation_column = annotation_column[,c("LEU","CT")]
rownames(annotation_column)
df_row <- data.frame(barplott)
rownames(df_row) = df_row$Term

#reorder df columns as in annotation_column
df <- df[,match(rownames(annotation_column), colnames(df))]


#reorder by PATHWAY the rownames
df_row = df_row[order(df_row$PATHWAY),]
#or by custom order
category_order = c("ECM", "WNT", "Growth factor beta", "Proliferation","Chemotaxis",  "TNF", "IFN", "Hypoxia", "Cytokine", "Migration", "Lipopolysaccharide")

head(df_row)
df_row <- df_row %>% 
  arrange(factor(PATHWAY, levels = category_order))

nrow(df_row[2])

nrow(df)


head(ann_colors)
df_row = subset(df_row, Term %in% rownames(df))
nrow(df_row)

df_1 <- df[match(df_row$Term, rownames(df)), ]  

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
minH = -1; maxH=1
myb = seq(minH, maxH, by = 0.01)
#crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))

#remove the description in brackets (GO:0071346)
rownames(df_1) = str_replace(rownames(df_1), "\\s*\\([^\\)]+\\)", " ")

rownames(df_row[2]) = str_replace(rownames(df_row[2]), "\\s*\\([^\\)]+\\)", " ")

pdf("LB_vs_CTRL.pdf")
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


####################################








#############################################
#############################################










LEUA = as.data.frame(LEUA)
library(reshape)
LA = LEUA[,c("Term","Combined.Score", "CT")]
LA = unique(LA)

final = acast(unique(LA), Term ~ CT, value.var = 'Combined.Score')
ciaooo <- final %>% replace(is.na(.), 0)
ciaooo = as.data.frame(ciaooo)
ciaooone <- mutate_all(ciaooo, function(x) as.numeric(as.character(x)))


#neg reg Wnt signaling pathway
negative regulation of canonical Wnt signaling pathway (GO:0090090)
negative regulation of Wnt signaling pathway (GO:0030178)
Wnt signaling pathway, planar cell polarity pathway (GO:0060071)

#pos regulation of Wnt
positive regulation of canonical Wnt signaling pathway (GO:0090263)
regulation of canonical Wnt signaling pathway (GO:0060828)
# one more

# interferon
regulation of type I interferon production (GO:0032479)
cellular response to type I interferon (GO:0071357)
type I interferon signaling pathway (GO:0060337)

# Tcell
T cell chemotaxis (GO:0010818)
T cell migration (GO:0072678)
# one more

# leukocyte
positive regulation of leukocyte chemotaxis (GO:0002690)
positive regulation of leukocyte migration (GO:0002687)

# transforming growth factor beta
##cellular response to growth factor stimulus (GO:0071363)
cellular response to transforming growth factor beta stimulus (GO:0071560)
transforming growth factor beta receptor signaling pathway (GO:0007179)
regulation of transforming growth factor beta receptor signaling pathway (GO:0017015)

# TNF
cellular response to tumor necrosis factor (GO:0071356)
negative regulation of tumor necrosis factor superfamily cytokine production (GO:1903556)
response to tumor necrosis factor (GO:0034612)


# Hypoxia
regulation of transcription from RNA polymerase II promoter in response to hypoxia (GO:0061418)
cellular response to hypoxia (GO:0071456)
#

### regulation of establishment of planar polarity (GO:0090175) No

# ECM
extracellular matrix disassembly (GO:0022617)
extracellular matrix organization (GO:0030198)
positive regulation of chemotaxis (GO:0050921)
#extracellular matrix assembly (GO:0085029)

# Cytokine
response to cytokine (GO:0034097)
cellular response to cytokine stimulus (GO:0071345)
#

## regulation of mesenchymal cell proliferation (GO:0010464) No

cellular response to interleukin.6 (GO:0071354)
interleukin.6.mediated signaling pathway (GO:0070102)
interleukin.1.mediated signaling pathway (GO:0070498)
cellular response to interleukin.1 (GO:0071347)
interleukin.12.mediated signaling pathway (GO:0035722)
cellular response to interleukin.12 (GO:0071349)
regulation of interleukin.8 secretion (GO:2000482)
regulation of interleukin.2 production (GO:0032663)


pdf('Cicciiiiii.pdf', 5, 5)
pheatmap::pheatmap(as.matrix(df), 
                   #cellheight = 25,
                   #cellwidth = 25,
                   #cluster_rows = T,  
                   #cluster_cols = T, 
                   scale = 'row', 
                   fontsize = 5)
dev.off()

nonblood = myData[!grepl("0071222", row.names(df)), ]

dff = df[!grepl("0071222", row.names(df)), ]
dff = dff[!grepl("GO:0048247", row.names(dff)), ]


pdf('Mandiiiiiiii.pdf', 5, 5)
pheatmap::pheatmap(as.matrix(dff), 
                   #cellheight = 25,
                   #cellwidth = 25,
                   #cluster_rows = T,  
                   #cluster_cols = T, 
                   scale = 'row', 
                   fontsize = 5)
dev.off()





ciaooo[ciaooo==0] <- NA
ciaooo


data2<-ciaooo[complete.cases(ciaooo),]

nrow(ciaooo)
nrow(data2)


pdf("vdjqvjsj.pdf", 20,10)
pheatmap(
  mat               = ciaooo,
  main              = "Default Heatmap",
  scale = "row"
)
dev.off()



head(ciaooo)

############

regulation of insulin-like growth factor receptor singling pathway


