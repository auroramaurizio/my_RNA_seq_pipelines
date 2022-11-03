#####################################################################################
#####################################################################################

setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/_/_A_1190_scRNAseq/9_bioinfo/COND_A_COND_B_merged_ctrls_EYFP+Marker_genes/figures/STROMA/CAMs/CAMs_final/enrichments_cherry_picked/dotplot_CAM_cherry_picked")

enrichR.file = "myCAMs_COND_B_vs_COND_Afdr_up.xlsx"


s = c("GO_Biological_Process_2021")
enrichR.table = data.frame()

for (dat in s) {
  
  Table <- openxlsx::read.xlsx(xlsxFile=enrichR.file, sheet=dat, startRow=1, colNames=T, rowNames=T,
                               
                               detectDates=F, skipEmptyRows=T, skipEmptyCols=T, na.strings="NA", fillMergedCells=F)
  
  enrichR.table = rbind(enrichR.table, Table)
  
}

# Considering only significant pathways

#p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < 0.05,])

p = row.names(enrichR.table)

if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
  
  pathways.dataframe = data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value)
  
  # Formatting the dataframe for the plot
  
  pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
  
  pathways.dataframe$Pathway.num = as.factor(dim(pathways.dataframe)[1]:1) # as.factor is somehow necessary to write the pathway names
  
}

myCAMs = pathways.dataframe
myCAMs$Condition = "myCAMs" 

nrow(myCAMs)

enrichR.file = "iCAMs_COND_B_vs_COND_Afdr_up.xlsx"

getwd()

s = c("GO_Biological_Process_2021")
enrichR.table = data.frame()

for (dat in s) {
  
  Table <- openxlsx::read.xlsx(xlsxFile=enrichR.file, sheet=dat, startRow=1, colNames=T, rowNames=T,
                               
                               detectDates=F, skipEmptyRows=T, skipEmptyCols=T, na.strings="NA", fillMergedCells=F)
  
  enrichR.table = rbind(enrichR.table, Table)
  
}

# Considering only significant pathways

p = row.names(enrichR.table)



if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
  
  pathways.dataframe = data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value)
  
  # Formatting the dataframe for the plot
  
  pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
  
  pathways.dataframe$Pathway.num = as.factor(dim(pathways.dataframe)[1]:1) # as.factor is somehow necessary to write the pathway names
  
}


iCAMs = pathways.dataframe
iCAMs$Condition = "iCAMs"

nrow(iCAMs)
head(iCAMs)



enrichR.file = "apCAMs_COND_B_vs_COND_Afdr_up.xlsx"

getwd()

s = c("GO_Biological_Process_2021")
enrichR.table = data.frame()

for (dat in s) {
  
  Table <- openxlsx::read.xlsx(xlsxFile=enrichR.file, sheet=dat, startRow=1, colNames=T, rowNames=T,
                               
                               detectDates=F, skipEmptyRows=T, skipEmptyCols=T, na.strings="NA", fillMergedCells=F)
  
  enrichR.table = rbind(enrichR.table, Table)
  
}

# Considering only significant pathways

p = row.names(enrichR.table)

fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))


if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
  
  pathways.dataframe = data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value)
  
  # Formatting the dataframe for the plot
  
  pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
  
  pathways.dataframe$Pathway.num = as.factor(dim(pathways.dataframe)[1]:1) # as.factor is somehow necessary to write the pathway names
  
}


apCAMs = pathways.dataframe
apCAMs$Condition = "apCAMs"

nrow(apCAMs)
nrow(iCAMs)
nrow(myCAMs)

merged = rbind(iCAMs,apCAMs)
merged = rbind(merged,myCAMs)
dati = unique(merged)

nrow(dati)

dati$Log10Adj.P.value = -log10(dati$p.value.adj)

colnames(dati)
library(stringr)
dati$Pathways = str_replace(dati$Pathway, "\\s*\\([^\\)]+\\)", "")
head(dati$Pathways)





prova = c("regulation of microtubule nucleation",
          "positive regulation of tumor necrosis factor-mediated signaling pathway",
          "negative regulation of inclusion body assembly",
          "cytoplasmic translation",
          "protein targeting to ER",
          "rRNA processing",
          "peptide biosynthetic process",
          "cotranslational protein targeting to membrane",
          "protein targeting to ER",
          "SRP dependent protein targeting to ER",
          "extracellular matrix organization",
          "smooth muscle contraction",
          "integrin−mediated signaling pathway",
          "cellular protein metabolic process", 
          "SRP-dependent cotranslational protein targeting to membrane",
          "neutrophil activation involved in immune response",
          "neutrophil mediated immunity",
          "neutrophil degranulation")



#x <- parse(text=dati$Overlap)
#for(i in seq_along(x))
#  dati$gene.ratio[i] <- eval(x[[i]])



dati_sub = dati[dati$Pathways %in% prova,]

head(dati_sub)

#sizes <- factor(dati_sub$Pathways, levels = prova)

#dati$pathways = sizes



unique(dati_sub$Condition)

library(dplyr)
#dati_sub = dati_sub %>% arrange(factor(Pathways, levels = prova))

#dati_sub = unique(dati_sub)

#dati_sub$Path <- factor(dati_sub$Pathways, levels=unique(dati_sub$Pathways))
dati_sub$Path <- factor(dati_sub$Pathways, levels=unique(prova))



col.range = c(0,11)
p <- ggplot(dati_sub, aes(Condition,Path)) +
  # + geom_point() 
  geom_point(mapping = aes_string(size = "gene.ratio", color = "Log10Adj.P.value")) + 
  scale_color_gradientn(colours = c("blue", "red"), limits=col.range, name = "-Log10(p.adj.value)")+
  #scale_size_continuous(range=c(1,7), name="Gene ratio")
  scale_size_binned(limits = c(0.1, 0.5), n.breaks = 5)+
  theme(axis.text.x = element_text( color='black', face="bold.italic", size=8)) +
  labs(title = "UP COND_B") 


pdf("pathway_COND_B_vs_COND_A_stroma_sc_fdr_up_dotplot_CAMs_test_6.pdf",9, 7)
p
dev.off()

write.xlsx(dati, "pathway_COND_B_vs_COND_A_stroma_sc_fdr_up_dotplot_CAMs_test_6.xlsx")
