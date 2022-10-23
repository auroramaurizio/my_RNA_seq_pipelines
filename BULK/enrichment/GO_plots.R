#######################

### Enrichment analysis 


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



#############################################################################
# enrichment with EnrichR example
#############################################################################


# full list of pathway source

databases <- listEnrichrDbs()

# selected databases to make the enrichment of (chosen from databases)

enrich.databases <- c("GO_Biological_Process_2021",
                      #"GO_Cellular_Component_2021",
                      "GO_Molecular_Function_2021",
                      "MSigDB_Hallmark_2020")



# Importing genes

my_genes = unique(as.character(read.delim("diff.genes_30.txt",h=F)$V2))

# Performing enrichment

my_enrichr = enrichR::enrichr(genes=my_genes, databases=enrich.databases)

# Saving results

openxlsx::write.xlsx(x=my_enrichr, file="diff.genes_30.enrichR.xlsx")


#############################################################################
# Plot the top N most significant pathways found with EnrichR with ggplot
#############################################################################


#################################
# Up regulated DEGs e.g. FC > 0 , padj < 0.05
#################################

# Setup parameters

fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))

N=10 # Plot the top 10 most significant pathways with ggplot

enrichR.file = "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Brendolan/BrendolanA_1190_scRNAseq/9_bioinfo/Bulk_Simona/Task5/bulk/enrichment/enrichR_bis/LEUB_vs_LEUA_fdr_up_.xlsx"


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

up = head(pathways.dataframe[1:N,], 10)
up$Condition = "UP"




###############################################
# Down regulated DEGs e.g. FC < 0 , padj < 0.05
###############################################


enrichR.file = "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Brendolan/BrendolanA_1190_scRNAseq/9_bioinfo/Bulk_Simona/Task5/bulk/enrichment/enrichR_bis/LEUB_vs_LEUA_fdr_down_.xlsx"


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


dw = head(pathways.dataframe[1:N,], 10)
dw$Condition = "DW"



##############################################
# Combine the top UP and DW regulated pathways
##############################################

dati = rbind(up,dw)


dati$Log10Adj.P.value = -log10(dati$p.value.adj)

#x <- parse(text=dati$Overlap)
#for(i in seq_along(x))
#  dati$gene.ratio[i] <- eval(x[[i]])

library(stringr)

#remove the description in brackets (GO:0071346)
dati$Pathway = str_replace(dati$Pathway, "\\s*\\([^\\)]+\\)", " ")

p <- ggplot(dati, aes(Condition,Pathway)) +
  # + geom_point() 
  geom_point(mapping = aes_string(size = "gene.ratio", color = "Log10Adj.P.value")) + 
  scale_color_gradientn(colours = c("blue", "red"), name = "-Log10(p.adj.value)")+
  #scale_size_continuous(range=c(1,7), name="Gene ratio")
  scale_size_binned(limits = c(0.1, 0.7), n.breaks = 5)+
  theme(axis.text.x = element_text( color='black', face="bold.italic", size=14)) +
  labs(title = "LEUB vs LEUA") 

getwd()
pdf("pathway_bulk.pdf",5, 5)
p
dev.off()

write.xlsx(dati, "pathway_LEUB_vs_LEUA_bulk_toplot.xlsx")

###############################################################################################################
###############################################################################################################
