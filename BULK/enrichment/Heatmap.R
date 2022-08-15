
library(ggplot2)
library(edgeR)
library(grid)
library(gridExtra)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(openxlsx)
library(viridis)
library(stringr)
library(pheatmap)
library(ComplexHeatmap)



rpkm = read.xlsx("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Task5/bulk/filtered.deseq2.toptable_clean.ALL_contrast.mark_seqc.header_added.notlog.rpkm_copy.xlsx",
                 sheet = 1,
                 startRow = 1,
                 colNames = TRUE,
                 rowNames = TRUE,
)


rpkm = rpkm[c("sample1","sample2",
              "sample3","sample4")]



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


Cadherins = read.xlsx("/Users/maurizio.aurora/Desktop/Cadherins.xlsx")
converted_Cadherins = convert_human_to_mouse(Cadherins$Approved.symbol)
converted_Cadherins = unique(converted_Cadherins)

annotation_row = converted_Cadherins

rpkm_CAM = rpkm[rownames(rpkm) %in% annotation_row,]

counts_filtered_df <- rpkm_CAM[apply(rpkm_CAM, MARGIN = 1, FUN = function(x) sd(x) != 0),]

pdf('Cadherins_bulk_names_new.pdf', 5, 10)
pheatmap::pheatmap(as.matrix(rpkm_CAM), 
                   main = "Cadherins",
                   show_rownames = T, 
                   show_colnames = F,
                   border_color = NA,
                   #cellheight = 25,
                   #cellwidth = 25,
                   cluster_rows = T,  
                   cluster_cols = T, 
                   col=myc,
                   scale = 'row', 
                   fontsize = 5, 
                   #annotation_colors = ann_colors, 
                   annotation_col = Metadata)
dev.off()



