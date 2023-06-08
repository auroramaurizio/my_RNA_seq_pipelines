

library(ggplot2)
library(openxlsx)
library(Seurat)
library(EWCE)



KO_vs_WT = read.xlsx(
  "/Users/maurizio.aurora/DGE_results_Cherry.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)


head(KO_vs_WT)
KO_vs_WT$MGI.symbol = rownames(KO_vs_WT)


UP = KO_vs_WT[KO_vs_WT$log2FoldChange > 0.58 & KO_vs_WT$pvalue < 0.05, ]

DW = KO_vs_WT[KO_vs_WT$log2FoldChange < -0.58 & KO_vs_WT$pvalue < 0.05, ] 

custom_crush_KO_vs_WT = rbind(UP, DW )


head(custom_crush_KO_vs_WT)

custom_tt_results_crush_KO_vs_WT = ewce_expression_data(sct_data=ctd,tt=custom_crush_KO_vs_WT,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange", 
                                                        reps = 10000)


head(custom_tt_results_crush_KO_vs_WT)
write.xlsx(custom_tt_results_crush_KO_vs_WT$joint_results,
           file = "custom_KO_vs_WT_barplot.xlsx", 
           row.names = T,
           asTable = T)


custom_tt_results_crush_KO_vs_WT$joint_results$SD_from_mean = custom_tt_results_crush_KO_vs_WT$joint_results$sd_from_mean
custom_tt_results_crush_KO_vs_WT$joint_results$SD_from_mean [which(custom_tt_results_crush_KO_vs_WT$joint_results$SD_from_mean < 0)] = 0
custom_tt_results_crush_KO_vs_WT$joint_results$pvalue = NA
custom_tt_results_crush_KO_vs_WT$joint_results$pvalue[custom_tt_results_crush_KO_vs_WT$joint_results$p<0.05]<-'*'


head(custom_tt_results_crush_KO_vs_WT)


S3_test <- ggplot(custom_tt_results_crush_KO_vs_WT$joint_results)+ 
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

pdf("custom_crush_KO_vs_WT_dotplot_minSDzero_pvalue0.05.pdf",8,10)
S3_test
dev.off()

16/2776
