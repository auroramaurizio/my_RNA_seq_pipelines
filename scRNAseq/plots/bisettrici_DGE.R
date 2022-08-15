testM = readRDS(file = "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/9_bioinfo/Seurat_objects/EYFP+_comb_ctrl_MRC_TRC_FDC.Rds")

DefaultAssay(testM) = "RNA"

TRC1_LEUA_vs_TRC1_combined_ctrl_mergedMRC_TRC_FDC =read.xlsx("TRC1_LEUA_vs_TRC1_combined_ctrl_mergedMRC_TRC_FDC.xlsx", rowNames = T)


TRC1 <- subset(testM, idents = "TRC1")
Idents(TRC1) <- "stim"
TRC1 <- log1p(AverageExpression(TRC1, verbose = FALSE)$RNA)
TRC1 = as.data.frame(TRC1)
TRC1$gene <- rownames(TRC1)
head(TRC1)


pval = subset(TRC1_LEUA_vs_TRC1_combined_ctrl_mergedMRC_TRC_FDC, p_val_adj < 0.05)

up10 <- pval %>% top_n(n = 10, wt = avg_logFC)
dw10 <- pval %>% top_n(n = -10, wt = avg_logFC)

up = up10$gene
down = dw10$gene

cc = c(up,down)

pdf("bisett_LEUA_vs_ctrl_merged_TRC1.pdf", 5,5)
plot=ggplot(TRC1, aes(LEUA,combined_ctrl), label=genes) + geom_point() + ggtitle("TRC1") + theme(
  plot.title = element_text(color="black", size=14, face="bold.italic"))
LabelPoints(plot = plot, points = cc, color="red", repel = TRUE) 
dev.off()


