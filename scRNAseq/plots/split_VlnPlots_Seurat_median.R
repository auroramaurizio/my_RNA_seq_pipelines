suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))


St_c = readRDS("S_c.RDS")

S_c = subset(St_c, idents = c("myC"))


#customize Seurat "split.by" condition VlnPlots (in this case split.by = "stim") to show median as a line


pmc = VlnPlot(mc, "myC", assay = "RNA", group.by = "stim", cols = c("#619CFF","#00BA38","#F8766D"), pt.size = 0)+scale_y_continuous(limits = c(-0.5,1))+
  stat_summary(fun=median, geom="point", color="black", shape = 95,  size = 15) + theme(axis.line.x=element_blank(),
                                                                                        axis.title.x=element_blank(),
                                                                                        legend.position="none",
                                                                                        axis.text.x=element_blank(),
                                                                                        axis.ticks.x=element_blank())+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")+ ggtitle("")+
  labs(caption="myC") +
  theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))




apc = subset(S_c, idents = c("apC"))

papc = VlnPlot(apc, "myC", assay = "RNA", group.by = "stim", cols = c("#619CFF","#00BA38","#F8766D"), pt.size = 0)+ scale_y_continuous(limits = c(-0.5,1))+
  stat_summary(fun=median, geom="point", color="black", shape = 95,  size = 15) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                        axis.title.x=element_blank(),
                                                                                        legend.title = element_blank(),
                                                                                        axis.title.y=element_blank())+ ggtitle("")+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")+ labs(x = "")+
  labs(caption="apC") +
  theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))



ic = subset(S_c, idents = c("iC"))

pic = VlnPlot(ic, "myC", assay = "RNA", group.by = "stim", cols = c("#619CFF","#00BA38","#F8766D"), pt.size = 0)+ scale_y_continuous(limits = c(-0.5,1))+
  stat_summary(fun=median, geom="point", color="black", shape = 95,  size = 15)  + theme(legend.position = 'none')+ theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                                          axis.title.x=element_blank(),
                                                                                                                          legend.title = element_blank(),
                                                                                                                          axis.title.y=element_blank(),legend.position="none")+ ggtitle("myCAFs")+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")+ labs(x = "")+
  labs(caption="iC") +
  theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))





library(patchwork)

pdf("myC_Violin_split.pdf", 7, 5)
pmc | pic | papc
dev.off()

#######################################################


St_c = readRDS("S_c.RDS")

md = St_c@meta.data
## change its rownames according to your wish
md$cond = gsub("combined_ctrl", "Ctrl", md$stim)
md$cond = gsub("LEUA", "leuk-A", md$cond)
md$cond = gsub("LEUB", "leuk-B", md$cond)
## check that md has the desired format of rownames
## insert it back
St_c@meta.data = md

my_comparisons <- list( c("Ctrl", "leuk-A"), c("Ctrl", "leuk-B"), c("leuk-A", "leuk-B") )


VlnPlot(St_c, "myGene", assay = "RNA", group.by = "cond", cols = c("#F8766D","#00BA38","#619CFF"), pt.size = 0)+ 
  stat_summary(fun=median, geom="point", color="black", shape = 95,  size = 15)  +
  labs(caption="myGene") + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)+
  theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))





