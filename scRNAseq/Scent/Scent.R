# run SCENT injury

injury <- as.data.frame(as.matrix(sample2@assays$RNA@data))


pt_idents_sample2 <- as.data.frame(Idents(sample2))
colnames(pt_idents_sample2) = "CT"
pheno = pt_idents_sample2$CT

length(pheno)
gene_list = rownames(injury)



#prova = c("TdTomato","Xkr4","Gm19938","Rp1","Sox17","Mrpl15")

#output = data.frame("mouse"="topo", "human"="uomo")

#pp = convert_mouse_to_human(prova)

human_mouse_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

convert_mouse_to_human <- function(gene_list){
  for(gene in gene_list){
    class_key = (human_mouse_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']] 
    if(!identical(class_key, integer(0)) ){
      human_genes = (human_mouse_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"EntrezGene.ID"]
      for(human_gene in human_genes){
        output <- rbind(output, list(gene, human_gene))
      }    }
    else{
      output <- rbind(output, list(gene, 'NULL'))
    }
  }
  
  return (output)
}

converted_injury = convert_mouse_to_human(gene_list)

nrow(converted_injury)
converted_injury_ = converted_injury[grep("topo", converted_injury$mouse), ]

converted_injury_ = subset(converted_injury, mouse!="topo")
converted_injury__ = subset(converted_injury_, human!="NULL")
head(converted_injury__)
converted_injury___ <- converted_injury__[!duplicated(converted_injury__$mouse),]
rownames(converted_injury___) = converted_injury___$mouse

merged_injury=merge(converted_injury___,injury,by="row.names")
merged_injury <- merged_injury[-c(1,2,3)]
colnames(merged_injury)

nrow(merged_injury)


#Log-transforming the library-size normalized scRNA-Seq data matrix is important to stabilize the variance of highly expressed genes. Because SCENT uses ratios of expression values it is important to avoid 0 values after log-transformation. For this reason, we use a pseudocount of +1.1 instead of +1. This particular offset ensures that the minimum value after log-transformation is greater than zero. However, for CCAT, 0 values in the data matrix are allowed, so below we also define a separate log-transformed data matrix using the usual pseudocount of +1. This is convenient to retain the sparse nature of the data matrix if it happens to be of the dgC-class (see `Matrix` package).


merged_injury <- log2(merged_injury+1.1)
range(merged_injury)

merged_injury_bis <- log2(merged_injury+1)
range(merged_injury_bis)

###########

integ.l <- DoIntegPPI(exp.m = merged_injury, ppiA.m = net13Jun12.m)
str(integ.l)

sr.o <- CompSRana(integ.l, local = FALSE, mc.cores = 4)

rownames(sr.o)

nrow(integ.l)

x_axis_labels=function(labels,every_nth=1,...) {
  axis(side=1,at=seq_along(labels),labels=F)
  text(x=(seq_along(labels))[seq_len(every_nth)==1],
       y=par("usr")[3]-0.075*(par("usr")[4]-par("usr")[3]),
       labels=labels[seq_len(every_nth)==1],xpd=TRUE,...)
}
# axis() draws the axis with ticks at positions specified by at.  Again, we don't plot the labels yet.
# text() plots the labels at positions given by x and y.
# We estimate the y-positions from the values of the y-axis (using par("usr")),

getwd()


dataframeSR = data.frame("SR" = sr.o$SR, "CT" = pheno )
nrow(pheno)
nrow(dataframeSR)
#saveRDS(dataframeSR, "dataframeSR_intact.RDS")
#saveRDS(sr.o$SR, "sr.oSR_injury.RDS")
saveRDS(dataframeSR, "dataframeSR_injury.RDS")
length(pheno)
col = c('VENOUS_PLVAP+' = '#F8766D',
        'VENOUS_PLVAP-' = 'brown',
        'ARTERIAL' = '#00BFC4', 
        'TIP_1' = '#CD9600',
        'TIP_2' = 'pink',
        'BARR_END_CAP' = '#7CAE00',
        'CAPILLARY_PLVAP-' = '#00BE67',
        'CAPILLARY_PLVAP+' = 'grey',
        'TIP_3' = 'purple')

pdf("injury_SCENT_SR_3TIP.pdf", 15, 5)
p = ggplot(dataframeSR , aes(x=CT, y=SR, fill=CT)) +
  geom_boxplot(color="black") +
  labs(title="Injury",x="CellType", y = "Entropy")

p+scale_fill_manual(values=col)
dev.off()




pot.o <- InferPotencyStates(potest.v=dataframeSR$SR, pheno.v = pheno)
pot.o$distr



###########################
###########################
###########################
