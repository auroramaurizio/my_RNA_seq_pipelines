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



# Biomart does not work anymore but it was:

require("biomaRt")

convertMouseGeneList <- function(x){
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = human, attributesL = c("entrezgene_id"), martL = human, uniqueRows=T)
  return(genesV2)
}

converted=convertMouseGeneList(gene_list)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hypeR")


library("hypeR")
