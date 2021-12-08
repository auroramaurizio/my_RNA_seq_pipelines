---
  title: "Bulk"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# set working diretory
prj = 'RNASeq'
PI = 'B'
#setwd(paste("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/",
#            PI,"/",prj,
#            "/7_bioinfo/",
#            sep=''))
#setwd("/Users/maurizio.aurora/Documents/combined_14")
setwd("/Users/maurizio.aurora/Documents/")


getwd()


# load libraries 
suppressMessages(library("edgeR"))
suppressMessages(library(data.table))
suppressMessages(library("DESeq2"))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library("RColorBrewer"))
suppressMessages(library(enrichR))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
#install.packages('VennDiagram')
suppressMessages(library('VennDiagram'))
#install.packages("venn")
suppressMessages(library(venn))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library("viridis"))
suppressMessages(library('pals'))
suppressMessages(library(RColorBrewer))
suppressMessages(library(wesanderson))
suppressMessages(library(patchwork))
suppressMessages(library(magick))
suppressMessages(library('philentropy'))
suppressMessages(library('IntClust'))
suppressMessages(library('pheatmap'))
suppressMessages(library(assertr))
suppressMessages(library("remotes"))
suppressMessages(library(GeneOverlap))
suppressMessages(library('stringr'))
suppressMessages(library(ggrepel))
#install.packages("multcompView", repos="http://R-Forge.R-project.org")
suppressMessages(library(multcompView))
suppressMessages(library(sva))

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

# path to the count file generated with pypette
#my_filecount = '/Users/maurizio.aurora/Documents/181020/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged/all.counts.gz'
#giulys_filecount='/Users/maurizio.aurora/Documents/BOOB/dataset/20191211/counts.gz'
#giulys_filecount='/Users/maurizio.aurora/Documents/041120/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged/all.counts.gz'


# path to the count file generated with pypette
my_filecount = '/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/raw_not_corrected_counts_tables_1133_1292/1292/all.counts.gz'
#giulys_filecount='/Users/maurizio.aurora/Documents/BOOB/dataset/20191211/counts.gz'
giulys_filecount='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/raw_not_corrected_counts_tables_1133_1292/1133/all.counts.gz'



#head(metadata_file_auri)
# Min number of replica in each group
#Nreplica= 4
Nreplica= 3

# import metadata

metadata = read.xlsx("/Users/maurizio.aurora/Documents/combine1/metadata.xlsx")
#remove outliers from metadata
metadata <- subset(metadata, sample!= c("72","DB_2_130"))


# import counts
annotation <- c('GeneID','Chr','Start','End','Strand','Length')
fCounts_auri <- read.delim(file=my_filecount, header=TRUE, check.names = F)
fCounts_giuly <- read.delim(file=giulys_filecount, header=TRUE, check.names = F)
#colnames(fCounts_giuly) =c("Geneid","101","104","72","74","77","80","83","92","95","98")

colnames(fCounts)

nrow(fCounts_auri) #55376 Geneid
nrow(fCounts_giuly) #53278 Geneid #55376

head(fCounts_giuly)
head(fCounts_auri)


intersected=intersect(fCounts_giuly$Geneid,fCounts_auri$Geneid)
length(intersected) #49579 #55376
head(intersected)

#fCounts <- merge(fCounts_auri, fCounts_giuly, by.x="Geneid", by.y="Geneid")
fCounts <- merge(fCounts_auri, fCounts_giuly, by=c("Geneid","Chr","Start","End","Length","Strand"))
nrow(fCounts) #49579 #55376
colnames(fCounts)

fCountsData_not_adj <- fCounts[
  , 
  -which(
    tolower(names(fCounts))
    %in% 
      tolower(annotation))]

fCountsAnnotation <- fCounts[
  , 
  which(
    tolower(names(fCounts))
    %in% 
      tolower(annotation))]

fCounts_bis= select(fCounts, "Geneid","Chr","Start","End","Length","Strand" )
colnames(fCounts_bis)

colnames(fCounts)
is.matrix(fCounts)
head(fCountsAnnotation)

geneidColname <- 'Geneid'
geneidIdx <- which(tolower(annotation) %in% tolower(geneidColname))
rownames(fCountsData_not_adj) <- fCounts[[geneidIdx]]

# Reordering counts matrix to have samples ordered as in metadata
# it assumes that the first column in metadata is sample name
fCountsData_not_adj <- fCountsData_not_adj[,match(metadata[,1], colnames(fCountsData_not_adj))] 




metadata_d = metadata[metadata$reporter == 'KDR_cherry',]
#metadata$Plxnd1 = as.factor(metadata$Plxnd1)
row.names(metadata_d) = metadata_d$sample
cov1=as.factor(metadata_d$nerve)
cov2=factor(metadata_d$Plxnd1, levels = c('WT', 'KO'))
cov3=as.factor(metadata_d$reporter)
exp_batch=as.factor(metadata_d$exp)
#covar_mat <- cbind(cov1, cov2, cov3)
covar_mat <- cbind(cov2, cov3)



#####################
## ONLY CHERRY SAMPLE
#####################
head(metadata)

metadata_d = metadata[metadata$reporter == 'KDR_cherry',]
#remove samples 72 and 130


getwd()

#metadata$Plxnd1 = as.factor(metadata$Plxnd1)
row.names(metadata_d) = metadata_d$sample
cov1=as.factor(metadata_d$nerve)
cov2=factor(metadata_d$Plxnd1, levels = c('WT', 'KO'))
cov3=as.factor(metadata_d$reporter)
exp_batch=as.factor(metadata_d$exp)
#covar_mat <- cbind(cov1, cov2, cov3)
covar_mat <- cbind(cov2, cov3)

fCountsData_adjusted_p = fCountsData_not_adj[,row.names(metadata_d)]

fCountsData_adjusted_adjusted <- ComBat_seq(as.matrix(fCountsData_adjusted_p), batch=exp_batch, group=cov2 ) #combined
fCountsData_adjusted_d=fCountsData_adjusted_adjusted

colnames(fCountsData_adjusted_d)


class(fCountsData_adjusted_p)


SAVE_variable <- list()
filename_xls <- paste('COUNTS_cherry',prj,'.xlsx', sep='')
variable2save_names <- c('all_counts', 'expGenes_counts','expGenes_LogCPM', 'expGenes_LogRPKM', 'expGenes_RPKM')

head(fCountsAnnotation)

# all_counts 
y <- DGEList(counts=fCountsData_adjusted_d, genes = fCountsAnnotation)

SAVE_variable[[variable2save_names[1]]] <- as.data.frame(y$counts)

# expGENES_counts
keep <- rowSums(cpm(y)>1)>=Nreplica
table(keep)
yf <- y[keep,]
SAVE_variable[[variable2save_names[2]]] <- as.data.frame(yf$counts)

#CPM
SAVE_variable[[variable2save_names[3]]] <- as.data.frame(cpm(yf, log=T))

#RPKM log
SAVE_variable[[variable2save_names[4]]] <- as.data.frame(rpkm(yf, log=T, gene.length =yf$genes$Length))

#RPKM not log
SAVE_variable[[variable2save_names[5]]] <- as.data.frame(rpkm(yf, log=F, gene.length =yf$genes$Length))


SAVE_variable <- list()
prj = "logRpkm_all_genes_cherry"
filename_xls <- paste('COUNTS_',prj,'.xlsx', sep='')
variable2save_names <- c('all_genes_log_rpkm')


#RPKM log all genes
SAVE_variable[[variable2save_names[1]]] <- as.data.frame(rpkm(y, log=T, gene.length =y$genes$Length))




write.xlsx(SAVE_variable,
           file = filename_xls, 
           row.names = T,
           asTable = T, 
           sheetName =variable2save_names)

getwd()


#####################
## ONLY CHERRY SAMPLE
#####################


metadata_d = metadata[metadata$reporter == 'KDR_cherry',]
#metadata$Plxnd1 = as.factor(metadata$Plxnd1)
row.names(metadata_d) = metadata_d$sample
cov1=as.factor(metadata_d$nerve)
cov2=factor(metadata_d$Plxnd1, levels = c('WT', 'KO'))
cov3=as.factor(metadata_d$reporter)
exp_batch=as.factor(metadata_d$exp)
#covar_mat <- cbind(cov1, cov2, cov3)
covar_mat <- cbind(cov2, cov3)

fCountsData_adjusted_p = fCountsData_not_adj[,row.names(metadata_d)]

fCountsData_adjusted_adjusted <- ComBat_seq(fCountsData_adjusted_p, batch=exp_batch, group=cov2 ) #combined
fCountsData_adjusted_d=fCountsData_adjusted_adjusted

#fCountsData_adjusted_d=fCountsData_adjusted[ , c("DB_1_127", "DB_2_130", "DB_5_148", "DB_6_151")]  


y <- DGEList(counts=fCountsData_adjusted_d, genes = fCountsAnnotation)
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$Length)

provaA=rowSums(cpm(y)>1)>=Nreplica
head(provaA)
keep <- rowSums(cpm(y)>1)>=Nreplica
length(provaA)
yf <- y[keep,]
nrow(yf)
N=500
vary <- apply(fCountsRPKM[keep,],1,var)
vary_s <- sort(vary, decreasing = T)
TOP_N <- names(vary_s[1:N])
yTOP <-  y[TOP_N,]
fCountsRPKMTOP <- fCountsRPKM[TOP_N,]
#PCA parameters
pcx = 1
pcy = 2
centering = TRUE
scaling = TRUE
# PCA
pca = prcomp(t(fCountsRPKMTOP), center=centering, scale=scaling)
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
score = as.data.frame(pca$x)
# plot paramters
xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
cum = var[pcx]+var[pcy]
names = rownames(pca$x)



score$exp = metadata_d$exp
score$reporter = metadata_d$reporter
score$nerve = metadata_d$nerve
score$Plxnd1 = metadata_d$Plxnd1
score$sampleID = metadata_d$sample

pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy], 
                        color=Plxnd1, shape = exp))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy], 
                                    color=Plxnd1, shape=exp, label = str_sub(x <- sampleID,-5,-1)),
                   size = 5,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
                   segment.color = 'grey50') +
  #geom_point(aes(size= IL7))+
  geom_point(size= 7)+
  labs(x=xlab, y=ylab, title=paste("PC",pcx," vs PC",pcy," scoreplot",sep="")) + 
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 18),
        legend.position="right", 
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_color_manual(values = brewer.pal(n = 4, name = "Spectral"))  

options(repr.plot.width=13, repr.plot.height=7)
pdf('PCA_top500rpkm_Cherry.pdf', width = 14, height = 10)
pca
dev.off()

annotation_column <- metadata_d[,2:(dim(metadata_d)[2])]


score$reporter = metadata_d$reporter
score$nerve = metadata_d$nerve
score$Plxnd1 = metadata_d$Plxnd1
score$sampleID = metadata_d$sample


row.names(annotation_column) <- metadata_d[,1]
#annotation_column = annotation_column[,c('condition','bay','IL7','cycle','Donor')]
options(repr.plot.width=12, repr.plot.height=10)

#colors <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255)
#colors = magma(255)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         #annotation_colors = ann_colors, 
                         cluster_rows = T, 
                         cluster_cols = T,                       
                         #main = 'Heatmap: 500 most variable genes - RPKM',
                         show_rownames = F,
                         #cutree_rows = 3,
                         #cutree_cols = 5,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                         display_numbers = F, 
                         col=colors,
                         filename = 'Heatmap_500rpkm_Cherry.pdf',
                         width = 14, height = 11)


##########

#metadata

f = "Plxnd1"
seqc_pvalue = 0.01
dgeResults = list()

metadata_d[metadata_d$Plxnd1 %in% comparison[[1]],]
fCountsData_adjusted_s = fCountsData_adjusted_d[,row.names(metadata_d)]
colnames(fCountsData_adjusted_d)
row.names(metadata_d) = metadata_d$sample


comparison = list(KO_vs_WT = c("KO","WT"))



for (i in names(comparison)) {
  
  metadata_s = metadata_d[metadata_d$Plxnd1 %in% comparison[[i]],]
  fCountsData_adjusted_s = fCountsData_adjusted_d[,row.names(metadata_s)]
  dds <- DESeqDataSetFromMatrix(
    countData = fCountsData_adjusted_s,
    colData  = metadata_s,
    #design   = as.formula('~ exp + Plxnd1'))
    design   = as.formula('~Plxnd1'))
  filter <- rowSums(cpm(counts(dds)) >= 1) >= Nreplica
  table(filter)
  ddsFiltered <- dds[filter,]
  dga <- DESeq(
    object = ddsFiltered,
    test = "Wald",
    fitType = "parametric",
    betaPrior = FALSE,
    minReplicatesForReplace = Inf)
  
  alpha = 0.05
  print(paste(comparison[[i]][1],"_vs_",comparison[[i]][2],sep=''))
  dgeResults.tmp <- results(dga,
                            contrast             = c(f,comparison[[i]][1],comparison[[i]][2]),
                            cooksCutoff          = Inf,
                            independentFiltering = TRUE,
                            alpha                = alpha,
                            pAdjustMethod        = "BH")
  summary(dgeResults.tmp)
  head(dgeResults.tmp$pvalue)
  dgeResults[[i]] <- dgeResults.tmp[order(dgeResults.tmp$pvalue, decreasing = F),]  
  
  #PCA
  vsd <- vst(dga, blind=FALSE)
  main.factor = "Plxnd1"
  
  
  
  pcaData <- plotPCA(vsd, intgroup=c(main.factor),returnData=TRUE)
  pcaData$exp=metadata_s$exp
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  #ggplot(pcaData, aes(PC1, PC2, color=Plxnd1, shape=exp))
  PCA = ggplot(pcaData, aes(PC1, PC2, color=Plxnd1, shape=exp)) +
    geom_label_repel(data= pcaData, aes(PC1, PC2, color=Plxnd1, label = str_sub(name,-5,-1)),
                     size = 6,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
                     segment.color = 'grey50') +
    geom_point(size=6) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    ggtitle(paste("PCA",i)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    scale_color_manual(values = c('salmon','grey'))
  pdf(paste('pca_',i,'.pdf',sep=''),width=8, height=6)
  print(PCA)
  dev.off()
  
  seqcUP = row.names(dgeResults[[i]])[dgeResults[[i]]$pvalue <= seqc_pvalue & 
                                        !is.na(dgeResults[[i]]$padj)&
                                        dgeResults[[i]]$log2FoldChange > 1]
  
  if (length(seqcUP) > 0) {
    # print heatmap
    annotation_column <- metadata_s[,2:(dim(metadata_s)[2])]
    row.names(annotation_column) <- metadata_s[,1]
    options(repr.plot.width=12, repr.plot.height=10)
    crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
    colors = crp(255)
    head(cpm(counts(dga)))
    length(seqcUP)
    HP <- pheatmap::pheatmap(cpm(counts(dga))[seqcUP,],
                             scale = 'row',
                             annotation_col = annotation_column,
                             annotation_colors = ann_colors, 
                             cluster_rows = T, 
                             cluster_cols = T,                       
                             #main = 'Heatmap: 500 most variable genes - RPKM',
                             show_rownames = F,
                             #cutree_rows = 3,
                             cutree_cols = 2,
                             fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                             display_numbers = F, 
                             col=colors,
                             filename = paste('Heatmap_seqcUP',i,'.pdf'),
                             width = 10, height = 11 )
  }
  
}


f = 'DGE_results_Cherry'
dir.create(f, showWarnings=TRUE, recursive=TRUE)

lapply(
  names(dgeResults),
  function(x) write.table(
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname),
    file.path(f, paste(x, ".tsv", sep="")),
    append=F,
    row.names=F,
    col.names=T,
    quote=F,
    sep="\t"))


dgeResults_table = list()    
dgeResults_table = lapply(
  names(dgeResults),
  function(x) 
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))

names(dgeResults_table) = names(dgeResults)


write.xlsx(dgeResults_table,
           file = 'DGE_results_Cherry.xlsx', 
           row.names = F,
           asTable = T, 
           startRow = 1,
           sheetName = str_sub(names(dgeResults),1,31)) 


##############



#plots
n.label = 20
FDR = T
pvalue = 0.01
for (Plxnd1 in names(dgeResults)) {
  results = as.data.frame(dgeResults[[Plxnd1]])
  results$DE = 'unm'
  if (!FDR) {
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'up'}
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < .1,]$DE = 'down' }
  } else {
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'up'}
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'down' }        
  }
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'SEQCup'}
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]$DE = 'SEQCdown' }
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'FDRup'}
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'FDRdown' }        
  
  results$DE <- factor(x = results$DE, levels = c("unm", "FDRdown","FDRup", 'SEQCdown','SEQCup'))
  mycolors = c('grey','dodgerblue4','darkred','dodgerblue2','coral'); names(mycolors) = levels(results$DE)
  results$DE2 = 'unm'; results[results$DE!='unm',]$DE2 = 'mod'
  results$DE2 <- factor(x = results$DE2, levels = c("mod","unm"))
  mysize = c(3,2); names(mysize) = levels(results$DE2)
  myalpha = c(1,0.2); names(mysize) = unique(results$DE2)
  
  # label N genes
  N = min(n.label, length(rownames(results[results$DE == 'FDRup',])))
  up_label = rownames(results[results$DE == 'FDRup',])[1:N]
  N = min(n.label, length(rownames(results[results$DE == 'FDRdown',])))
  down_label = rownames(results[results$DE == 'FDRdown',])[1:N]
  
  MAplot = ggplot(results) +
    geom_point(aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    xlim(c(0,1.e5)) +
    scale_x_continuous(trans='log10') +
    ggtitle(paste("MAPlot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    #scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +  
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Mean expression", y = "log2 fold change")  
  
  print(MAplot)
  pdf(paste(f,'/','MAplot_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(MAplot)
  dev.off()
  
  # Vulcano plot 
  VP = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    #geom_label_repel(data= results[up_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
    geom_label_repel(data= results[down_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #geom_label_repel(data= results["Plxnd1",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results[down_label,]), size = 3,
                     #label = row.names(results[up_label,]), size = 3,
                     #label = row.names(results["Plxnd1",]), size = 3,
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")   
  
  #pdf(paste('VulcanoPlot_Up',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  pdf(paste('VulcanoPlot_Dw_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(VP)
  dev.off()  
  
  
  
  
  # Vulcano plot 
  VP0 = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[up_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #geom_label_repel(data= results[down_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #geom_label_repel(data= results["Plxnd1",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #label = row.names(results[down_label,]), size = 3,
                     label = row.names(results[up_label,]), size = 3,
                     #label = row.names(results["Plxnd1",]), size = 3,
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")   
  
  print(VP0)
  #pdf(paste(f,'/','VulcanoPlot_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Plnd1',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  pdf(paste('VulcanoPlot_Up_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Dw',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(VP0)
  dev.off()
  
  
  
  
  # Vulcano plot 
  VP1 = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    #geom_label_repel(data= results[up_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
    #geom_label_repel(data= results[down_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
    geom_label_repel(data= results["Plxnd1",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #label = row.names(results[down_label,]), size = 3,
                     #label = row.names(results[up_label,]), size = 3,
                     label = row.names(results["Plxnd1",]), size = 3,
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")   
  
  print(VP1)
  #pdf(paste(f,'/','VulcanoPlot_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  pdf(paste('VulcanoPlot_Plnd1_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Up',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Dw',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(VP1)
  dev.off()
  
  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;   
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  
  pdf(paste(f,'/','MA_VP_',Plxnd1,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  dev.off()
  
  #A1 <- image_read_pdf(paste(f,'/','MA_VP_',Plxnd1,'.pdf',sep=''), density = 70)
  #image_write(A1, path = paste(f,'/','MA_VP_',Plxnd1,'.tiff',sep=''), format = "tiff")
  
}



####### FDR ########

fdrUP = list()

alpha = 0.05
fdrUP = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange > 0])
names(fdrUP)= names(dgeResults)       

fdrDW = list()
fdrDW = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange < 0])
names(fdrDW)= names(dgeResults)

####### SEQC ########
seqcUP = list()

pvalue = 0.01
seqcUP = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange > 1])
names(seqcUP)= names(dgeResults)               

seqcDW = list()
seqcDW = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange < -1])
names(seqcDW)= names(dgeResults)

print('FDRup');lengths(fdrUP); print('FDRdw');lengths(fdrDW)
print('SEQCup');lengths(seqcUP); print('SEQCdw');lengths(seqcDW)

print('FDR_tot');lengths(fdrUP) + lengths(fdrDW)
print('SEQC_tot'); lengths(seqcUP) + lengths(seqcDW)


######################

nrow(yf)



###### ENRICHR ########
databases <- listEnrichrDbs()

dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
# -------------------------
# enrichment Parameters
# -------------------------
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016",
                      "Jensen_TISSUES",
                      "Jensen_COMPARTMENTS",
                      "Jensen_DISEASES")
# alpha used in DGE
padj.cutoff = alpha; 

# -------------------------
# Perform Enrichment
# -------------------------
enrichr.list <- list()
for (i in 1:length(dgeResults)){
  print(names(dgeResults)[i])
  .res <- dgeResults[[i]]
  up.genes   <- fdrUP[[i]]
  down.genes <- fdrDW[[i]]
  both.genes <- c(up.genes, down.genes)
  write.table(up.genes, paste('./enrichR/FDRup_',names(dgeResults)[i],
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  write.table(down.genes, paste('./enrichR/FDRdw_',names(dgeResults)[i],
                                '.txt', sep =''), 
              quote = F, row.names = F, col.names = F)
  write.table(both.genes, paste('./enrichR/FDRboth_',names(dgeResults)[i],
                                '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  
  enrichr.list[[i]] <- lapply(list(up.genes,down.genes,both.genes),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down","fdr_both")
}
names(enrichr.list) <- names(dgeResults)

# -----------------------------
# Write excels files
# -----------------------------


for (i in 1:length(dgeResults)){
  for (j in c("fdr_up","fdr_down","fdr_both")){
    filename = paste(
      file.path('enrichR', 
                names(dgeResults)[[i]]),
      j,
      ".xlsx",
      sep="_")
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
  }
}

write.xlsx(x = metadata, file= "metadata.xlsx")




###### enrichment #######

require("biomaRt")

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genesV2)
}


require("biomaRt")

# GSEA list 
outdir = 'GSEA/'
dir.create(outdir)
for (c in names(dgeResults)) {
  l_ranked = data.frame(GeneID= row.names(dgeResults[[c]]), LogFC = dgeResults[[c]][,'log2FoldChange'])
  l_ranked = l_ranked[order(l_ranked$LogFC, decreasing = T),]
  musGenes <-  l_ranked$GeneID
  converted=convertMouseGeneList(musGenes)
  result = merge(l_ranked, converted, by.x = "GeneID", by.y = "MGI.symbol")
  result1 <- result[c(3,2)]
  colnames(result1)<- c("GeneID","LogFC")
  result1 = result1[order(result1$LogFC, decreasing = T),]
  write.table(result1, file = paste(outdir, c,'_ranked_list.rnk', sep =''), 
              quote = F, row.names= F, col.names = F, sep ='\t')
  
}



#########################

### Jaccard distances ###
# -------------------------
# Jaccard distances
# -------------------------
dir = paste("./enrichR/JaccardPlots/", sep = '')
dir.create(dir)
#'WikiPathways_2016'
#database = c('KEGG_2016','Reactome_2016','GO_Biological_Process_2018')
p.value.thr = 0.01

breaksList = seq(0, 1, by = 0.001)

#cutree_rows_values = c(5,12,3,10,5,1)



lf = list.files(paste("./enrichR/", sep=''), pattern=glob2rx("*.xlsx"))
print(lf)
for (file in lf) {
  enrichR.file = paste("./enrichR/",file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  Pathways.Table = data.frame()
  for (dat in s[1:length(s)]) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  length(pathways)
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';')) 
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)
  
  if (length(gene.all) !=0) {
    # MAtrix
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways 
    colnames(M) = gene.all 
    
    for (pat in pathways) {
      for (gene in gene.all) {
        if (gene %in% gene.list[[pat]]) {
          M[pat,gene] <- 1 
        }            
      }    
    }
    
    if (length(pathways) >1) {
      # Jaccard dist
      Jacard.Matrix <- distance(M, method = "jaccard")
      if (length(pathways)==2) {
        Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
        Jacard.Matrix = Jacard.Matrix_new
      }
      
      row.names(Jacard.Matrix) <- pathways
      colnames(Jacard.Matrix) <- pathways
      w=30; h=80; fs = 5; #cutree_rows_N = 10
      myb = seq(0,1,by = 0.01)
      myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(myb))
      pheatmap(Jacard.Matrix,
               border_color = 'darkgrey',
               color = myc, 
               breaks = myb,
               cluster_rows = TRUE,
               cluster_cols = TRUE, 
               cutree_rows = cutree_rows_N,
               show_colnames = FALSE,
               main = paste(file,'- Jaccard distance heatmap'),
               fontsize = 8,
               fontsize_row = fs,
               filename = paste(dir,file,'cherry_JaccardDist.pdf', sep=''),
               width=w, height=h)
    }
  }
}




######
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genesV2)
}


musGenes <-  l_ranked$GeneID
converted=convertMouseGeneList(musGenes)
result = merge(l_ranked, converted, by.x = "GeneID", by.y = "MGI.symbol")
result1 <- result[c(3,2)]
colnames(result1)<- c("GeneID","LogFC")
result1 = result1[order(result1$LogFC, decreasing = T),]
write.table(result1, file = paste(outdir, c,'_ranked_list.rnk', sep =''), 
            quote = F, row.names= F, col.names = F, sep ='\t')

VEGFR2 <- c("DUSP5","DUSP2","VWF","DUSP1","EGF","CUL3","RASAL1","DUSP8","DUSP16","DUSP6","RAP1B","DUSP10","RASGEF1A","PPP2R5E","RASA1","UBC","FGFR2","HBEGF")
converted=convertHumanGeneList(VEGFR2)
humanx <- unique(converted[, 2])

head(humanx)


VEGF_Biocarta <- c("ARNT","EIF1","EIF1AX","EIF2B1","EIF2B2","EIF2B3","EIF2B4","EIF2B5","EIF2S1","EIF2S2","EIF2S3","ELAVL1",
                   "FLT1","FLT4","HIF1A","HRAS","KDR","NOS3","PIK3CA","PIK3CG","PIK3R1","PLCG1","PRKCA","PRKCB","SHC1","VEGFA","VHL")
converted=convertHumanGeneList(VEGF_Biocarta)
humanx <- unique(converted[, 2])
humanx1 <-  humanx[humanx %in% row.names(dgeResults[[i]])]

TGFbetaHSWP366 <- c("TGIF1","MEF2C","JUN","ROCK1","FOS","DAB2","MYC","EP300","FOSB","SIK1","JUNB","MET","ATF3")					
converted=convertHumanGeneList(TGFbetaHSWP366)
humanx <- unique(converted[, 2])
#humanx1 <-  humanx[humanx %in% row.names(dgeResults[[i]])]

TGFbetaMM <- c("TGIF1","BMP4","JUN","EGF","FST","EP300","FOS")
converted=convertHumanGeneList(TGFbetaMM)
humanx <- unique(converted[, 2])


if (length(humanx) > 0) {
  # print heatmap
  annotation_column <- metadata_s[,2:(dim(metadata_s)[2])]
  row.names(annotation_column) <- metadata_s[,1]
  options(repr.plot.width=12, repr.plot.height=10)
  crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
  colors = crp(255)
  head(cpm(counts(dga)))
  length(seqcUP)
  HP <- pheatmap::pheatmap(cpm(counts(dga))[humanx,],
                           scale = 'row',
                           annotation_col = annotation_column,
                           annotation_colors = ann_colors, 
                           cluster_rows = T, 
                           cluster_cols = T,                       
                           main = 'TGFbetaMM',
                           show_rownames = T,
                           #cutree_rows = 3,
                           cutree_cols = 2,
                           fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                           display_numbers = F, 
                           col=colors,
                           filename = paste('TGFbetaMM',i,'.pdf'),
                           width = 10, height = 11 )
}


TNFalphaNFkB <- c("RPL30","PSMD12","MCM7","RIPK3","STAT1","CAV1","ACTL6A","PEBP1","CD3EAP","RPL8","PML","PSMB5","POLR1B","PSMD3","CSNK2B","MCM5","POLR2H")							
converted=convertHumanGeneList(TNFalphaNFkB)
humanx <- unique(converted[, 2])
humanx1 <-  humanx[humanx %in% row.names(dgeResults[[i]])]


if (length(humanx1) > 0) {
  # print heatmap
  annotation_column <- metadata_s[,2:(dim(metadata_s)[2])]
  row.names(annotation_column) <- metadata_s[,1]
  options(repr.plot.width=12, repr.plot.height=10)
  crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
  colors = crp(255)
  head(cpm(counts(dga)))
  length(seqcUP)
  HP <- pheatmap::pheatmap(cpm(counts(dga))[humanx1,],
                           scale = 'row',
                           annotation_col = annotation_column,
                           annotation_colors = ann_colors, 
                           cluster_rows = T, 
                           cluster_cols = T,                       
                           main = 'TNFalphaNFkB',
                           show_rownames = T,
                           #cutree_rows = 3,
                           cutree_cols = 2,
                           fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                           display_numbers = F, 
                           col=colors,
                           filename = paste('TNFalphaNFkB',i,'.pdf'),
                           width = 10, height = 11 )
}


#surfaceome
setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/")
f = 'Surfaceome'
dir.create(f, showWarnings=TRUE, recursive=TRUE)
setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/Surfaceome")

getwd()

S2_File_CSPA = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/literature/S2_File_CSPA.xlsx",
  sheet = 2,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)

head(S2_File_CSPA)

#result = merge(l_ranked, converted, by.x = "GeneID", by.y = "MGI.symbol")


crush_d14_vs_crush_d7 = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)


added_surface_CSPA <- merge(x=crush_d14_vs_crush_d7, y=S2_File_CSPA, by.x=c("Geneid"), by.y=c("ENTREZ.gene.symbol"),all.x= TRUE,no.dups)
nrow(added_surface_CSPA)
added_surface_CSPA = added_surface_CSPA[order(added_surface_CSPA$padj, decreasing = F),]

write.xlsx(added_surface_CSPA, 
           'crush_d14_vs_crush_d7_surfaceome_CSPA.xlsx', 
           asTable = T, row.names= F)


crush_d14_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 2,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)

#crush_d14_vs_intact

added_surface_CSPA <- merge(x=crush_d14_vs_intact, y=S2_File_CSPA, by.x=c("Geneid"), by.y=c("ENTREZ.gene.symbol"),all.x= TRUE,no.dups)
nrow(added_surface_CSPA)
added_surface_CSPA = added_surface_CSPA[order(added_surface_CSPA$padj, decreasing = F),]

write.xlsx(added_surface_CSPA, 
           'crush_d14_vs_intact_surfaceome_CSPA.xlsx', 
           asTable = T, row.names= F)



crush_d7_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 3,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)

#crush_d14_vs_intact

added_surface_CSPA <- merge(x=crush_d7_vs_intact, y=S2_File_CSPA, by.x=c("Geneid"), by.y=c("ENTREZ.gene.symbol"),all.x= TRUE,no.dups)
nrow(added_surface_CSPA)
added_surface_CSPA = added_surface_CSPA[order(added_surface_CSPA$padj, decreasing = F),]

write.xlsx(added_surface_CSPA, 
           'crush_d7_vs_intact_surfaceome_CSPA.xlsx', 
           asTable = T, row.names= F)


head(crush_d7_vs_intact)
head(added_surface_CSPA)

nrow(crush_d7_vs_intact)
nrow(added_surface_CSPA)


test = added_surface_CSPA %>% count(Geneid, sort = TRUE)
head(test)
added_surface_CSPA["H2-K1" %in% added_surface_CSPA$Geneid,]



#surfaceome
setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/only_KDRCherry_Combat_no_outliers/")
f = 'Surfaceome'
dir.create(f, showWarnings=TRUE, recursive=TRUE)
setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/only_KDRCherry_Combat_no_outliers/Surfaceome")


KO_vs_WT = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/only_KDRCherry_Combat_no_outliers/DGE_results_Cherry/DGE_results_Cherry.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)


added_surface_CSPA <- merge(x=KO_vs_WT, y=S2_File_CSPA, by.x=c("Geneid"), by.y=c("ENTREZ.gene.symbol"),all.x= TRUE,no.dups)
nrow(added_surface_CSPA)
added_surface_CSPA = added_surface_CSPA[order(added_surface_CSPA$padj, decreasing = F),]

write.xlsx(added_surface_CSPA, 
           'KO_vs_WT_surfaceome_CSPA.xlsx', 
           asTable = T, row.names= F)

head(added_surface_CSPA)
nrow(KO_vs_WT)
nrow(added_surface_CSPA)

########################
########################

#Surfaceome surfy

#Load samples
crush_d14_vs_crush_d7 = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)

crush_d14_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 2,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)


crush_d7_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 3,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)


KO_vs_WT = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/only_KDRCherry_Combat_no_outliers/DGE_results_Cherry/DGE_results_Cherry.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)


#Load Surfy



#Convert to mouse


require("biomaRt")

convertMouseGeneList <- function(x){
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

converted=convertMouseGeneList(rownames(count_raw))

#Intersect



########################
########################













head(surface.proteins_surfy)




head(added_surface_CSPA)

head(added_surface_CSPA)

head(added_surface)



print(table(added_surface_CSPA$Geneid))
test = added_surface_CSPA %>% count(Geneid, sort = TRUE)
head(test)

colnames(added_surface)
write.xlsx(added_surface_CSPA, 
           'Surfaceome_details_1.xlsx', 
           asTable = T, row.names= F)

#629



surface.proteins_surfy <- added_surface[added_surface$Surfaceome.Label %in% c('surface'), ]
nrow(surface.proteins_surfy)
surface.proteins_CSPA_yes = added_surface_CSPA[added_surface_CSPA$UniProt.cell.surface %in% c('yes'),] 
head(surface.proteins_CSPA_yes)

added_surface_CSPA_bus <- merge(x=surface.proteins_surfy, y=surface.proteins_CSPA_yes, by=c("Geneid"))

nrow(added_surface_CSPA_bus)
head(added_surface_CSPA_bus)


head(surface.proteins_CSPA_yes)

Surfy=surface.proteins_surfy$Geneid
CSPA=surface.proteins_CSPA_yes$Geneid


########################################








### import libraries 

#Configuration
setwd('/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE')
# biomart dictionary
#hg19
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                 host="grch37.ensembl.org", 
                 path="/biomart/martservice", 
                 dataset="hsapiens_gene_ensembl")
#hg38
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# Surfaceome table
surfaceome_table = '/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Alfano/958_surfaceome/table_S3_surfaceome.xlsx'
ST = read.xlsx(surfaceome_table, startRow = 2)
ST = ST[!is.na(ST$UniProt.name),]
ST = ST[!is.na(ST$GeneID),]
row.names(ST) <- ST$UniProt.name

fdrUP_filtered = intersect(row.names(gene_dictionary_NEW), fdrUP)
fdrDW_filtered = intersect(row.names(gene_dictionary_NEW), fdrDW)
entrez_fdrUP_filtered = gene_dictionary_NEW[fdrUP_filtered,'entrez']
geneid = unlist(strsplit(as.character(entrez_fdrUP_filtered),', '))
geneid = geneid[!is.na(geneid)]
geneid = geneid[geneid != 'NA']


head(ST$UniProt.gene)
head(ST$Ensembl.gene)

require("biomaRt")

humGenes <- c("FER1L5", "TMEM129", "ESYT2")


convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  print(genesV2)
  #humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  return(genesV2)
}

genes <- convertHumanGeneList(ST$UniProt.gene)
colnames(genes) = c("UniProt.gene","MGI.symbol")
#pheno1 <- merge(ST, genes, c("UniProt.gene"), all.x = TRUE)
pheno1 <- merge(ST, genes, c("UniProt.gene"))

priva =pheno1[pheno1$MGI.symbol %in% crush_d14_vs_crush_d7$Geneid,]

nrow(priva)

nrow(crush_d14_vs_crush_d7)
genes["Mbp" %in% genes$MGI.symbol,]

colnames(genes)
"Gdpd5" %in% surface.proteins$MGI.symbol,]
nrow(ST)
nrow(pheno1)
nrow(genes)
colnames(pheno1)









x = list(Surfy, CSPA)
length(x)

venn.plot <- venn.diagram(
  x = x,
  category.names = c("Surfy","CSPA"),
  main = 'Surface proteins d14 vs d7 ',
  filename = "crush_d14_vs_crush_d7_Surfy_vs_CSPA.png",
  scaled = TRUE,
  ext.text = TRUE,
  fill = c("#FC8D5980","#FFFFBF80"),
  alpha = 0.6,
  cex = 2,
  cat.cex = 1,
  cat.pos = c(-20,20),
  lty = 'blank',
  #cat.col = mycolors_c[1:2],
  main.cex = 2,
  sub.cex = 1)


mag24_GSE97239$study = 'GSE97239'
mag24_ArchS4$study = 'ArchS4'
mag24_GSE133624$study = 'GSE133624'
mag24_plot = rbind(mag24_GSE97239,mag24_ArchS4,mag24_GSE133624)

mag24_plot <- ggplot(mag24_plot, aes(x = ProteinID, y = Log2FC, fill = study)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, vjust = .5, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "Protein", y = "LogFC")

mag24_plot
pdf(paste('magic24_logFC.pdf',sep=''),width=13, height=8)
plot(mag24_plot)
dev.off()











Surfaceome_details.xlsx

#how many surface proteins in 1
#how many surface proteins in 2

colnames(added_surface_CSPA)


#How many in common?

nrow(added_surface_CSPA)

nrow(crush_d14_vs_crush_d7)


head(crush_d14_vs_crush_d7)

colnames(crush_d14_vs_crush_d7)
MGI.symbol


1296 mouse surfaceome proteins and their annotation.



head(surface.proteins,2)

sub.surface.proteins=surface.proteins[,c("UniProt.gene","UniProt.description","Surfaceome.Label","topology.source","Membranome.Almen.main-class","Membranome.Almen.sub-class","MGI.symbol")]
head(sub.surface.proteins)
added_surface <- merge(x=crush_d14_vs_crush_d7, y=sub.surface.proteins, by.x=c("Geneid"), by.y=c("MGI.symbol"),all.x= TRUE,no.dups)

added_surface = added_surface[order(added_surface$log2FoldChange, decreasing = T),]


length(unique(added_surface$Geneid))

colnames(added_surface)
write.xlsx(added_surface, 
           'Surfaceome_details.xlsx', 
           asTable = T, row.names= F)


colnames(added_surface)
test = added_surface %>% count(Geneid, sort = TRUE)
table(added_surface$MGI.symbol)


write.xlsx(proteins, 
           'proteins_details.xlsx', 
           asTable = T, row.names= F)


nrow(added_surface)

proteins = pheno1[pheno1$MGI.symbol %in% crush_d14_vs_crush_d7$Geneid,]










colnames(surface.proteins)
nrow(added_surface)
nrow(crush_d14_vs_crush_d7)
nrow(crush_d14_vs_crush_d7)

head(crush_d14_vs_crush_d7$Geneid)
nrow(crush_d14_vs_crush_d7)
nrow(pheno1)
nrow(proteins)

colnames(pheno1)
proteins = pheno1[pheno1$MGI.symbol %in% crush_d14_vs_crush_d7$Geneid,]
proteins = proteins[!is.na(proteins$Surfaceome.Label),]
print(table(proteins$Surfaceome.Label))
surface.proteins = proteins[proteins$Surfaceome.Label=='surface',] 

proteins = proteins[!is.na(proteins$Surfaceome.Label),]
new = rownames(surface.proteins)
nrow(surface.proteins)
head(surface.proteins$MGI.symbol)
#nonsurface    surface 
#1447       1254 
nrow(proteins)


surfaceome_DGEresults = surfaceome_DGEresults[order(surfaceome_DGEresults$log2FoldChange, decreasing = T),]
write.xlsx(surfaceome_DGEresults,
           paste('DGE_Surfaceome_',i,'_',m,'_crush_d14_vs_crush_d7.xlsx',sep=''), 
           asTable = T, row.names= T)
write.xlsx(surface.proteins[surfaceome_DGEresults$ProteinID,], 
           paste('Surfaceome_details',i,'_',m,'_crush_d14_vs_crush_d7.xlsx', sep =''), 
           asTable = T, row.names= F)





dgeResults = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/7_bioinfo/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
)






# biomart dictionary
#hg19
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                 host="grch37.ensembl.org", 
                 path="/biomart/martservice", 
                 dataset="hsapiens_gene_ensembl")
#hg38
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# Surfaceome table

#Surfaceome
  fdrUP <- row.names(dgeResults)[dgeResults$padj <= alpha & 
                                        !is.na(dgeResults$padj)&
                                        dgeResults$log2FoldChange > 0]
  
  fdrDW <- row.names(dgeResults)[dgeResults$padj <= alpha & 
                                        !is.na(dgeResults$padj) &
                                        dgeResults$log2FoldChange < 0]  
  
  gene_dictionary = getBM(attributes = c( 'ensembl_gene_id',"entrezgene_id", "hgnc_symbol"),
                          filter = 'hgnc_symbol',
                          values = c(fdrUP, fdrDW),
                          mart = ensembl)
  
  gene_dictionary_NEW = data.frame()
  
  for (ii in unique(gene_dictionary$hgnc_symbol)) {
    gd.entry = data.frame(hgnc_symbol = ii,
                          entrez = paste(gene_dictionary[gene_dictionary$hgnc_symbol == ii,]$entrezgene_id, 
                                         collapse=', '))
    gene_dictionary_NEW = rbind(gene_dictionary_NEW, gd.entry)
  }
  row.names(gene_dictionary_NEW)  <-  gene_dictionary_NEW$hgnc_symbol
  
  
  modulation = c('UP','DW')
  fdr = list(UP = fdrUP, DW = fdrDW)
  for (m in modulation) {
    print(m) 
    fdr_filtered = intersect(row.names(gene_dictionary_NEW), fdr[[m]])
    entrez_fdr_filtered = gene_dictionary_NEW[fdr_filtered,'entrez']
    geneid = unlist(strsplit(as.character(entrez_fdr_filtered),', '))
    geneid = geneid[!is.na(geneid)]
    geneid = geneid[geneid != 'NA']
    proteins = ST[ST$GeneID %in% geneid,]
    proteins = proteins[!is.na(proteins$Surfaceome.Label),]
    print(table(proteins$Surfaceome.Label))
    surface.proteins = proteins[proteins$Surfaceome.Label=='surface',] 
    Surface_dictionary = gene_dictionary[gene_dictionary$entrezgene_id %in% surface.proteins$GeneID,]
    hgnc_Surface = Surface_dictionary$hgnc_symbol
    Surface_dictionary$SP = ' '
    for (j in 1:dim(Surface_dictionary)[1]) {
      Surface_dictionary$SP[j] = row.names(surface.proteins[surface.proteins$GeneID == 
                                                              Surface_dictionary$entrezgene_id[j],]) 
    }
    surfaceome_DGEresults = as.data.frame(dgeResults[[i]][hgnc_Surface,])
    surfaceome_DGEresults$ProteinID = Surface_dictionary$SP
    surfaceome_DGEresults = surfaceome_DGEresults[order(surfaceome_DGEresults$log2FoldChange, decreasing = T),]
    write.xlsx(surfaceome_DGEresults,
               paste('DGE_Surfaceome_',i,'_',m,'_crush_d14_vs_crush_d7.xlsx',sep=''), 
               asTable = T, row.names= T)
    write.xlsx(surface.proteins[surfaceome_DGEresults$ProteinID,], 
               paste('Surfaceome_details',i,'_',m,'_crush_d14_vs_crush_d7.xlsx', sep =''), 
               asTable = T, row.names= F)        
  }
}




####################################



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 10
cols = gg_color_hue(n)

cols


