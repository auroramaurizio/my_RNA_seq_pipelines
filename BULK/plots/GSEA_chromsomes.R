library(openxlsx)
setwd("/Users/maurizio.aurora/Downloads")
dgeResults <- read.xlsx("Supplemental_Table_S1.xlsx", sheet = "LeuA vs Ctrl")
head(dgeResults)
rownames(dgeResults) <- dgeResults$GeneID
forGSEA <- as.data.frame(dgeResults$logFC)
colnames(forGSEA) <- "logFC"
head(forGSEA)
rownames(forGSEA) = rownames(dgeResults)
head(forGSEA)
forGSEA$gene <- rownames(forGSEA)
head(forGSEA)
forGSEA <- forGSEA[order(forGSEA$logFC), ]
head(forGSEA)

dfList<-list(forGSEA)
names(dfList) = c("LeuA_vs_Ctrl")

#######
dir.create('GSEA/', showWarnings=FALSE, recursive=TRUE)

# -------------------------
# enrichment Parameters
# -------------------------
# databases to make the enrichment of
GSEA.databases <- c("MSigDB")


library(clusterProfiler)
library(msigdbr)
library(clusterProfiler)
library(msigdbr)
library(openxlsx)
library(enrichplot)

msigdb_df <- read.gmt("/Users/maurizio.aurora/Downloads/m1.all.v2024.1.Mm.symbols.gmt")

# Sample dataframe
df <- msigdb_df
head(df)
df$term <- gsub("^(chr[0-9]+).*", "\\1", df$term)
# Remove suffixes to match canonical names
df$term <- gsub("^(chr[0-9]+).*", "\\1", df$term)
df$term <- gsub("chrX[[:alnum:]]+", "chrX", df$term)
df$term <- gsub("chrY[[:alnum:]]+", "chrY", df$term)
# Print the updated dataframe
unique(df$term)

msigdb_df <- df


#unique(msigdb_df$gs_name)
colnames(msigdb_df) <- c("gs_name","gene_symbol")
#msigdb_df <- msigdbr(species = "Mus musculus", category = "M2")
#msigdb_df[msigdb_df$gene_symbol == "Prdm15",]

tail(msigdb_df)
#unique(msigdb_df$gs_name)
#MSigDB = msigdbr_df[grepl("MSigDB", msigdbr_df$gs_name),]
geneset_msigdb_df = msigdb_df[,c("gs_name","gene_symbol")]

# Create output directory
dir.create('GSEA/', showWarnings = FALSE, recursive = TRUE)

# Enrichment parameters
GSEA.databases <- c("geneset_msigdb_df")  # Placeholder for database name

# Example TERM2GENE setup (replace 'C1' with the correct gene set data frame)
# TERM2GENE must be a data frame with two columns: term ID and gene name
# Example:
# C1 <- data.frame(term = c("Term1", "Term2"), gene = c("Gene1", "Gene2"))

# Loop through data frames in the list
for (i in names(dfList)) {
  df_obj <- dfList[[i]]
  
  # Ensure 'gene' column exists; replace this if needed
  if (!"gene" %in% colnames(df_obj)) {
    stop("The data frame must contain a 'gene' column.")
  }
  
  # Create gene list for GSEA
  original_gene_list <- df_obj$logFC
  names(original_gene_list) <- df_obj$gene
  
  # Omit NA values
  gene_list <- na.omit(original_gene_list)
  
  # Sort the list in decreasing order (required for clusterProfiler)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Perform GSEA analysis (ensure TERM2GENE, e.g., 'C1', is properly defined)
  gsea_results <- GSEA(gene_list, 
                       TERM2GENE = geneset_msigdb_df,  # Replace 'C1' with the actual gene set data frame
                       verbose = FALSE, 
                       pvalueCutoff = 0.05, 
                       maxGSSize = 5000, 
                       minGSSize = 5, 
                       pAdjustMethod = "BH", 
                       eps = 0)
  
  write.xlsx(gsea_results[gsea_results$Description], "chr_enrichment_LeuA_vs_Ctrl.xlsx")
  # Assuming 'gsea_results' is your GSEA result
  chr15_enrichment <- gsea_results[gsea_results$Description == "chr1" | 
                                     grepl("chr1", gsea_results$Description), ]
  write.xlsx(chr15_enrichment, "chr1_enrichment_LeuA_vs_Ctrl.xlsx")
  # Check the results for chromosome 1
  #print(chr1_enrichment)
  getwd()
  pdf(paste0('GSEA/MSigDBchr15_LeuA_vs_Ctrll', i, '.pdf'), width = 6, height = 5)
  gseaplot2(gsea_results, geneSetID = which(gsea_results$Description == "chr1"), 
            title = "Chromosome 1 Enrichment")
  dev.off()
  # Save GSEA plot to PDF
  library(enrichplot)
  pdf(paste0('GSEA/MSigDB_LeuA_vs_Ctrl', i, '.pdf'), width = 6, height = 5)
  gseaplot2(gsea_results, title = "MSigDB", geneSetID = 1)
  dev.off()
}



