library(openxlsx)
setwd("/Users/maurizio.aurora/Downloads")
dgeResults <- read.xlsx("Supplemental_Table_S1.xlsx", sheet = "LeuA vs Ctrl")
nrow(dgeResults)
#dgeResults  <- dgeResults[dgeResults$Pvalue_adj <= 0.05,]

head(dgeResults)
# Install required packages if not already installed
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db")
}
library(org.Mm.eg.db)

head(dgeResults)

# Map GeneID to Chromosome
dgeResults$Chromosome <- mapIds(
  org.Mm.eg.db,
  keys = dgeResults$GeneID,
  column = "CHR",
  keytype = "SYMBOL",
  multiVals = "first"
)

head(dgeResults)

# Select only LeukA and CTRL samples
leukB_samples <- grep("^TCL_Bcell3665_|^TCL_Bcell3666_|^TCL_Bcell3667_|^TCL_Bcell3668_", colnames(dgeResults))
leukA_samples <- grep("^TCL_Bcell3658_|^TCL_Bcell3660_|^TCL_Bcell3661_|^TCL_Bcell3662_", colnames(dgeResults))
#ctrl_samples <- grep("^TCL_Bcell3669_|^TCL_Bcell3670_|^TCL_Bcell3671_|^TCL_Bcell3672_", colnames(dgeResults))

#expression_cols <- grep("^TCL_Bcell", colnames(dgeResults))


# Combine the selected columns
expression_cols <- c(leukA_samples, leukB_samples)
#expression_cols <- c(leukA_samples, ctrl_samples)

head(expression_cols)

# Calculate variance for each gene
dgeResults$Variance <- apply(dgeResults[, expression_cols], 1, var)

dgeResults  <- dgeResults[!is.na(dgeResults$Chromosome),]

library(dplyr)

# Aggregate variance by chromosome
chromosome_variance <- dgeResults %>%
  filter(!is.na(Chromosome)) %>% # Remove genes without chromosome info
  group_by(Chromosome) %>%
  summarise(Total_Variance = sum(Variance, na.rm = TRUE)) %>%
  arrange(Chromosome)

nrow(dgeResults)

variance_data <- dgeResults %>%
  group_by(Chromosome) %>%
  summarise(
    TotalVariance = sum(Variance),        # Total variance per chromosome
    MeanVariance = mean(Variance),       # Mean variance
    SDVariance = sd(Variance)            # Standard deviation of variance
  )

write.xlsx(variance_data, "variance_data_BA.xlsx")

write.xlsx(dgeResults, "dgeResults_vBA.xlsx")
library(gtools)

# Ensure Chromosome column is a character vector
variance_data$Chromosome <- as.character(variance_data$Chromosome)

# Sort chromosomes naturally
variance_data$Chromosome <- factor(
  variance_data$Chromosome,
  levels = mixedsort(unique(variance_data$Chromosome))
)

head(variance_data)

library(gtools)
pdf("LeukA_vs_ctrl_total_variance.pdf")
# Plot variance explained with SD bars
ggplot(variance_data, aes(x = Chromosome, y = TotalVariance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_errorbar(
    aes(ymin = TotalVariance - SDVariance, ymax = TotalVariance + SDVariance),
    width = 0.2,
    color = "black"
  ) +
  labs(
    title = "Variance Explained by DEGs per Chromosome",
    x = "Chromosome",
    y = "Total Variance Explained"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

