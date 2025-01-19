# Load necessary libraries
library(dplyr)
library(ggplot2)
library(reshape2)



# Example Input Data
# Replace these data frames with your actual data
# DEG list with associated chromosome info
degs_up <- readLines("FDRup_YES_vs_NO.txt")
degs_up_df <- as.data.frame(degs_up)
colnames(degs_up_df) <- "Geneid"
tail(degs_up_df)

# Background list of all expressed genes with chromosome info
background <- read.table("/Users/maurizio.aurora//wgEncodeGencodeBasicVM22.txt")

background_bk <- background

background <-  unique(background[,c("V3","V13")])
colnames(background) <- c("chr", "Geneid")
background <- background %>%
  filter(!grepl("_random", chr))

background <- background %>%
  filter(!grepl("chrUn", chr))


degs <- unique(left_join(degs_up_df, background, by = "Geneid"))
nrow(degs)
nrow(degs_up_df)

# Step 1: Count DEGs and background genes per chromosome

deg_counts <- degs %>%
  group_by(chr) %>%
  summarise(deg_count = n())

background_counts <- background %>%
  group_by(chr) %>%
  summarise(total_genes = n())


# Step 2: Merge DEG and background counts
merged_data <- left_join(background_counts, deg_counts, by = "chr") %>%
  mutate(deg_count = ifelse(is.na(deg_count), 0, deg_count)) # Handle chromosomes without DEGs

# Step 3: Compute expected DEG counts based on background proportions
total_degs <- nrow(degs_up_df)
total_background_genes <- sum(merged_data$total_genes)
merged_data <- merged_data %>%
  mutate(expected_deg_count = (total_genes / total_background_genes) * total_degs)


# Order the data by chromosome (if not already ordered)
merged_data$chr <- factor(merged_data$chr, 
                          levels = c(paste0("chr", 1:19), "chrX", "chrY", "chrM"))

data <- merged_data
# Create the plot and save as PDF
pdf("chromosome_total_genes_hist_ordered.pdf", width = 14, height = 6)

ggplot(data, aes(x = chr, y = total_genes)) +
  geom_bar(stat = "identity", aes(fill = ifelse(chr == "chr1", "chr1", "other")), color = "black") + 
  scale_fill_manual(values = c("chr1" = "pink", "other" = "black")) +
  labs(
    title = "Distribution of Total Genes by Chromosome",
    x = "Chromosome",
    y = "Total Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

sum(background_counts$total_genes)

background_chr1 <- 1462  # Total genes on chr1 in background (replace with actual number)
background_other <- 41846 - background_chr1  # Total background genes minus chr1 genes

deg_counts <- degs %>%
  group_by(chr) %>%
  summarise(deg_count = n())

background_counts <- background %>%
  group_by(chr) %>%
  summarise(total_genes = n())

up_chr1 <- 11
up_other <- 95-11
# Create the contingency table
contingency_table <- matrix(c(up_chr1, up_other, 
                              background_chr1 - up_chr1, background_other - up_other),
                            nrow = 2, 
                            byrow = TRUE)

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)

# View results
print(fisher_test_result)


> print(fisher_test_result)

#Fisher's Exact Test for Count Data


pdf("chromosome_UP_DEGS_hist_ordered.pdf", width = 14, height = 6)

ggplot(data, aes(x = chr, y = deg_count)) +
  geom_bar(stat = "identity", aes(fill = ifelse(chr == "chr1", "chr1", "other")), color = "black") + 
  scale_fill_manual(values = c("chr1" = "pink", "other" = "black")) +
  labs(
    title = "Distribution of UP DEGs  by Chromosome",
    x = "Chromosome",
    y = "Total Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()





# Create a new column for coloring based on conditions
data$deg_color <- ifelse(data$chr == "chr1", "pink", "black")
data$expected_deg_color <- ifelse(data$chr == "chr1", "purple", "grey")


pdf("chromosome_UP_DEGS_hist_comp_ordered.pdf", width = 14, height = 6)
# Create the plot
ggplot(data, aes(x = chr)) +
  geom_bar(aes(y = deg_count, fill = deg_color), stat = "identity", position = position_dodge(width = 0.8), width = 0.4) +
  geom_bar(aes(y = expected_deg_count, fill = expected_deg_color), stat = "identity", position = position_dodge(width = 0.8), width = 0.4) +
  labs(
    title = "Chromosomal Distribution of DEGs vs Expected DEGs",
    x = "Chromosome",
    y = "Count",
    fill = "Group"
  ) +
  theme_minimal() +
  scale_fill_identity() +  # Use the defined fill colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()



head(data)

pdf("chromosome_UP_DEGS_hist_comp_ordered.pdf", width = 14, height = 6)
# Reshape Data to Long Format
data_long <- data %>%
  pivot_longer(
    cols = c(deg_count, expected_deg_count),
    names_to = "type",
    values_to = "count"
  )

# Assign Colors
data_long$type_color <- ifelse(data_long$type == "deg_count", data$deg_color, data$expected_deg_color)

# Barplot
ggplot(data_long, aes(x = chr, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(
    values = c("deg_count" = "black", "expected_deg_count" = "grey"),
    labels = c("DEG Count", "Expected DEG Count")
  ) +
  labs(x = "Chromosome", y = "Count", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#############


# Example Input Data
# Replace these data frames with your actual data
# DEG list with associated chromosome info
degs_dw <- readLines("/Users/maurizio.aurora/Dropbox (HSR Global)/Omics_2024/Integration_bulkRNAseq_WES_ATAC/Results/enrichR/FDRdw_YES_vs_NO.txt")
degs_dw_df <- as.data.frame(degs_dw)
colnames(degs_dw_df) <- "Geneid"
tail(degs_dw_df)


# Background list of all expressed genes with chromosome info
background <- read.table("/Users/maurizio.aurora//wgEncodeGencodeBasicVM22.txt")

background_bk <- background

background <-  unique(background[,c("V3","V13")])
colnames(background) <- c("chr", "Geneid")
background <- background %>%
  filter(!grepl("_random", chr))

background <- background %>%
  filter(!grepl("chrUn", chr))


degs <- unique(left_join(degs_dw_df, background, by = "Geneid"))
nrow(degs)
nrow(degs_dw_df)

# Step 1: Count DEGs and background genes per chromosome

deg_counts <- degs %>%
  group_by(chr) %>%
  summarise(deg_count = n())

background_counts <- background %>%
  group_by(chr) %>%
  summarise(total_genes = n())

tail(background_counts)

# Step 2: Merge DEG and background counts
merged_data <- left_join(background_counts, deg_counts, by = "chr") %>%
  mutate(deg_count = ifelse(is.na(deg_count), 0, deg_count)) # Handle chromosomes without DEGs

# Step 3: Compute expected DEG counts based on background proportions
total_degs <- nrow(degs_dw_df)
total_background_genes <- sum(merged_data$total_genes)
merged_data <- merged_data %>%
  mutate(expected_deg_count = (total_genes / total_background_genes) * total_degs)





library(ggplot2)
library(dplyr)

# Order the data by chromosome (if not already ordered)
merged_data$chr <- factor(merged_data$chr, 
                          levels = c(paste0("chr", 1:19), "chrX", "chrY", "chrM"))

data <- merged_data
# Create the plot and save as PDF
pdf("chromosome_total_genes_hist_ordered.pdf", width = 14, height = 6)

ggplot(data, aes(x = chr, y = total_genes)) +
  geom_bar(stat = "identity", aes(fill = ifelse(chr == "chr1", "chr1", "other")), color = "black") + 
  scale_fill_manual(values = c("chr1" = "pink", "other" = "black")) +
  labs(
    title = "Distribution of Total Genes by Chromosome",
    x = "Chromosome",
    y = "Total Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()



pdf("chromosome_dw_DEGS_hist_ordered.pdf", width = 14, height = 6)

ggplot(data, aes(x = chr, y = deg_count)) +
  geom_bar(stat = "identity", aes(fill = ifelse(chr == "chr1", "chr1", "other")), color = "black") + 
  scale_fill_manual(values = c("chr1" = "pink", "other" = "black")) +
  labs(
    title = "Distribution of UP DEGs  by Chromosome",
    x = "Chromosome",
    y = "Total Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

#################################################################




# Create a new column for coloring based on conditions
data$deg_color <- ifelse(data$chr == "chr1", "pink", "black")
data$expected_deg_color <- ifelse(data$chr == "chr1", "purple", "grey")

head(data)

pdf("chromosome_dw_DEGS_hist_comp_ordered.pdf", width = 14, height = 6)
# Create the plot
ggplot(data, aes(x = chr)) +
  geom_bar(aes(y = deg_count, fill = deg_color), stat = "identity", position = position_dodge(width = 0.8), width = 0.4) +
  geom_bar(aes(y = expected_deg_count, fill = expected_deg_color), stat = "identity", position = position_dodge(width = 0.8), width = 0.4) +
  labs(
    title = "Chromosomal Distribution of DEGs vs Expected DEGs",
    x = "Chromosome",
    y = "Count",
    fill = "Group"
  ) +
  theme_minimal() +
  scale_fill_identity() +  # Use the defined fill colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

#############



pdf("chromosome_dw_DEGS_hist_comp_ordered.pdf", width = 14, height = 6)
# Reshape Data to Long Format
data_long <- data %>%
  pivot_longer(
    cols = c(deg_count, expected_deg_count),
    names_to = "type",
    values_to = "count"
  )

# Assign Colors
data_long$type_color <- ifelse(data_long$type == "deg_count", data$deg_color, data$expected_deg_color)

# Barplot
ggplot(data_long, aes(x = chr, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(
    values = c("deg_count" = "black", "expected_deg_count" = "grey"),
    labels = c("DEG Count", "Expected DEG Count")
  ) +
  labs(x = "Chromosome", y = "Count", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

####



# Step 4: Perform Chi-square test
chisq_result <- chisq.test(merged_data$deg_count, p = merged_data$total_genes / total_background_genes)

# Print results
print("Chi-Square Test Results:")
print(chisq_result)

# Step 5: Visualization
# Normalize DEG and background counts for visualization
merged_data <- merged_data %>%
  mutate(
    normalized_deg_count = deg_count / total_genes,
    normalized_background = total_genes / total_background_genes
  )




pdf("merged_data_chr.pdf", 14, 6)
ggplot(merged_data, aes(x = chr)) +
  geom_bar(aes(y = normalized_deg_count, fill = "DEGs"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = normalized_background, fill = "Background"), stat = "identity", position = "dodge") +
  labs(
    title = "Chromosomal Distribution of DEGs",
    x = "Chromosome",
    y = "Proportion",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate labels by 45 degrees
  )
dev.off()




