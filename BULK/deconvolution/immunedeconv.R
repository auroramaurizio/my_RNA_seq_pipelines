library("immunedeconv")
library(ggplot2)
library(open.xlsx)
library(tidyverse)


#https://omnideconv.org/immunedeconv/articles/immunedeconv.html

TPM <- ovary_mat_sel_TPM
resEpic = deconvolute(TPM, "epic")
ovary_resEpic <- as.data.frame(resEpic)
write.xlsx(ovary_resEpic, "ovary_resEpic.xlsx", rowNames = TRUE)

resQs = deconvolute(TPM, "quantiseq", tumor=TRUE)
ovary_resQs <- as.data.frame(resQs)
write.xlsx(ovary_resQs, "ovary_resQs.xlsx", rowNames = TRUE)


estimate <- deconvolute_estimate(TPM)
ovary_estimate <- as.data.frame(estimate)
write.xlsx(ovary_estimate, "ovary_estimate.xlsx", rowNames = TRUE)




##############

pdf("Ovary_TPM_epic_nolab_sorted.pdf", 10,15)
resEpic %>%
  gather(sample, fraction, -cell_type) %>%
  group_by(sample) %>%
  mutate(uncell_fraction = ifelse(cell_type == "uncharacterized cell", fraction, NA)) %>%
  summarise(uncell_fraction = max(uncell_fraction, na.rm = TRUE)) %>%
  right_join(gather(resEpic, sample, fraction, -cell_type), by = "sample") %>%
  group_by(sample) %>%
  mutate(fraction = fraction / sum(fraction)) %>%
  ungroup() %>%
  arrange(desc(uncell_fraction)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()


pdf("Ovary_TPM_epic_yeslab_sorted.pdf", 25,25)
resEpic %>%
  gather(sample, fraction, -cell_type) %>%
  group_by(sample) %>%
  mutate(uncell_fraction = ifelse(cell_type == "uncharacterized cell", fraction, NA)) %>%
  summarise(uncell_fraction = max(uncell_fraction, na.rm = TRUE)) %>%
  right_join(gather(resEpic, sample, fraction, -cell_type), by = "sample") %>%
  group_by(sample) %>%
  mutate(fraction = fraction / sum(fraction)) %>%
  ungroup() %>%
  arrange(desc(uncell_fraction)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_brewer(palette = "Set3")
dev.off()


pdf("Ovary_TPM_Epic_yeslab_stacked.pdf", 20,20)
resEpic %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Set3") +
  scale_x_discrete(limits = rev(levels(resEpic)))
dev.off()
 


###########



# Reshape the data and sort by the "uncharacterized cell" fraction in descending order
pdf("Ovary_TPM_resQs_nolab_sorted.pdf", 10,10)
resQs %>%
  gather(sample, fraction, -cell_type) %>%
  group_by(sample) %>%
  mutate(uncell_fraction = ifelse(cell_type == "uncharacterized cell", fraction, NA)) %>%
  ungroup() %>%
  arrange(desc(uncell_fraction)) %>%
  ggplot(aes(x = factor(sample, levels = unique(sample)), y = fraction, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

pdf("Ovary_TPM_resQs_yeslab_sorted.pdf", 10,15)
resQs %>%
  gather(sample, fraction, -cell_type) %>%
  group_by(sample) %>%
  mutate(uncell_fraction = ifelse(cell_type == "uncharacterized cell", fraction, NA)) %>%
  ungroup() %>%
  arrange(desc(uncell_fraction)) %>%
  ggplot(aes(x = factor(sample, levels = unique(sample)), y = fraction, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_brewer(palette = "Set3")
dev.off()

pdf("Ovary_TPM_resQs_nolab_stacked.pdf", 10,15)
resQs %>%
  gather(sample, fraction, -cell_type) %>%
  filter(cell_type != "uncharacterized cell") %>%
  group_by(sample) %>%
  mutate(fraction = fraction / sum(fraction)) %>%
  ungroup() %>%
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.y = element_text(size = 10), axis.ticks.y = element_line()) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()


##############




colnames(TPM)

metadata <- read.xlsx("GerSom_SNParray_RNAseq_RCR22.xlsx", sheet = 4)
filtered_metadata <- metadata[metadata$Tumor.Type == "ovary", ]
selected_ids <- colnames(ovary_resQs)
filtered_metadata = unique(filtered_metadata[which(filtered_metadata$AAC_LIMS_BC %in% selected_ids),])
nrow(filtered_metadata)
ncol(estimate_ovary)
head(filtered_metadata)


rownames(filtered_metadata) <- filtered_metadata$AAC_LIMS_BC
filtered_metadata$AAC_LIMS_BC <- NULL
filtered_metadata$Tumor.Type <- NULL
colors <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)


rownames(ovary_resQs) <- ovary_resQs$cell_type
ovary_resQs_df <- as.data.frame(ovary_resQs$cell_type)
colnames(ovary_resQs_df) <- "ovary_res_QS"
rownames(ovary_resQs_df) <- ovary_resQs_df$ovary_res_QS

annotation_row <- ovary_resQs_df
head(annotation_row)
rownames(annotation_row) <- ovary_resQs_df$ovary_res_QS
filtered_metadata[filtered_metadata == "Missing Data"] <- "Missing"

#ovary_resQs[is.na(ovary_resQs)] <- 0

rownames(unique(ovary_resQs))

df <- ovary_resQs[row.names(ovary_resQs) != "uncharacterized cell", ]

df$cell_type <- NULL

filtered_metadata_somatic <- filtered_metadata[, !grepl("germline", names(filtered_metadata))]


pheatmap(df,
         show_colnames = F,
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = T,
         scale = "row",
         annotation_colors = annotation_colors,
         annotation_col = filtered_metadata_somatic,
         annotation_row = annotation_row,
         filename = 'resQs_ovary_heatmap_scaled.pdf',
         width = 40, height = 8,
         display_numbers = F, # Option to display correlation coefficients
         color = colorRampPalette(c("blue", "white", "red"))(50), # Customize colors
         main = "resQs ovary Mutations")


metadata <- read.xlsx("GerSom_SNParray_RNAseq_RCR22.xlsx", sheet = 4)
filtered_metadata <- metadata[metadata$Tumor.Type == "ovary", ]
selected_ids <- colnames(estimate_ovary)
filtered_metadata = unique(filtered_metadata[which(filtered_metadata$AAC_LIMS_BC %in% selected_ids),])
nrow(filtered_metadata)
ncol(estimate_ovary)
head(filtered_metadata)


rownames(filtered_metadata) <- filtered_metadata$AAC_LIMS_BC
filtered_metadata$AAC_LIMS_BC <- NULL
filtered_metadata$Tumor.Type <- NULL
colors <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)


filtered_metadata[filtered_metadata == "Missing Data"] <- "Missing"
colnames(filtered_metadata)


annotation_colors = list(
  MAP3K1_somatic = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  CDK12_somatic = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  RELN_somatic = c("FALSE" = "white",
                   "Missing" = "gray",
                   "TRUE" = "blue"),
  CDH1_somatic =  c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  GATA3_somatic = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  TOP2A_somatic = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  PIK3CA_somatic = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  LRP1B_somatic = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  BRCA1_somatic = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  BRCA2_somatic = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  KMT2C_somatic = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  NF1_somatic = c("FALSE" = "white",
                  "Missing" = "gray",
                  "TRUE" = "blue"),
  RB1_somatic = c("FALSE" = "white",
                  "Missing" = "gray",
                  "TRUE" = "blue"),
  TP53_somatic = c("FALSE" = "white",
                   "Missing" = "gray",
                   "TRUE" = "blue"),
  MAP3K1_germline = c("FALSE" = "white",
                      "Missing" = "gray",
                      "TRUE" = "blue"),
  CDK12_germline =  c("FALSE" = "white",
                      "Missing" = "gray",
                      "TRUE" = "blue"),
  RELN_germline = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"),
  CDH1_germline =  c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  GATA3_germline = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),                                                                                                                                                                        TOP2A_germline = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  PIK3CA_germline = c("FALSE" = "white",
                      "Missing" = "gray",
                      "TRUE" = "blue"),
  LRP1B_germline = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  BRCA1_germline = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  BRCA2_germline = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  KMT2C_germline = c("FALSE" = "white",
                     "Missing" = "gray",
                     "TRUE" = "blue"),
  NF1_germline = c("FALSE" = "white",
                   "Missing" = "gray",
                   "TRUE" = "blue"),
  RB1_germline = c("FALSE" = "white",
                   "Missing" = "gray",
                   "TRUE" = "blue"),
  TP53_germline = c("FALSE" = "white",
                    "Missing" = "gray",
                    "TRUE" = "blue"))





annotation_row <- data.frame(
  ScoreType = c("StromalScore", "ImmuneScore")
)
rownames(annotation_row) <- annotation_row$ScoreType

estimate_ovary <- estimate_ovary[row.names(estimate_ovary) != "ESTIMATEScore", ]
estimate_ovary <- estimate_ovary[row.names(estimate_ovary) != "TumorPurity", ]




pheatmap(estimate_ovary,
         show_colnames = F,
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = T,
         annotation_colors = annotation_colors,
         annotation_col = filtered_metadata_somatic,
         annotation_row = annotation_row,
         filename = 'estimate_ovary_heatmap_not_scaled.pdf',
         width = 40, height = 8,
         display_numbers = F, # Option to display correlation coefficients
         color = colorRampPalette(c("blue", "white", "red"))(50), # Customize colors
         main = "Estimate ovary Mutations")

pheatmap(estimate_ovary,
         show_colnames = F,
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = T,
         scale = "row",
         annotation_colors = annotation_colors,
         annotation_col = filtered_metadata_somatic,
         annotation_row = annotation_row,
         filename = 'estimate_ovary_heatmap_scaled.pdf',
         width = 40, height = 8,
         display_numbers = F, # Option to display correlation coefficients
         color = colorRampPalette(c("blue", "white", "red"))(50), # Customize colors
         main = "Estimate ovary Mutations")




