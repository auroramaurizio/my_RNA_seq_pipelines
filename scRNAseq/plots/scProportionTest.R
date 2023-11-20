
#devtools::install_github("rpolicastro/scProportionTest")
library("scProportionTest")
library(Seurat)
library(ggplot2)

stroma2 = readRDS("stroma2_281023.RDS")
#stroma2$customclassif = Idents(stroma2)
prop_test <- sc_utils(stroma2)


prop_test <- permutation_test(
  prop_test, cluster_identity = "customclassif",
  sample_1 = "leuk-B", sample_2 = "leuk-A",
  sample_identity = "cond")

pdf("scProportionTest_stroma_B_vs_A.pdf")
permutation_plot(prop_test)
dev.off()

prop_test <- sc_utils(stroma2)

prop_test <- permutation_test(
  prop_test, cluster_identity = "customclassif",
  sample_1 = "leuk-B", sample_2 = "Ctrl",
  sample_identity = "cond")

pdf("scProportionTest_stroma_B_vs_Ctrl.pdf")
permutation_plot(prop_test)
dev.off()

prop_test <- sc_utils(stroma2)

prop_test <- permutation_test(
  prop_test, cluster_identity = "customclassif",
  sample_1 = "leuk-A", sample_2 = "Ctrl",
  sample_identity = "cond")

pdf("scProportionTest_stroma_A_vs_Ctrl.pdf")
permutation_plot(prop_test)
dev.off()