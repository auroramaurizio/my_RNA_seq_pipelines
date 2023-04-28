#BiocManager::install("NormqPCR") this is an alternative to try
#I run NormFinder for R (https://moma.dk/files/newDocOldStab_v5.pdf and https://pubmed.ncbi.nlm.nih.gov/15289330/) on pcr CTs.
#library(openxlsx)

#all timepoints and tissues together

test = read.xlsx("/Users/maurizio.aurora/Downloads/sample_normfinder.xlsx")
head(test)
Normfinder(test)
?read.xlsx
# Version history: see bottom of file


setwd("/Users/maurizio.aurora/Downloads/Normfinder")
source("r.NormOldStab5.txt")
getwd()
Result=Normfinder("all.txt")
head(Result)
write.xlsx(Result, "all.xlsx", rowNames = TRUE)

#single tissues at different timepoints

Result=Normfinder("BrainB1.txt")
head(Result)
write.xlsx(Result, "BrainB1.xlsx", rowNames = TRUE)


Result=Normfinder("BrainB2.txt")
head(Result)
write.xlsx(Result, "BrainB2.xlsx", rowNames = TRUE)


Result=Normfinder("BrainB3.txt")
head(Result)
write.xlsx(Result, "BrainB3.xlsx", rowNames = TRUE)


Result=Normfinder("Eye.txt")
head(Result)
write.xlsx(Result, "Eye.xlsx", rowNames = TRUE)


Result=Normfinder("Tongue.txt")
head(Result)
write.xlsx(Result, "Tongue.xlsx", rowNames = TRUE)


Result=Normfinder("Cochlea.txt")
head(Result)
write.xlsx(Result, "Cochlea.xlsx", rowNames = TRUE)


Result=Normfinder("Olf.Epith.txt")
head(Result)
write.xlsx(Result, "Olf.Epith.xlsx", rowNames = TRUE)


