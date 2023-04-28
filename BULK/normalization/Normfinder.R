
#I run NormFinder for R (https://moma.dk/files/newDocOldStab_v5.pdf and https://pubmed.ncbi.nlm.nih.gov/15289330/) on pcr CTs.
library(openxlsx)

#BiocManager::install("NormqPCR") # this is an alternative to try



# install 
setwd("/Users/maurizio.aurora/Downloads/Normfinder")
source("r.NormOldStab5.txt")
getwd()

#place input files in the NormFinder directory together with r.NormOldStab5.txt
#try NormFinder on the test dataset sample_normfinder.txt

Result=Normfinder("sample_normfinder.txt")


#all timepoints and tissues together

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


