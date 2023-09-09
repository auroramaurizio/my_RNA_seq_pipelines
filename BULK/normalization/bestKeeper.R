library("ctrlGene")
library("openxlsx")

expression = read.xlsx("/Users/maurizio.aurora/Downloads/Raw_Ct_rep1_14.4.23_geNorm.xlsx", 
                      sheet = "Olf.Epith",  
                      rowNames = TRUE)

df1 =  as.data.frame(lapply(expression,as.numeric))
rownames(df1) = rownames(expression)

bk = bestKeeper(df1, ctVal = TRUE)



setwd("/Users/maurizio.aurora/bestKeeper")

write.xlsx(bk, "bk_Olf.Epith.xlsx", rowNames = T)
