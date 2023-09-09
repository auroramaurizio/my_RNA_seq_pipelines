library("ctrlGene")
library("openxlsx")

# input should be a CTs dataframe with genes on the columns
# samples on the rows
# head(df1)
#            ACTB       B2M      GAPD 
# FIB1  0.4993028 0.3278251 0.5162567
# FIB2  0.2387128 0.4780923 0.2877965



Olf.Epith = read.xlsx("/Users/maurizio.aurora/Downloads/Raw_Ct_rep1_14.4.23.xlsx", 
                sheet = "Olf.Epith",  
                rowNames = FALSE)
rownames(Olf.Epith) = Olf.Epith$Sample


transposed = as.data.frame(t(Olf.Epith))
rownames(transposed) = paste(rownames(transposed), transposed$Group, sep = "_")
transposed = transposed[-1,]
df <- subset(transposed, select = -Group )
write.xlsx(df, "/Users/maurizio.aurora/Downloads/Raw_Ct_rep1_14.4.23_Olf.Epith.xlsx", rowNames = T)
expression = df
head(expression)
df1 =  as.data.frame(lapply(expression,as.numeric))
rownames(df1) = rownames(expression)

#A sorted dataframe with two columns, 'Genes' and 'Avg.M'. 
#The last two genes are the two most stable control genes.

#Avg.M is average expression stability values (M) of remaining control genes 
#during stepwise exclusion of the least stable control gene.

geNorm = geNorm(df1, ctVal = TRUE)

setwd("/Users/maurizio.aurora/geNorm")

write.xlsx(geNorm, "geNorm_Olf.Epith.xlsx", rowNames = T)

Olf.Epith = read.xlsx("/Users/maurizio.aurora/geNorm/geNorm_Olf.Epith.xlsx", 
                      rowNames = TRUE)

pdf("geNorm_Olf.Epith.pdf")
plotM(Olf.Epith)
dev.off()



