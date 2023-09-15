library("ctrlGene")
library("openxlsx")

setwd("/Users/maurizio.aurora/")
## all

all = read.xlsx("/Users/maurizio.aurora/Downloads/all_AllTissues_input.xlsx",
                rowNames = FALSE)
rownames(all) = all$Metadata
write.xlsx(all, "alltissues_fem_male_input.xlsx", rowNames = F)
all <- subset(all, select = -Metadata )
df1 =  as.data.frame(lapply(all,as.numeric))
rownames(df1) = rownames(all)
all = df1
is.na(df1$ActB)
setwd("/Users/maurizio.aurora/geNorm")


which(is.na(df1))
print("Count of total missing values - ")
sum(is.na(df1))


geNorm = geNorm(df1, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_alltissues_MF.xlsx", rowNames = F)

pdf("geNorm_all_MF.pdf")
plotM(geNorm)
dev.off()

setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(all, ctVal = TRUE)
write.xlsx(bk, "bk_alltissues_MF.xlsx", rowNames = T)

## females

all = read.xlsx("/Users/maurizio.aurora/Downloads/all_AllTissues_input.xlsx", 
                rowNames = FALSE)

rownames(all) = all$Metadata
fem = all[grep("^F", all$Metadata), ]
write.xlsx(fem, "alltissues_F_input.xlsx", rowNames = F)
fem <- subset(fem, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(fem, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_alltissues_F.xlsx", rowNames = F)

pdf("geNorm_all_F.pdf")
plotM(geNorm)
dev.off()

setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(fem, ctVal = TRUE)
write.xlsx(bk, "bk_alltissues_F.xlsx", rowNames = T)


## Males

all = read.xlsx("/Users/maurizio.aurora/Downloads/all_AllTissues_input.xlsx", 
                rowNames = FALSE)

rownames(all) = all$Metadata
male = all[grep("^M", all$Metadata), ]
write.xlsx(male, "alltissues_M_input.xlsx", rowNames = F)
male <- subset(male, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(male, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_alltissues_M.xlsx", rowNames = F)
pdf("geNorm_alltissues_male.pdf")
plotM(geNorm)
dev.off()

setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(male, ctVal = TRUE)
write.xlsx(bk, "bk_alltissues_M.xlsx", rowNames = T)



## B1: Olf. Bulb all
all = read.xlsx("/Users/maurizio.aurora/Downloads/all_AllTissues_input.xlsx", 
                rowNames = FALSE)
rownames(all) = all$Metadata
B1 = all[grep("B1", all$Metadata), ]
write.xlsx(B1, "B1_FM_input.xlsx", rowNames = F)
B1 <- subset(B1, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B1, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_FM_B1.xlsx", rowNames = F)
pdf("geNorm_FM_B1.pdf")
plotM(geNorm)
dev.off()

setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B1, ctVal = TRUE)
write.xlsx(bk, "bk_FM_B1.xlsx", rowNames = T)


## B2: Cerebellum all
B2 = all[grep("B2", all$Metadata), ]
write.xlsx(B2, "B2_FM_input.xlsx", rowNames = F)
B2 <- subset(B2, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B2, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_FM_B2.xlsx", rowNames = F)
pdf("geNorm_FM_B2.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B2, ctVal = TRUE)
write.xlsx(bk, "bk_FM_B2.xlsx", rowNames = T)

## B3: mid-Brain 

B3 = all[grep("B3", all$Metadata), ]
write.xlsx(B3, "B3_FM_input.xlsx", rowNames = F)
B3 <- subset(B3, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B3, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_FM_B3.xlsx", rowNames = F)
pdf("geNorm_FM_B3.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B3, ctVal = TRUE)
write.xlsx(bk, "bk_FM_B3.xlsx", rowNames = T)

## E: Eye

temp = all[grep("E", all$Metadata), ]
E = temp[!grepl("OE", temp$Metadata), ]
write.xlsx(E, "E_FM_input.xlsx", rowNames = F)
E <- subset(E, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(E, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_FM_E.xlsx", rowNames = F)
pdf("geNorm_FM_E.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(E, ctVal = TRUE)
write.xlsx(bk, "bk_FM_E.xlsx", rowNames = T)


## T: Tongue

To = all[grep("T", all$Metadata), ]
write.xlsx(To, "To_FM_input.xlsx", rowNames = F)
To <- subset(To, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(To, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_FM_T.xlsx", rowNames = F)
pdf("geNorm_FM_T.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(To, ctVal = TRUE)
write.xlsx(bk, "bk_FM_To.xlsx", rowNames = T)

## C: Cochlea

Co = all[grep("C", all$Metadata), ]
write.xlsx(Co, "Co_FM_input.xlsx", rowNames = F)
Co <- subset(Co, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(Co, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_FM_C.xlsx", rowNames = F)
pdf("geNorm_FM_C.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(Co, ctVal = TRUE)
write.xlsx(bk, "bk_FM_C.xlsx", rowNames = T)

## OE: Olf. Epithelium

OE = all[grep("OE", all$Metadata), ]
write.xlsx(OE, "OE_FM_input.xlsx", rowNames = F)
OE <- subset(OE, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(OE, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_FM_OE.xlsx", rowNames = F)
pdf("geNorm_FM_OE.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(OE, ctVal = TRUE)
write.xlsx(bk, "bk_FM_OE.xlsx", rowNames = T)


##### single tissues only Females

fem = all[grep("^F", all$Metadata), ]
## B1: Olf. Bulb fem
B1 = fem[grep("B1", fem$Metadata), ]
write.xlsx(B1, "B1_fem_input.xlsx", rowNames = F)
B1 <- subset(B1, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B1, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_fem_B1.xlsx", rowNames = F)
pdf("geNorm_fem_B1.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B1, ctVal = TRUE)
write.xlsx(bk, "bk_fem_B1.xlsx", rowNames = T)


## B2: Cerebellum fem
B2 = fem[grep("B2", fem$Metadata), ]
write.xlsx(B2, "B2_fem_input.xlsx", rowNames = F)
B2 <- subset(B2, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B2, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_fem_B2.xlsx", rowNames = F)
pdf("geNorm_fem_B2.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B2, ctVal = TRUE)
write.xlsx(bk, "bk_fem_B2.xlsx", rowNames = T)

## B3: mid-Brain 

B3 = fem[grep("B3", fem$Metadata), ]
write.xlsx(B3, "B3_fem_input.xlsx", rowNames = F)
B3 <- subset(B3, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B3, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_fem_B3.xlsx", rowNames = F)
pdf("geNorm_fem_B3.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B3, ctVal = TRUE)
write.xlsx(bk, "bk_fem_B3.xlsx", rowNames = T)

## E: Eye

temp = fem[grep("E", fem$Metadata), ]
E = temp[!grepl("OE", temp$Metadata), ]
write.xlsx(E, "E_fem_input.xlsx", rowNames = F)
E <- subset(E, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(E, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_fem_E.xlsx", rowNames = F)
pdf("geNorm_fem_E.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(E, ctVal = TRUE)
write.xlsx(bk, "bk_fem_E.xlsx", rowNames = T)


## T: Tongue

To = fem[grep("T", fem$Metadata), ]
write.xlsx(To, "To_fem_input.xlsx", rowNames = F)
To <- subset(To, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(To, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_fem_T.xlsx", rowNames = F)
pdf("geNorm_fem_T.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(To, ctVal = TRUE)
write.xlsx(bk, "bk_fem_T.xlsx", rowNames = T)

## C: Cochlea

Co = fem[grep("C", fem$Metadata), ]
write.xlsx(Co, "Co_fem_input.xlsx", rowNames = F)
Co <- subset(Co, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(Co, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_fem_C.xlsx", rowNames = F)
pdf("geNorm_fem_C.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(Co, ctVal = TRUE)
write.xlsx(bk, "bk_fem_C.xlsx", rowNames = T)

## OE: Olf. Epithelium

OE = fem[grep("OE", fem$Metadata), ]
write.xlsx(OE, "OE_fem_input.xlsx", rowNames = F)
OE <- subset(OE, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(OE, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_fem_OE.xlsx", rowNames = F)
pdf("geNorm_fem_OE.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(OE, ctVal = TRUE)
write.xlsx(bk, "bk_fem_OE.xlsx", rowNames = T)

##### single tissues only Males

## B1: Olf. Bulb male
male = all[grep("^M", all$Metadata), ]
B1 = male[grep("B1", male$Metadata), ]
write.xlsx(B1, "B1_male_input.xlsx", rowNames = F)
B1 <- subset(B1, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B1, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_male_B1.xlsx", rowNames = F)
pdf("geNorm_male_B1.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B1, ctVal = TRUE)
write.xlsx(bk, "bk_male_B1.xlsx", rowNames = T)


## B2: Cerebellum male
B2 = male[grep("B2", male$Metadata), ]
write.xlsx(B2, "B2_male_input.xlsx", rowNames = F)
B2 <- subset(B2, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B2, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_male_B2.xlsx", rowNames = F)
pdf("geNorm_male_B2.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B2, ctVal = TRUE)
write.xlsx(bk, "bk_male_B2.xlsx", rowNames = T)

## B3: mid-Brain 

B3 = male[grep("B3", male$Metadata), ]
write.xlsx(B3, "B3_male_input.xlsx", rowNames = F)
B3 <- subset(B3, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(B3, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_male_B3.xlsx", rowNames = F)
pdf("geNorm_male_B3.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(B3, ctVal = TRUE)
write.xlsx(bk, "bk_male_B3.xlsx", rowNames = T)

## E: Eye

temp = male[grep("E", male$Metadata), ]
E = temp[!grepl("OE", temp$Metadata), ]
write.xlsx(E, "E_male_input.xlsx", rowNames = F)
E <- subset(E, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(E, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_male_E.xlsx", rowNames = F)
pdf("geNorm_male_E.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(E, ctVal = TRUE)
write.xlsx(bk, "bk_male_E.xlsx", rowNames = T)


## T: Tongue

To = male[grep("T", male$Metadata), ]
write.xlsx(To, "To_male_input.xlsx", rowNames = F)
To <- subset(To, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(To, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_male_T.xlsx", rowNames = F)
pdf("geNorm_male_T.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(To, ctVal = TRUE)
write.xlsx(bk, "bk_male_T.xlsx", rowNames = T)

## C: Cochlea

Co = male[grep("C", male$Metadata), ]
write.xlsx(Co, "Co_male_input.xlsx", rowNames = F)
Co <- subset(Co, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(Co, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_male_C.xlsx", rowNames = F)
pdf("geNorm_male_C.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(Co, ctVal = TRUE)
write.xlsx(bk, "bk_male_C.xlsx", rowNames = T)

## OE: Olf. Epithelium

OE = male[grep("OE", male$Metadata), ]
write.xlsx(OE, "OE_male_input.xlsx", rowNames = F)
OE <- subset(OE, select = -Metadata )
setwd("/Users/maurizio.aurora/geNorm")
geNorm = geNorm(OE, ctVal = TRUE)
write.xlsx(geNorm, "geNorm_male_OE.xlsx", rowNames = F)
pdf("geNorm_male_OE.pdf")
plotM(geNorm)
dev.off()
setwd("/Users/maurizio.aurora/bestKeeper")
bk = bestKeeper(OE, ctVal = TRUE)
write.xlsx(bk, "bk_male_OE.xlsx", rowNames = T)


########################

all = read.xlsx("/Users/maurizio.aurora/Downloads/all_AllTissues_input.xlsx",
                rowNames = FALSE)



library(dplyr)
mg_exp_df = (all %>% 
               mutate(tissue = ifelse(Metadata %in% c("C", grep("C", Metadata, value=T)), "C",
                                      ifelse(Metadata %in% c("T", grep("T", Metadata, value=T)), "T",
                                             ifelse(Metadata %in% c("OE", grep("OE", Metadata, value=T)), "OE",
                                                    ifelse(Metadata %in% c("C", grep("C", Metadata, value=T)), "C",
                                                           ifelse(Metadata %in% c("B1", grep("B1", Metadata, value=T)), "B1",
                                                                  ifelse(Metadata %in% c("B2", grep("B2", Metadata, value=T)), "B2",
                                                                         ifelse(Metadata %in% c("B3", grep("B3", Metadata, value=T)), "B3",
                                                   ifelse(Metadata %in% c(grep("E", Metadata, value=T)), "E", NA))))))))))




mg_exp_df$meta_new = gsub("\\#.*","",mg_exp_df$Metadata)
all_new = (mg_exp_df %>% 
               mutate(age = ifelse(meta_new %in% c("12", grep("12", meta_new, value=T)), "12",
                                      ifelse(meta_new %in% c("6", grep("6", meta_new, value=T)), "6",
                                             ifelse(meta_new %in% c("2", grep("2", meta_new, value=T)), "2",
                                                    ifelse(meta_new %in% c("18", grep("18", meta_new, value=T)), "18", NA))))))



all_new_saved = all_new
all_new <- subset(all_new, select = -Metadata )
all_new <- subset(all_new, select = -meta_new )

tra = t(all_new)



write.xlsx(tra, "NormFinder_alltissues_MF_input.xlsx", rowNames = T)

# install 
setwd("/Users/maurizio.aurora/Downloads/Normfinder")
source("r.NormOldStab5.txt")
getwd()


Result=Normfinder("NormFinder_alltissues_MF_input.txt")
head(Result)
write.xlsx(Result, "NormFinder_alltissues_MF.xlsx", rowNames = TRUE)



###

fem = all_new_saved[grep("^F", all_new_saved$Metadata), ]
fem <- subset(fem, select = -Metadata )
fem <- subset(fem, select = -meta_new )

tra = t(fem)
tail(tra)
write.xlsx(tra, "NormFinder_alltissues_F_input.xlsx", rowNames = T)

Result=Normfinder("NormFinder_alltissues_F_input.txt")

write.xlsx(Result, "NormFinder_alltissues_F.xlsx", rowNames = T)


#####

male = all_new_saved[grep("^M", all_new_saved$Metadata), ]
male <- subset(male, select = -Metadata )
male <- subset(male, select = -meta_new )

tra = t(male)
tail(tra)
write.xlsx(tra, "NormFinder_alltissues_M_input.xlsx", rowNames = T)

Result=Normfinder("NormFinder_alltissues_M_input.txt")

write.xlsx(Result, "NormFinder_alltissues_M.xlsx", rowNames = T)

head(all_new_saved)
