
#suppressPackageStartupMessages(library(monocle3)) # OK
suppressPackageStartupMessages(library(monocle)) # OK
suppressPackageStartupMessages(library(SeuratWrappers)) # NO
suppressPackageStartupMessages(library(Seurat)) # OK
#suppressPackageStartupMessages(library(SeuratData)) # NO
suppressPackageStartupMessages(library(ggplot2)) # OK
suppressPackageStartupMessages(library(patchwork)) # OK
suppressPackageStartupMessages(library(magrittr)) # OK
suppressPackageStartupMessages(library(future)) # OK
suppressPackageStartupMessages(library(cowplot)) # OK
suppressPackageStartupMessages(library(dplyr)) # OK
library("ggalluvial")
###################################################################

#saveRDS(integrated, file = "integrated_subset_clust_renamed_30pc.Rds")

setwd(" ")



#load the integrated object
integrated = readRDS(file ="integrated.Rds")
#select the desired assy
DefaultAssay(integrated) = "integrated"
#select the district
i.combined <- subset (integrated,  subset = `district` == 'iLN')


pt <- table(Idents(i.combined), i.combined$development)
pt <- as.data.frame(pt)
pt$Cluster <- as.character(pt$Var1)
pt$Group=paste(pt$Cluster,pt$Var2)

head(pt)

getwd()

library(scales)
show_col(hue_pal()(9))

pdf("inguinal_alluvial_test.pdf")
ggplot(data = pt,
       aes(axis1 = Cluster, axis2 = Var2,
           y = Freq, fill = Cluster, label = Cluster)) +
  scale_x_discrete(limits = c("Cluster", "Development"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill=Cluster)) +
  geom_stratum(aes(fill=Cluster)) +
  theme(legend.position = "none") +
  geom_text(stat = "stratum", size = 5, aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Inguinal")+
  #scale_fill_manual(values=c(embryo='red',P2='green', adult='blue', ciao = 'green'))+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5, face = "bold")) +
  #scale_fill_manual(values=c('0' = '#F8766D', '1' = '#7CAE00', '2' = '#00BFC4', '3' = '#C77CFF')) +
  scale_fill_manual(values=c('0' = '#F8766D', '1' = '#D89000', '2' = '#A3A500', '3' = '#39B600', '4' = '#00BF7D','5' = '#00BFC4', '6' ='#00B0F6', '7' ='#9590FF', '8' = '#E76BF3', '9'="#FF62BC")) +
  theme(legend.title = element_text(size=20, face="bold"))
dev.off()


data(vaccinations)

head(vaccinations)
vaccinations <- transform(vaccinations,
                          response = factor(response, rev(levels(response))))



ggplot(vaccinations,
       aes(x = survey, stratum = response, alluvium = subject,
           y = freq,
           fill = response, label = response)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  ggtitle("vaccination survey responses at three points in time")


getwd()

help(transform)
suppressMessages(library(openxlsx))


intact = read.table("D7_vs_intact.tsv")
embryo = read.table("embryo_vs_adult.tsv")
D7 = read.table("D7_vs_adult.tsv")

embryob = transform(embryo, freq.loc = ave(seq(nrow(embryo)), state, FUN=length))
intactb = transform(intact, freq.loc = ave(seq(nrow(intact)), state, FUN=length))
D7b = transform(D7, freq.loc = ave(seq(nrow(D7)), state, FUN=length))


intactb$Freq = round((intactb$freq.loc/(nrow(intactb)))*100)
embryob$Freq = round((embryob$freq.loc/(nrow(embryob)))*100)
D7b$Freq = round((D7b$freq.loc/(nrow(D7b)))*100)





test = rbind(embryob,D7b)





pdf("embryo_D7_adult.pdf")
ggplot(test,
       aes(x = comparison, stratum = state, alluvium = gene,
           y = Freq,
           fill = state, label = state)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none")
dev.off()


test = rbind(embryob,intactb)

pdf("embryo_intact_adult.pdf")
ggplot(test,
       aes(x = comparison, stratum = state, alluvium = gene,
           y = Freq,
           fill = state, label = state)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none")
dev.off()




integrated_rn = readRDS("integrated_rn.Rds")

table(integrated_rn$stim)
pt <- table(Idents(integrated_rn), integrated_rn$stim)
pt <- as.data.frame(pt)
pt$Cluster <- as.character(pt$Var1)
colnames(pt) = c("Clusters", "Var2", "Freq", "Cluster")
pt =pt[!grepl("8", pt$Cluster),]


pt0 =pt[grepl("A", pt$Var2),]

pt1 =pt[grepl("B", pt$Var2),]

pt2 =pt[grepl("C", pt$Var2),]

pt3 =pt[grepl("D", pt$Var2),]

pt4 =pt[grepl("E", pt$Var2),]

pt5 =pt[grepl("F", pt$Var2),]

pt0$freq  = pt0$Freq/(sum(pt0$Freq)*100)
pt2$freq  = pt2$Freq/(sum(pt2$Freq)*100)
pt1$freq  = pt1$Freq/(sum(pt1$Freq)*100)
pt3$freq  = pt3$Freq/(sum(pt3$Freq)*100)
pt4$freq  = pt4$Freq/(sum(pt4$Freq)*100)
pt5$freq  = pt5$Freq/(sum(pt5$Freq)*100)

pt6 = rbind(pt0, pt1, pt2, pt3, pt4, pt5)
positions <- c("A","B", "C", "D", "E", "F")


pdf("alluvial_stim.pdf", 20, 10)
ggplot(pt6,
       aes(x = Var2, stratum = Cluster, alluvium = Cluster,
           y = freq,
           fill = Cluster)) +
  #scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha = .9) +
  theme_bw(base_size = 15) +
  scale_x_discrete(limits = positions)+
  geom_flow(stat = "alluvium", lode.guidance = "forward") +
  scale_fill_manual(values=c('ARTERIAL' = '#0066FF',
                             'BARR_END_CAP' = '#336666',
                             'CAPILLARY_PLVAP-' = '#399933',
                             'CAPILLARY_PLVAP+' = '#99CC33',
                             'TIP_1' = '#6600CC',
                             'TIP_2' = '#FF99CC',
                             'TIP_3' = '#FF00FF',
                             '8' = 'grey',
                             'VENOUS_PLVAP-' = '#990000',
                             'VENOUS_PLVAP+' = '#FF6666'))

dev.off()
