library(CellChat)
options(stringsAsFactors = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = './'
)
knitr::opts_chunk$set(eval = FALSE)


CellChatDB <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo
write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
write.csv(complex_input, file = "complex_input_CellChatDB.csv")
write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")



#######
options(stringsAsFactors = FALSE)
interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")


#########


setwd("/Users/maurizio.aurora/miniconda3/lib/R/library/CellChat") # This is the folder of CellChat package downloaded from Github
CellChatDB.mouse <- CellChatDB
usethis::use_data(CellChatDB.mouse, overwrite = TRUE)

interact = CellChatDB.mouse$interaction
interact

interact[grep("Plxnd1", interact$interaction_name_2), ]




# set the used database in the object
cellchat_injury@DB <- CellChatDB.mouse
# subset the expression data of signaling genes for saving computation cost
cellchat_injury <- subsetData(cellchat_injury) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat_injury <- identifyOverExpressedGenes(cellchat_injury)
cellchat_injury <- identifyOverExpressedInteractions(cellchat_injury)
# project gene expression data onto PPI network (optional)
cellchat_injury <- projectData(cellchat_injury, PPI.mouse)
cellchat_injury <- computeCommunProb(cellchat_injury)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_injury <- filterCommunication(cellchat_injury, min.cells = 10)
cellchat_injury <- computeCommunProbPathway(cellchat_injury)
cellchat_injury <- aggregateNet(cellchat_injury)


cellchat_intact@DB <- CellChatDB.mouse
# subset the expression data of signaling genes for saving computation cost
# subset the expression data of signaling genes for saving computation cost
cellchat_intact <- subsetData(cellchat_intact) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat_intact <- identifyOverExpressedGenes(cellchat_intact)
cellchat_intact <- identifyOverExpressedInteractions(cellchat_intact)
# project gene expression data onto PPI network (optional)
cellchat_intact <- projectData(cellchat_intact, PPI.mouse)
cellchat_intact <- computeCommunProb(cellchat_intact)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_intact <- filterCommunication(cellchat_intact, min.cells = 10)
cellchat_intact <- computeCommunProbPathway(cellchat_intact)
cellchat_intact <- aggregateNet(cellchat_intact)



object.list <- list(intact = cellchat_intact, injury = cellchat_injury)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
 