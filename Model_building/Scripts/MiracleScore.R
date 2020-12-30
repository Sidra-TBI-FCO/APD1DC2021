library(yaGST)
library(Miracle)

setwd("DreamChallenge_PD1/APD1DC2021/Model_building/")

#Run Miracleon Riaz datasets

#Load Riaz data
load("Processed_data/Riaz etal/002.Riaz_Data_normalized.gene.counts.Rdata")

##Perform Miracle
Mir_res_ALL <- Calculate_Miracle(dataNorm, platform = "gene")  #available platforms: ens", "u133p2", "entrez", "gene"

## Return Miracle results
Mir_res <- Mir_res_ALL$Miracle
names(Mir_res) <- rownames(Mir_res_ALL)
