suppressMessages(library(yaGST))
suppressMessages(library(Miracle))

setwd("DreamChallenge_PD1/APD1DC2021/Model_building/")

#Run Miracleon Riaz datasets

#Load normalized data
load("Required_Files/normalized-log2-count.RData")

##Perform Miracle
Mir_res_ALL <- Calculate_Miracle(normalized.log2.count, platform = "gene")  #available platforms: ens", "u133p2", "entrez", "gene"

## Return Miracle results
Mir_res <- Mir_res_ALL$Miracle
names(Mir_res) <- rownames(Mir_res_ALL)
