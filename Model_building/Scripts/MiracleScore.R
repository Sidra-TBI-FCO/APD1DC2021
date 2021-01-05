suppressMessages(library(yaGST))
suppressMessages(library(Miracle))

#setwd("DreamChallenge_PD1/APD1DC2021/Model_building/")

print("Starting calculating Miracle Scores")
#Load normalized data
load("/data/normalized-log2-count.Rdata")

##Perform Miracle
Mir_res_ALL <- Calculate_Miracle(normalized.log2.count, platform = "gene")  #available platforms: ens", "u133p2", "entrez", "gene"

## Return Miracle results
Mir_res <- Mir_res_ALL$Miracle
names(Mir_res) <- rownames(Mir_res_ALL)

miraclescore_sig = data.frame("patientID" = rownames(Mir_res_ALL) ,"prediction"= Mir_res)

#write.csv(miraclescore_sig, file = "./Model_building/Processed_data/predictions_miraclescore.csv", quote = F, row.names = F); 
write.csv(miraclescore_sig, file = "/output/predictions_miraclescore.csv", quote = F, row.names = F); 

print("Done writing out miracle score signature")