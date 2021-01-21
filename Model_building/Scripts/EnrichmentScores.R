suppressMessages(library(yaGST))

#setwd("DreamChallenge_PD1/APD1DC2021/Model_building/")

print("Starting calculating Miracle Scores")
#Load normalized data
load("/data1/Selelected_Path_VariousGeneIDs.RData")

## Calculate enrichment scores using gene Symbols. Note that you can use other gene annotations
resMWW <- c()
for(i in 1:ncol(normalized.log2.count)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(normalized.log2.count[,i], x)$nes )))
}
## Return enrichments 
resMWW

write.csv(resMWW, file = "/output/enrichments_scores.csv", quote = F, row.names = F); 

print("Done writing out enrichment scores")
