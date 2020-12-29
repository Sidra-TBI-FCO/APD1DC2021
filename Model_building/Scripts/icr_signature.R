#!/usr/bin/env Rscript

#Use the ICR signature to make the score to be used as prediction
suppressMessages(library(data.table))
suppressMessages(library(NOISeq))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(clue))

# get the input rna-seq gene level count data
print("Processing RNA-Seq data")
#load("./Model_building/Required_Files/normalized-log2-count.Rdata") # gene count normalized matrix
load("/data/normalized-log2-count.Rdata")
#load("./Model_building/Required_Files/ICR_genes.RData")
load("/data/ICR_genes.RData")
print("Done reading in counts")

print("Computing ICR Score")
# Subset to get ICR genes 
ICR_subset_RNAseq = normalized.log2.count[which(rownames(normalized.log2.count) %in% ICR_genes),]

ICR_subset_RNAseq = t(ICR_subset_RNAseq) # transpose
print("Done computing ICR Score")

# calculating ICR score (row mean)
ICR_score = rowMeans(ICR_subset_RNAseq) 

icrscore_sig = data.frame("patientID" = colnames(normalized.log2.count) ,"prediction"=ICR_score)

# write out ICRScore signature to prediciton file
#write.csv(icrscore_sig, file = "./Model_building/Processed_data/predictions_ICR.csv", quote = F, row.names = F); 
write.csv(icrscore_sig, file = "/output/predictions_ICR.csv", quote = F, row.names = F); 

print("Done writing out signature")                          
rm(ICR_subset_RNAseq,ICR_score,normalized.log2.count)
