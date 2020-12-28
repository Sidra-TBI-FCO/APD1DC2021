#!/usr/bin/env Rscript

#Use the ICR signature to make the score to be used as prediction
suppressMessages(library(data.table))
suppressMessages(library(NOISeq))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(clue))

# get the input rna-seq gene level count data
print("Processing RNA-Seq data")
counts <- fread("/data/GRCh37ERCC_refseq105_genes_count.csv",data.table=F)
#counts <- fread("/export/cse02/rmall/AntiPD1_Challenge/Anti-PD1-DREAM-Examples/baseline/support_files/GRCh37ERCC_refseq105_genes_count.csv",data.table=F)

rownames(counts) <- counts[,1]
counts <- counts[,-1];
counts <- as.data.frame(counts)
print("Done reading in counts")

#Calculating the ICR Score
load("/data/ICR_genes.RData")
#load("/export/cse02/rmall/AntiPD1_Challenge/Anti-PD1-DREAM-Examples/icrscore/support_files/ICR_genes.RData")

print("Computing ICR Score")
tmm <- NOISeq::tmm(counts); 
ICR_subset_RNAseq = t(counts[which(rownames(tmm) %in% ICR_genes), ])  #t=transpose :columns will be the rows and rows will be the columns # which rownames are in the ICR genes (the position)
ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)
clustering = rowMeans(ICR_subset_RNAseq_log2)  # ICR score is the mean expressions of the rows(ICR genes), create a dataframe.
print("Done computing ICR Score")

icrscore_sig  <- data.frame("patientID" = colnames(counts) ,"prediction"=clustering)

# write out ICRScore signature to prediciton file
write.csv(icrscore_sig, file = "/output/predictions_ICR.csv", quote = F, row.names = F); 
#write.csv(icrscore_sig, file = "/export/cse02/rmall/AntiPD1_Challenge/Anti-PD1-DREAM-Examples/output/predictions_2.csv", quote=F, row.names = F);

print("Done writing out signature")
rm(ICR_subset_RNAseq,ICR_subset_RNAseq_log2,clustering,counts)

