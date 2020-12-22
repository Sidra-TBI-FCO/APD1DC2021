#!/usr/bin/env Rscript

#Use the ICR signature to make the score to be used as prediction
suppressMessages(library(data.table))
suppressMessages(library(NOISeq))
suppressMessages(library(GSVA))
suppressMessages(library(gclus))

# get the input rna-seq gene level count data
print("Processing RNA-Seq data")
counts <- fread("/data/GRCh37ERCC_refseq105_genes_count.csv",data.table=F)
#counts <- fread("/export/cse02/rmall/AntiPD1_Challenge/Anti-PD1-DREAM-Examples/baseline/support_files/GRCh37ERCC_refseq105_genes_count.csv",data.table=F)
rownames(counts) <- counts[,1]
counts <- counts[,-1];
counts <- as.data.frame(counts)
print("Done reading in counts")

#Calculating the Pathways Score
load("/Required_Files/Selected.pathways.3.4.RData")
#load("/export/cse02/rmall/AntiPD1_Challenge/Anti-PD1-DREAM-Examples/pathways_score/support_files/Selected.pathways.3.4.RData")

print("Computing Pathways Score")
# Set parameter
Gene.set = "Selected.pathways"
tmm <- NOISeq::tmm(counts); 
Expression.data = log(tmm +1, 2)
available_genes = rownames(Expression.data)
Gene.list = Selected.pathways
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Expression.data))]

cat(paste0(Gene.set," ssGSEA ", ". Total number of genes is ", length(unlist(Gene.list)), ".",
           " Of which ", length(unlist(Gene.list)[unlist(Gene.list) %in% available_genes]), 
           " genes are available in expression data."), append = TRUE, sep = "\n")

## ssGSEA
ES = gsva(Expression.data,Gene.list,method="ssgsea")   # ES is the enrichment score
signatures = c("[HM] TGF beta signaling","[LM] Proliferation")  
ES  = ES[which(rownames(ES) %in% signatures),]
print("Done computing ES scores for TGF beta and Proliferation Pathways")

pathway_score  <- data.frame("patientID" = colnames(ES) ,"prediction"=as.numeric(t(ES)[,1]))

# write out TGF-beta pathways enrichment score to prediciton file
write.csv(pathway_score, file = "/output/predictions.csv", quote = F, row.names = F); 
#write.csv(icrscore_sig, file = "/export/cse02/rmall/AntiPD1_Challenge/Anti-PD1-DREAM-Examples/output/predictions_2.csv", quote=F, row.names = F);
print("Done writing out signature")

