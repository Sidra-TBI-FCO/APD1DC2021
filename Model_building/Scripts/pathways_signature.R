#!/usr/bin/env Rscript

#Use the TGF-Beta and proliferation to make the enrichment score to be used as prediction
suppressMessages(library(data.table))
suppressMessages(library(NOISeq))
suppressMessages(library(GSVA))
suppressMessages(library(gclus))

# get the input rna-seq gene level count data
print("Processing RNA-Seq data")
#load("./Model_building/Required_Files/normalized-log2-count.RData") # gene count normalized matrix
load("/data1/normalized-log2-count.RData")
#load("./Model_building/Required_Files/Selected.pathways.3.4.RData")
load("/data1/Selected.pathways.3.4.RData")
load("/data1/ICR_genes.RData")
#add ICR to pathways
Selected.pathways$`[TBI] ICR` <- ICR_genes
print("Done reading in counts")


print("Computing Pathways Score")
# Set parameter
Gene.set = "Selected.pathways"
available_genes = rownames(normalized.log2.count)
Gene.list = Selected.pathways
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(normalized.log2.count))]

cat(paste0(Gene.set," ssGSEA ", ". Total number of genes is ", length(unlist(Gene.list)), ".",
           " Of which ", length(unlist(Gene.list)[unlist(Gene.list) %in% available_genes]), 
           " genes are available in expression data."), append = TRUE, sep = "\n")

## ssGSEA
ES = gsva(normalized.log2.count,Gene.list,method="ssgsea")   # ES is the enrichment score
#signatures = c("[HM] TGF beta signaling","[LM] Proliferation")  
#ES  = ES[which(rownames(ES) %in% signatures),]
print("Done computing ES scores")

#pathway_score  <- data.frame("patientID" = colnames(ES) ,"prediction"=as.numeric(t(ES)[,1]))

# write out TGF-beta pathways enrichment score to prediciton file
#write.csv(pathway_score, file = "Model_building/Processed_data/predictions_pathways.csv", quote = F, row.names = F);
#dir.create("./output",showWarnings = FALSE)
#write.csv(pathway_score, file = "./output/predictions_pathways.csv", quote = F, row.names = F); 
save(ES,file = "/data1/patway_score.RData")
print("Done writing out signature")

