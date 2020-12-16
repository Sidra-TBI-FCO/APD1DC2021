# Script for the Deconvolution and GSEA
# Setup environment
rm(list=ls())

setwd("~/Sidra Medicine - Research Division/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET")                                                                    # Setwd to location were output files have to be saved.
source("~/Sidra Medicine - Research Division/Sidra Medicine - Research Division/TBI-LAB - General/Bioinformatics tools/R scripts/ipak.function.R")

required.bioconductor.packages = c("GSVA","gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameter
Gene.set = "Selected.pathways"
# Load data 
load(paste0("~/Sidra Medicine - Research Division/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/Selected Pathways/Selected.pathways.3.4.RData"))
normalized.data = load("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/Synthetic_gene.expression.matrix.Rdata")

#Setting the working directory to the location you want 
setwd("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data")

Expression.data = log(data +1, 2)
available_genes = rownames(Expression.data)
Gene.list = Selected.pathways
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Expression.data))]

cat(paste0(Gene.set," ssGSEA ", ". Total number of genes is ", length(unlist(Gene.list)), ".",
           " Of which ", length(unlist(Gene.list)[unlist(Gene.list) %in% available_genes]), 
           " genes are available in expression data."), append = TRUE, sep = "\n")

#Expression.data = as.matrix(Expression.data)
## ssGSEA
ES = gsva(Expression.data,Gene.list,method="ssgsea")   # ES is the enrichment score

signatures = c("[HM] TGF beta signaling","[LM] Proliferation")  
ES  = ES[which(rownames(ES) %in% signatures),]

## Create folders
dir.create("./Analysis/",showWarnings = FALSE)
dir.create("./Analysis/Signature_Enrichment", showWarnings = FALSE)

## Save Enrichment scores
save(ES, 
     file = paste0("./Analysis/004.ssGSEA", Gene.set,".Rdata"))