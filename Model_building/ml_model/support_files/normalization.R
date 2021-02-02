#!/usr/bin/env Rscript

# Normalizing the raw data 

# Load Library
suppressMessages(library(EDASeq))
suppressMessages(library(base64enc))
suppressMessages(library(preprocessCore))

# load input data
print("Loading Required Data")
#gene.count = read.csv("./Synthetic_Data/GRCh37ERCC_refseq105_genes_count.csv") # raw gene count data
gene.count = read.csv("/data/GRCh37ERCC_refseq105_genes_count.csv")
#load("./Model_building/Required_Files/geneInfo.July2017.RData") # Load geneInfo file
load("/data1/geneInfo.July2017.RData")

# Rename rownames
rownames(gene.count) = gene.count$X
gene.count$X = NULL

# convert gene.count into matrix 
count.matrix = as.matrix(gene.count)

# Remove Na 
geneInfo = geneInfo[!is.na(geneInfo[,1]),]

#Extract genes from geneInfo that matches genes in samples.matrix
# unique = returns a vector, data frame or array like x but with duplicate elements/rows removed.
common.genes = unique(rownames(count.matrix)[which(rownames(count.matrix) %in% rownames(geneInfo))])

# Extract the common genes for geneInfo and samples.matrix
geneInfo = geneInfo[common.genes,]
genes.filtered = count.matrix[common.genes,]

mode(genes.filtered) = "numeric"

count.normalized = newSeqExpressionSet(genes.filtered, featureData = geneInfo)
fData(count.normalized)[,"gcContent"] = as.numeric(geneInfo[,"gcContent"])

print("Normalization")
#within and between lane normalization
# removes lane gene-specific effects, for example effects related to gene length or GC content
count.normalized = withinLaneNormalization(count.normalized, "gcContent", which = "upper", offset = T)

# removes effect related to between lane distributional differences, as sequencing depth
count.normalized = betweenLaneNormalization(count.normalized, which = "upper", offset = T)

print("Log2 Transformation")
#Take= log (unnormalized + .1) + offst(normalized)
count.norm.log = log(genes.filtered +.1) + offst(count.normalized)
count.norm.log = floor(exp(count.norm.log) - .1)  #return non decimal values

#Quantile Normalization
count.quantiles.norm = normalize.quantiles(count.norm.log)
count.quantiles.norm = floor(count.quantiles.norm)

rownames(count.quantiles.norm) = rownames(count.norm.log)
colnames(count.quantiles.norm) = colnames(count.norm.log)

#Log transformation (final normalized-transformed file)
normalized.log2.count = log(count.quantiles.norm+1,2) #log base 2

# final normalized-log2 transformed RData 
#save(normalized.log2.count, file = paste0("./Model_building/Required_Files/normalized-log2-count.RData"))
save(normalized.log2.count, file = paste0("/data1/normalized-log2-count.RData"))
print("Done Normalization")
