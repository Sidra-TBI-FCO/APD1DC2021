#script for ICR clustering of RNAseq data from subreads per cancer 
#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("ConsensusClusterPlus","clue")
ipak(required.packages)

# load data
data = read.csv(paste0("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/GRCh37ERCC_refseq105_genes_count.csv"))
load("~/Sidra Medicine - Research Division/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/ICR genes/ICR_genes.RData")

# Create directories
dir.create("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/Analysis/",showWarnings = FALSE)

rownames(data) = data$X
data$X = NULL
data = as.matrix(data)
# Save the matrix in the correct format for further analysis in which gene names in the row names and samples are column names 
save(data,file = "~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/Synthetic_gene.expression.matrix.Rdata")

#Calculating the ICRscore:
# subset expression matix for ICR
ICR_subset_RNAseq = t(data[which(row.names(data) %in% ICR_genes), ])  #t=transpose :columns will be the rows and rows will be the columns # which rownames are in the ICR genes (the position)
ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)
clustering = data.frame(ICRscore = rowMeans(ICR_subset_RNAseq_log2))  # ICR score is the mean expressions of the rows(ICR genes), create a dataframe.

save(clustering, file = paste0("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/synthetic_table_cluster.RData"))
