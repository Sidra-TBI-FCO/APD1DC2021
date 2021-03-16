#!/usr/bin/env Rscript

# Normalizing the raw data 
suppressMessages(library(data.table))
suppressMessages(library(yaGST))
suppressMessages(library(Miracle))

load("/data1/revised_normalized-log2-count.RData")
load("/data1/ICR_genes.RData")

resMWW <- c()
for(i in 1:ncol(normalized.log2.count)){
  resMWW <- c(resMWW,mwwGST(normalized.log2.count[,i],ICR_genes)$nes) 
}

final_df <- data.frame(patientID = colnames(normalized.log2.count), prediction = -resMWW)

write.table(final_df,file="/output/predictions.csv",quote=F,sep=",",row.names=F,col.names=T)

