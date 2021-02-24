#!/usr/bin/env Rscript

# Normalizing the raw data 
suppressMessages(library(data.table))
suppressMessages(library(caret))
suppressMessages(library(yaGST))
suppressMessages(library(Miracle))

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/Model_building/")

deg_df <- read_excel("Required_Files/predictive genes.xlsx")

load("Required_Files/revised_normalized-log2-count.RData")

resMWW <- c()
for(i in 1:ncol(normalized.log2.count)){
  resMWW <- c(resMWW,mwwGST(normalized.log2.count[,i],deg_df$gene)$nes) 
}
resMWW
