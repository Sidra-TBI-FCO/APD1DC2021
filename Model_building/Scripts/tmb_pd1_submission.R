#!/usr/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(yaGST))

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/")

#Load the log2 normalized dataset
#load("/data1/normalized-log2-count.RData")
load("Model_building/Required_Files/revised_normalized-log2-count.RData")

#Load the pathways
load("Model_building/Required_Files/Selelected_Path_VariousGeneIDs.RData")
load("Model_building/Required_Files/ICR_genes.RData")


SelPaths <- list(SelPath_Symb$TMB_Proliferation, ICR_genes)

resMWW <- c()
for(i in 1:ncol(normalized.log2.count)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPaths,function(x) mwwGST(normalized.log2.count[,i], x)$nes )))
}
resMWW <- as.data.frame(resMWW)
colnames(resMWW) <- c("TMB_Proliferation","ICR")

#Get clinical data
clinical_df <- read.table("Synthetic_Data/clinical_data.csv",header=TRUE,sep=",")
clinical_df$patientID <- as.character(as.vector(clinical_df$patientID))


all_tmb <- log2(clinical_df[!is.na(clinical_df$TMB),]$TMB+1)
clinical_df$TMB_Scaled <- log2(clinical_df$TMB+1)
tertiles_tmb <- quantile(all_tmb, probs = seq(0,1,1/3))
tmb_low_cutoff <- tertiles_tmb[[2]]
tmb_high_cutoff <- tertiles_tmb[[3]]

tertiles_icr <- quantile(resMWW$ICR, probs = seq(0,1,1/3))
icr_cutoff <- tertiles_icr[[3]]

tertiles_proliferation <- quantile(resMWW$TMB_Proliferation, probs = seq(0,1,1/3))
prol_low_cutoff <- tertiles_proliferation[[2]]
prol_high_cutoff <- tertiles_proliferation[[3]]

prediction_df <- NULL
for (i in 1:nrow(clinical_df))
{
  tmb_var <- clinical_df[i,]$TMB_Scaled
  pdl1_var <- clinical_df[i,]$PDL1
  
  #Convert the PDL1 score into category
  if (is.na(pdl1_var))
  {
    if (resMWW[i,]$ICR >= icr_cutoff)
    {
      pdl1_var2 <- "PDL1_High"
    }
  }
  else
  {
      if (pdl1_var>=50)
      {
        pdl1_var2 <- "PDL1_High"
      }
      else if (pdl1_var<=49)
      {
        pdl1_var2 <- "PDL1_Low"
      }
  }
  
  #Convert TMB score into category
  if (is.na(tmb_var))
  {
    if (resMWW[i,]$TMB_Proliferation>=prol_high_cutoff)
    {
      tmb_var2 <- "TMB_High"
    }
    else if (resMWW[i,]$TMB_Proliferation<=prol_low_cutoff)
    {
      tmb_var2 <- "TMB_Low"
    }
    else{
      tmb_var <- "TMB_Medium"
    }
  }
  else{
    if (tmb_var>=tmb_high_cutoff)
    {
      tmb_var2 <- "TMB_High"
    }
    else if (tmb_var<=tmb_low_cutoff)
    {
      tmb_var2 <- "TMB_Low"
    }
    else{
      tmb_var2 <- "TMB_Medium"
    }
  }
  temp <- cbind(clinical_df[i,]$patientID,tmb_var2,pdl1_var2)
  prediction_df <- rbind(prediction_df,temp)
}
prediction_df <- as.data.frame(prediction_df)
colnames(prediction_df) <- c("patientID","TMB_Signature","PDL1_Signature")
