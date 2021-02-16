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
resMWW$patientID <- colnames(normalized.log2.count)

#Get ids of all patients
get_all_patients <- colnames(normalized.log2.count)

#Get ICR tertiles
#tertiles_icr <- quantile(resMWW$ICR, probs = seq(0,1,1/3))
tertiles_icr <- summary(resMWW$ICR)
icr_cutoff <- tertiles_icr[[4]]

#Get proliferation tertiles for imputing TMB
tertiles_proliferation <- quantile(resMWW$TMB_Proliferation, probs = seq(0,1,1/3))
prol_high_cutoff <- tertiles_proliferation[[3]]

prediction_df <- NULL
for (i in 1:length(get_all_patients))
{
  patient_id <- get_all_patients[i]
  tmb_var <- NA
  pdl1_var <- NA
  
  #Convert the PDL1 score into category
  resMMW_patient_idx <- which(resMWW$patientID==patient_id)
  if (is.na(pdl1_var))
  {
    if (resMWW[resMMW_patient_idx,]$ICR >= icr_cutoff)
    {
      pdl1_var2 <- "PDL1_High"
    }
    else{
      pdl1_var2 <- "PDL1_Low"
    }
  }
  
  #Convert TMB score into category
  if (is.na(tmb_var))
  {
    if (resMWW[resMMW_patient_idx,]$TMB_Proliferation>=prol_high_cutoff)
    {
      tmb_var2 <- "TMB_High"
    }
    else 
    {
      tmb_var2 <- "TMB_Low"
    }
  }
  temp <- cbind(patient_id,tmb_var2,pdl1_var2)
  prediction_df <- rbind(prediction_df,temp)
}
prediction_df <- as.data.frame(prediction_df)
colnames(prediction_df) <- c("patientID","TMB_Signature","PDL1_Signature")
prediction_df$patientID <- as.character(as.vector(prediction_df$patientID))
prediction_df$TMB_Signature <- as.character(as.vector(prediction_df$TMB_Signature))
prediction_df$PDL1_Signature <- as.character(as.vector(prediction_df$PDL1_Signature))

final_df <- NULL
for (i in 1:nrow(prediction_df))
{
  if (prediction_df[i,]$TMB_Signature=="TMB_High" & prediction_df[i,]$PDL1_Signature=="PDL1_High")
  {
    score <- 3
  }
  else if (prediction_df[i,]$TMB_Signature=="TMB_Low" & prediction_df[i,]$PDL1_Signature=="PDL1_Low")
  {
    score <- 1
  }
  else{
    score <- 2 
  }
  temp <- cbind(prediction_df[i,]$patientID,score)
  final_df <- rbind(final_df,temp)
}
final_df <- as.data.frame(final_df)
colnames(final_df) <- c("patientID","prediction")
final_df$patientID <- as.character(as.vector(final_df$patientID))
final_df$prediction <- as.numeric(as.vector(final_df$prediction))

write.table(final_df,"output/predictions.csv",row.names=F,col.names=T,quote=F,sep=",")