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

#Get clinical data
clinical_df <- read.table("Synthetic_Data/clinical_data.csv",header=TRUE,sep=",")
clinical_df$patientID <- as.character(as.vector(clinical_df$patientID))

#Get ids of all patients
get_all_patients <- colnames(normalized.log2.count)

#Get TMB tertiles
all_tmb <- log2(clinical_df[!is.na(clinical_df$TMB),]$TMB+1)
clinical_df$TMB_Scaled <- log2(clinical_df$TMB+1)
tertiles_tmb <- quantile(all_tmb, probs = seq(0,1,1/3))
tmb_low_cutoff <- tertiles_tmb[[2]]
tmb_high_cutoff <- tertiles_tmb[[3]]

#Get ICR tertiles
#tertiles_icr <- quantile(resMWW$ICR, probs = seq(0,1,1/3))
tertiles_icr <- summary(resMWW$ICR)
icr_cutoff <- tertiles_icr[[3]]

#Get proliferation tertiles for imputing TMB
tertiles_proliferation <- quantile(resMWW$TMB_Proliferation, probs = seq(0,1,1/3))
prol_low_cutoff <- tertiles_proliferation[[2]]
prol_high_cutoff <- tertiles_proliferation[[3]]

prediction_df <- NULL
for (i in 1:length(get_all_patients))
{
  patient_id <- get_all_patients[i]
  tmb_var <- NA
  pdl1_var <- NA
  
  if (patient_id %in% clinical_df$patientID)
  {
    tmb_var <- clinical_df[clinical_df$patientID==patient_id,]$TMB_Scaled
    #pdl1_var <- clinical_df[clinical_df$patientID==patient_id,]$PDL1
  }
  
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
  else{
    if (tmb_var>=tmb_high_cutoff)
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

final_df2 <- NULL
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
  final_df2 <- rbind(final_df2,temp)
}
final_df2 <- as.data.frame(final_df2)
colnames(final_df2) <- c("patientID","prediction")
final_df2$patientID <- as.character(as.vector(final_df2$patientID))
final_df2$prediction <- as.numeric(as.vector(final_df2$prediction))

write.table(final_df2,"output/predictions.csv",row.names=F,col.names=T,quote=F,sep=",")