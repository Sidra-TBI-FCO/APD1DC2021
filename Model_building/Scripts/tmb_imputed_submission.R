#!/usr/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(yaGST))
suppressMessages(library(kernlab))

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

#Get Proliferation matrix
prol_matrix <- normalized.log2.count[SelPath_Symb$TMB_Proliferation,]
rbf <- rbfdot(sigma=0.5)
sim_matrix <- kernelMatrix(rbf,t(prol_matrix))

#Get ids of all patients
get_all_patients <- colnames(normalized.log2.count)

#Creat a TMB based data frame
tmb_df <- NULL
for (i in 1:length(get_all_patients))
{
  patientID <- get_all_patients[i]
  if (patientID %in% clinical_df$patientID)
  {
    #If patient has tmb already don't do anything
    tmb_var <- clinical_df[clinical_df$patientID==patientID,]$TMB
    
    #If tmb is na then impute with weighted tmb of all patients with TMB present
    if (is.na(tmb_var))
    {
      get_weights <- sim_matrix[,patientID]
      temp <- 0
      for (j in 1:length(names(get_weights)))
      {
        t_pID <- names(get_weights)[j]
        t_TMB <- clinical_df[clinical_df$patientID==t_pID,]$TMB
        if (!is.na(t_TMB))
        {
          temp <- temp + t_TMB*as.numeric(get_weights[j])
        }
      }
      tmb_var <- temp/sum(get_weights)
    }
  }
  else{
    #Impute TMB as weighted tmb of all patients with tmb present based on distance
    get_weights <- sim_matrix[,patientID]
    temp <- 0
    for (j in 1:length(names(get_weights)))
    {
      t_pID <- names(get_weights)[j]
      t_TMB <- clinical_df[clinical_df$patientID==t_pID,]$TMB
      if (!is.na(t_TMB))
      {
        temp <- temp + t_TMB*as.numeric(get_weights[j])
      }
    }
    tmb_var <- temp/sum(get_weights)
  }
  output <- cbind(patientID,tmb_var)
  tmb_df <- rbind(tmb_df,output)
}
tmb_df <- as.data.frame(tmb_df)
colnames(tmb_df) <- c("patientID","TMB_Imputed")
tmb_df$patientID <- as.character(as.vector(tmb_df$patientID))
tmb_df$TMB_Imputed <- as.numeric(as.vector(tmb_df$TMB_Imputed))

tmb_df
