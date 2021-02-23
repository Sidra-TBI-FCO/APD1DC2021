#!/usr/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(yaGST))

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/")

source("Model_building/Scripts/knn_graph.R")

get_imputed_tmb <- function(sim_matrix, clinical_df, patientID, k){
  
  #If tmb is na then impute with weighted tmb of all patients with TMB present
  get_weights <- sim_matrix[,patientID]
  get_weights <- sort(get_weights, decreasing = TRUE)
  temp <- 0
  count <- 0
  total_weights <- 0
  for (j in 1:length(names(get_weights)))
  {
    t_pID <- names(get_weights)[j]
    t_TMB <- clinical_df[clinical_df$patientID==t_pID,]$TMB
    if (!is.na(t_TMB))
    {
      #temp <- temp + t_TMB
      temp <- temp + t_TMB*as.numeric(get_weights[j])
      total_weights <- total_weights + as.numeric(get_weights[j])
      count <- count+1
    }
    if (count==k) break;
  }
  tmb_var <- temp/sum(total_weights)
  #tmb_var <- temp/k
  return(tmb_var)
}

normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}


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
clinical_df$TMB <- log2(clinical_df$TMB+1)

#Get the distance, gaussianized and similarity matrix
k <- 5
prol_matrix <- as.matrix(normalized.log2.count)
dist_matrix <- dist(t(prol_matrix),"euclidean")
W <- affinity_matrix(as.matrix(dist_matrix),k)
K <- W/rowSums(W)
diag(K) <- 0
sim_matrix <- K

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
    if (is.na(tmb_var))
    {
      tmb_var <- get_imputed_tmb(sim_matrix,clinical_df,patientID,k)
    }
  }
  else{
    tmb_var <- get_imputed_tmb(sim_matrix,clinical_df,patientID,k)
  }
  output <- cbind(patientID,tmb_var)
  tmb_df <- rbind(tmb_df,output)
}
tmb_df <- as.data.frame(tmb_df)
colnames(tmb_df) <- c("patientID","TMB_Imputed")
tmb_df$patientID <- as.character(as.vector(tmb_df$patientID))
tmb_df$TMB_Imputed <- as.numeric(as.vector(tmb_df$TMB_Imputed))

tmb_df$TMB_Imputed_Scaled <- normalize(tmb_df$TMB_Imputed)
tmb_df$ICR <- resMWW$ICR

final_df <- data.frame(patientID = tmb_df$patientID, prediction = tmb_df$TMB_Imputed_Scaled+0.2*max(tmb_df$ICR-0.5,0))
final_df$patientID <- as.character(as.vector(final_df$patientID))
final_df$prediction <- as.numeric(as.vector(final_df$prediction))

write.table(final_df, "output/predictions.csv", sep=",", quote=F, row.names=F, col.names=T)