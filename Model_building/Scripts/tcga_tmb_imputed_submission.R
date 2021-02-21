#!/usr/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(yaGST))
suppressMessages(library(kernlab))

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/")
load("Model_building/Required_Files/NSCLC_TMB_Proliferation.Rdata")

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
      temp <- temp + t_TMB
      #temp <- temp + t_TMB*as.numeric(get_weights[j])
      #total_weights <- total_weights + as.numeric(get_weights[j])
      count <- count+1
    }
    if (count==k) break;
  }
  #tmb_var <- temp/sum(total_weights)
  tmb_var <- temp/k
  return(tmb_var)
}

NSCLC_without_NA <- NSCLC[!is.na(NSCLC$TMB),]
orig_NSCLC_without_NA <- NSCLC_without_NA
N_without_NA <- nrow(NSCLC_without_NA)


seed_ids <- sample(1:1000,10)
val_cor <- NULL
full_cor <- NULL
full_median_cor <- NULL
for (seed_id in seed_ids)
{
  
  set.seed(seed_id)
  percent <- 0.2
  validation_sample_ids <- sample(N_without_NA)[1:ceiling(percent*N_without_NA)]
  
  #Get the validation samples
  validation_samples <- NSCLC_without_NA[validation_sample_ids,]
  
  #Replace validation TMB with NA
  NSCLC_without_NA[validation_sample_ids,]$TMB <- NA
  
  #Make similarity matrix
  k <- 15
  prol_matrix <- as.matrix(log2(NSCLC_without_NA[,c(1:5)]+1))
  dist_matrix <- dist(prol_matrix)
  W <- affinity_matrix(as.matrix(dist_matrix),k)
  K <- W/rowSums(W)
  diag(K) <- 0
  sim_matrix <- K
  
  #Get ids of all patients
  get_all_patients <- rownames(NSCLC_without_NA)
  NSCLC_without_NA$patientID <- rownames(NSCLC_without_NA)
  
  #Create a TMB based data frame
  tmb_df <- NULL
  median_tmb <- median(NSCLC_without_NA$TMB,na.rm = T)
  for (i in 1:length(get_all_patients))
  {
    patientID <- get_all_patients[i]
    if (patientID %in% NSCLC_without_NA$patientID)
    {
      #If patient has tmb already don't do anything
      tmb_var <- NSCLC_without_NA[NSCLC_without_NA$patientID==patientID,]$TMB
      if (is.na(tmb_var))
      {
        tmb_var <- get_imputed_tmb(sim_matrix,NSCLC_without_NA,patientID,k)
        median_tmb_var <- median_tmb
      }
      else{
        median_tmb_var <- tmb_var
      }
    }
    else{
      tmb_var <- get_imputed_tmb(sim_matrix,NSCLC_without_NA,patientID,k)
      median_tmb_var <- median_tmb
    }
    output <- cbind(patientID,tmb_var,median_tmb_var)
    tmb_df <- rbind(tmb_df,output)
  }
  tmb_df <- as.data.frame(tmb_df)
  colnames(tmb_df) <- c("patientID","TMB_Imputed","TMB_Imputed_Median")
  tmb_df$patientID <- as.character(as.vector(tmb_df$patientID))
  tmb_df$TMB_Imputed <- as.numeric(as.vector(tmb_df$TMB_Imputed))
  tmb_df$TMB_Imputed_Median <- as.numeric(as.vector(tmb_df$TMB_Imputed_Median))
  
  
  correlation_validation_set <- cor(tmb_df[validation_sample_ids,]$TMB_Imputed,validation_samples$TMB)
  correlation_full_set_imputed <- cor(tmb_df$TMB_Imputed,orig_NSCLC_without_NA$TMB)
  correlation_full_set_median_imputed <- cor(tmb_df$TMB_Imputed_Median, orig_NSCLC_without_NA$TMB)
  plot(validation_samples$TMB,tmb_df[validation_sample_ids,]$TMB_Imputed)
  plot(tmb_df$TMB_Imputed,orig_NSCLC_without_NA$TMB)
  plot(tmb_df$TMB_Imputed_Median, orig_NSCLC_without_NA$TMB)
  
  print(paste0("K = ",k," Percentage Missing is: ",percent*100,"% validation correlation = ",round(correlation_validation_set,3),
               " full correlation = ",round(correlation_full_set_imputed,3)," full correlation median = ",
               round(correlation_full_set_median_imputed,3)))
  
  
  val_cor <- c(val_cor,correlation_validation_set)
  full_cor <- c(full_cor,correlation_full_set_imputed)
  full_median_cor <- c(full_median_cor,correlation_full_set_median_imputed)
  NSCLC_without_NA <- orig_NSCLC_without_NA
}
boxplot(val_cor,ylab="Validation Correlation at 20% missingness",xlab="Imputed NAs")
boxplot(full_cor, full_median_cor, ylab="Full Correlation at 20% missingess", xlab="Left is Imputed and Right is Median")