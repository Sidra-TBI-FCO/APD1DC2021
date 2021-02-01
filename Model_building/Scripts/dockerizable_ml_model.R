#!/usr/bin/env Rscript

# Normalizing the raw data 
suppressMessages(library(data.table))
suppressMessages(library(caret))
suppressMessages(library(yaGST))
suppressMessages(library(Miracle))
suppressMessages(library(xgboost))

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/")

#Load the log2 normalized dataset
#load("/data1/normalized-log2-count.RData")
load("Model_building/Required_Files/normalized-log2-count.Rdata")

#Load normalized data
#load("/data1/Selelected_Path_VariousGeneIDs.RData")
load("Model_building/Required_Files/Selelected_Path_VariousGeneIDs.RData")

#Perform Miracle
Mir_res_ALL <- Calculate_Miracle(normalized.log2.count, platform = "gene")  #available platforms: ens", "u133p2", "entrez", "gene"

## Calculate enrichment scores using gene Symbols. Note that you can use other gene annotations
resMWW <- c()
for(i in 1:ncol(normalized.log2.count)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(normalized.log2.count[,i], x)$nes )))
}
resMWW

output_df <- cbind(Mir_res_ALL,resMWW)
output_df <- as.data.frame(output_df)
for (i in 1:ncol(output_df))
{
  output_df[,i]  <- as.numeric(as.vector(output_df[,i]))
}

#Load the best models 
#load("/data1/CV_ML_models.Rdata")
load("Model_building/Required_Files/CV_ML_models.Rdata")
model_list <- list(gpFit, rdaFit, rfFit, gbmFit, xgbFit)
names(model_list) <- c("GP","RDA","RF","GBM","XGB")
all_predictions <- predict(model_list, newdata = output_df[,req_columns], type="prob")

get_response_info <- NULL
for (i in 1:length(model_list))
{
  get_response_info <- cbind(get_response_info,all_predictions[[i]]$Response)
}
get_response_info <- as.data.frame(get_response_info)
predictions_df <- data.frame("patientID"=rownames(output_df),"prediction"=rowMeans(get_response_info))

#Get clinical data
clinical_df <- read.table("Synthetic_Data/clinical_data.csv",header=TRUE,sep=",")
clinical_df$patientID <- as.character(as.vector(clinical_df$patientID))

par(mfrow=c(2,1))
plot(clinical_df$TMB)
plot(1/(1+exp(-log(clinical_df$TMB)+log(as.numeric(quantile(clinical_df$TMB, na.rm=T)[3])))))
clinical_df$TMB_Scaled <- 1/(1+exp(-log(clinical_df$TMB)+log(as.numeric(quantile(clinical_df$TMB, na.rm=T)[3]))))

final_df <- NULL
for (i in 1:nrow(predictions_df))
{
  patientID <- as.character(predictions_df[i,]$patientID)
  tmb_var <- clinical_df[clinical_df$patientID==patientID,]$TMB_Scaled
  if (!is.na(tmb_var))
  {
    prediction <- (1*predictions_df[i,]$prediction+1*tmb_var)/2
  }
  else{
    prediction <- predictions_df[i,]$prediction
  }
  final_df <- rbind(final_df,cbind(patientID,prediction))
}
final_df <- as.data.frame(final_df)
colnames(final_df) <- c("patientID","prediction")
final_df$patientID <- as.character(as.vector(final_df$patientID))
final_df$prediction <- as.numeric(as.vector(final_df$prediction))

write.csv(final_df, file="output/mix_predictions.csv",quote=F,row.names=F)

print("Done writing ML model predictions")
