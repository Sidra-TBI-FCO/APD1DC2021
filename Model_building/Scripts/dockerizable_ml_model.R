#!/usr/bin/env Rscript

# Normalizing the raw data 
suppressMessages(library(data.table))
suppressMessages(library(caret))
suppressMessages(library(yaGST))
suppressMessages(library(Miracle))
suppressMessages(library(xgboost))
suppressMessages(library(matrixStats))

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/")

#Load the log2 normalized dataset
load("Model_building/Required_Files/revised_normalized-log2-count.RData")

#Load normalized data
#load("/data1/Selelected_Path_VariousGeneIDs.RData")
load("Model_building/Required_Files/Selelected_Path_VariousGeneIDs.RData")
load("Model_building/Required_Files/c2 Lance Signatures.Rdata")

#Perform Miracle
Mir_res_ALL <- Calculate_Miracle(normalized.log2.count, platform = "gene")  #available platforms: ens", "u133p2", "entrez", "gene"

#Add Lance Signature to SelPath_Symb
SelPath_Symb$Lance_Signature <- c2.signatures

## Calculate enrichment scores using gene Symbols. Note that you can use other gene annotations
resMWW <- c()
for(i in 1:ncol(normalized.log2.count)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(normalized.log2.count[,i], x)$nes )))
  LanceScore <- mwwExtGST(as.matrix(normalized.log2.count)[,i], SelPath_Symb$Lance_Signature[1:10],
                          SelPath_Symb$Lance_Signature[11:12])
  resMWW[nrow(resMWW),"Lance_Signature"] <- LanceScore$nes
}
resMWW

output_df <- cbind(Mir_res_ALL,resMWW)
output_df <- as.data.frame(output_df)
for (i in 1:ncol(output_df))
{
  output_df[,i]  <- as.numeric(as.vector(output_df[,i]))
}

#Load the best models 
load("Model_building/Required_Files/CV_ML_models_submission_q3_v3.Rdata")
model_list <- list(gpFit, rfFit, rdaFit, gbmFit, svmFit)
names(model_list) <- c("GP","RF","RDA","GBM","SVM")
all_predictions <- predict(model_list, newdata = output_df[,req_columns], type="prob")

get_response_info <- NULL
for (i in 1:length(model_list))
{
  get_response_info <- cbind(get_response_info,all_predictions[[i]]$NonResponse)
}
get_response_info <- as.data.frame(get_response_info)
predictions_df <- data.frame("patientID"=rownames(output_df),"prediction"=rowMedians(as.matrix(get_response_info)))

#Get clinical data
clinical_df <- read.table("Synthetic_Data/clinical_data.csv",header=TRUE,sep=",")
clinical_df$patientID <- as.character(as.vector(clinical_df$patientID))

par(mfrow=c(2,1))
plot(clinical_df$TMB)
plot((1/(1+exp(-log(clinical_df$TMB)+log(243)))))
#clinical_df$TMB_Scaled <- 1/(1+exp(-log(clinical_df$TMB)+log(as.numeric(quantile(clinical_df$TMB, na.rm=T)[3]))))
clinical_df$TMB_Scaled <- (1/(1+exp(-log(clinical_df$TMB)+log(243))))

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

library(ROCR)
library(pROC)
library(PRROC)

load("Model_building/Required_Files/clinical_HNSC.Rdata")
load("Model_building/Required_Files/predictions_df_HNSC.Rdata")
clinical$Cat_Response <- 1
clinical[clinical$simple.response!="Responder",]$Cat_Response <- 0

pred <- prediction(predictions_df$prediction, clinical$Cat_Response)
perf <- performance(pred, "prec", "rec")
plot(perf, avg= "threshold", colorize=TRUE, lwd= 3, ylim = c(0,1),
     main= "... Precision/Recall graphs ...")

fg <- predictions_df[which(clinical$Cat_Response==0),]$prediction
bg <- predictions_df[which(clinical$Cat_Response==1),]$prediction
pr <- pr.curve(scores.class0 = fg , scores.class1 = bg, curve = T)
roc <- roc.curve(scores.class0 = fg , scores.class1 = bg, curve=T)
plot(pr)
