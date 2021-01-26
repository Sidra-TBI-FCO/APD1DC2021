#!/usr/bin/env Rscript

# Normalizing the raw data 
suppressMessages(library(data.table))
suppressMessages(library(caret))
suppressMessages(library(yaGST))
suppressMessages(library(Miracle))

#Load the log2 normalized dataset
load("Model_building/Required_Files/normalized-log2-count.Rdata")

#Load normalized data
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
load("Model_building/Required_Files/CV_ML_models.Rdata")
model_list <- list(gpFit, svmFit, rfFit, gbmFit)
names(model_list) <- c("GP","SVM","RF","GBM")
all_predictions <- predict(model_list, newdata = output_df[,req_columns], type="prob")

get_nonresponse_info <- NULL
for (i in 1:length(model_list))
{
  get_nonresponse_info <- cbind(get_nonresponse_info,all_predictions[[i]]$NonResponse)
}
get_nonresponse_info <- as.data.frame(get_nonresponse_info)

predictions_df <- data.frame("patientID"=rownames(output_df),"prediction"=rowMeans(get_nonresponse_info))
write.csv(predictions_df, file="output/predictions.csv",quote=F,row.names=F)

print("Done writing ML model predictions")
