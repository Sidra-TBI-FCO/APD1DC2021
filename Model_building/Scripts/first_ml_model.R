library(data.table)
library(caret)
library(AppliedPredictiveModeling)
library(MLmetrics)
library(doParallel)
library(gbm)
library(randomForest)
library(kernlab)
library(bst)
library(xgboost)

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/")

#Set no of parallel cpus
cl <- makePSOCKcluster(16)
registerDoParallel(cl)

#Load all the data
load("Model_building/Required_Files/Master_Datasets_Selected_response.RData")

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

#Datasets with response variable
dataset_ids <- NULL
col_ids <- NULL
for (i in 1:length(dfs))
{
  temp_df <- dfs[[i]]
  colnames_temp_df <- colnames(temp_df)
  if ("Response" %in% colnames_temp_df ||
      "Status" %in% colnames_temp_df || "ResponseOverall" %in% colnames_temp_df)
  {
    dataset_ids <- c(dataset_ids,i)
    col_id <- c(which(colnames_temp_df=="Response"),which(colnames_temp_df=="Status"),which(colnames_temp_df=="ResponseOverall"))
    col_ids <- c(col_ids, col_id)
  }
}

# #Make a uniform "Status" column in same scale in all datasets
# for (i in 1:length(dataset_ids))
# {
#   dataset_id <- dataset_ids[i]
#   temp_df <- dfs[[dataset_id]]
#   col_id <- col_ids[i]
#   response_variable <- colnames(temp_df)[col_id]
#   if (response_variable == "ResponseOverall")
#   {
#     temp_df$Status <- -1
#     temp_df[temp_df[,col_id]=="R",]$Status <- 1
#     temp_df[temp_df[,col_id]=="NR",]$Status <- 0
#   }
#   else if (response_variable == "Response")
#   {
#     if (class(temp_df[,col_id])=="character" || class(temp_df[,col_id])=="factor")
#     {
#       temp_df$Response <- as.character(as.vector(temp_df$Response))
#       temp_df$Status <- -1
#       temp_df[temp_df[,col_id]=="Response",]$Status <- 1
#       temp_df[temp_df[,col_id]=="Nonresponse",]$Status <- 0
#     }
#     else
#     {
#       temp_df$Status <- -1
#       temp_df[temp_df[,col_id]==0,]$Status <- 0
#       temp_df[temp_df[,col_id]==1,]$Status <- 1
#     }
#   }
#   dfs[[dataset_id]] <- temp_df
# }

#Get all the common colnames (remove Chen dataset)
common_colnames <- colnames(dfs[[dataset_ids[[2]]]])
for (i in 3:length(dataset_ids))
{
#  if (i!=4)
#  {
    common_colnames <- intersect(common_colnames,colnames(dfs[[dataset_ids[i]]]))
#  }
}

#Make the combined df and remove rows with NA in Status (remove Chen)
combined_df <- NULL
for (i in 1:length(dataset_ids))
{
  if (i!=1)
  {
    temp_df <- dfs[[dataset_ids[i]]]
    temp_df <- temp_df[,common_colnames]
    rev_temp_df <- temp_df[which(!is.na(temp_df$response1)),]
    combined_df <- rbind(combined_df,rev_temp_df)
  }
}
combined_df <- as.data.frame(combined_df)
colnames(combined_df) <- common_colnames
combined_df$response1 <- as.numeric(as.vector(combined_df$response1))
combined_df$response1 <- as.factor(as.vector(combined_df$response1))

#Make Scatterplot with "ICR", "IE_Specific" , "Miracle", "ID_Specific"
transparentTheme(trans = .4)
featurePlot(x = combined_df[, c(1:4,58)], 
            y = combined_df$response1, 
            plot = "ellipse",
            ## Add a key at the top
            auto.key = list(columns = 2))


#Make the density plots
transparentTheme(trans = .9)
featurePlot(x = combined_df[, c(1:4,58)], 
            y = combined_df$response1,
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(5, 1), 
            auto.key = list(columns = 2))

#Make the boxplots
featurePlot(x = combined_df[, c(1:4,58)], 
            y = combined_df$response1, 
            plot = "box", 
            ## Pass in options to bwplot() 
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)),  
            layout = c(5,1 ), 
            auto.key = list(columns = 2))


#Check for nonzero variance 
nzv <- nearZeroVar(combined_df, saveMetrics= TRUE)

#Find correlated descriptors
descrCor <-  cor(combined_df[,c(1:58)])
highlyCorDescr <- findCorrelation(descrCor, cutoff = 0.99)
summary(descrCor[upper.tri(descrCor)])
combined_df$response1 <- as.numeric(as.vector(combined_df$response1))

#Revised Data Frame after removal of highly correlated features
rev_df <- combined_df[,-highlyCorDescr]
req_columns <- setdiff(colnames(rev_df),"response1")
rev_df$response1 <- as.factor(as.vector(rev_df$response1))
levels(rev_df$response1) <- c("Response","NonResponse")

#Construct a simple ML model
set.seed(998)
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           #search = "random",
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = prSummary)


gbmGrid <-  expand.grid(interaction.depth = c(1,2,3), 
                        n.trees = (1:30)*50, 
                        shrinkage = c(0.01,0.05,0.1),
                        n.minobsinnode = c(10,15))

gbmFit <- train(response1 ~ ., data = rev_df, 
                 method = "gbm", 
                 trControl = fitControl,
                 metric = "AUC",
                 verbose = FALSE,
                 tuneGrid = gbmGrid
                 )
gbmFit

#Make plot of the hyper-parameters
trellis.par.set(caretTheme())
plot(gbmFit)

#Best model id and parameters and performance
best_id <- best(gbmFit$results, metric="AUC", maximize = T)
gbmFit$results[best_id,]

#Get variable importance for GBM model
varImp(gbmFit)

###########################################################################################################
#Make hyper-parameter grid for SVM
svmGrid <-  expand.grid(C = 2^seq(-5,5),
                        sigma=2^seq(-20,2))

#Make an SVM model
svmFit <- train(response1 ~ ., data = rev_df, 
                method = "svmRadial", 
                trControl = fitControl, 
                preProc = c("center", "scale"),
                tuneGrid = svmGrid,
                metric = "AUC")
svmFit                 

#Make plot of the hyper-parameters
trellis.par.set(caretTheme())
plot(svmFit)

#Get variable importance for SVM model
varImp(svmFit)

############################################################################################################
#Make a Regularized Discriminant Analysis model
rdaFit <- train(response1 ~ ., data = rev_df, 
                method = "rda", 
                trControl = fitControl, 
                preProc = c("center","scale"),
                tuneLength = 10,
                metric = "AUC")
rdaFit    

#Make plot of hyper-parameters
trellis.par.set(caretTheme())
plot(rdaFit)

#Get variable importance for RDA model
varImp(rdaFit)

###########################################################################################################
#Fit a Gaussian process 
gpGrid <- expand.grid(.sigma = 2^seq(-20,2))

gpFit <- train(response1 ~ ., data = rev_df,
               method = "gaussprRadial",
               trControl = fitControl,
               preProc = c("center", "scale"),
               tuneGrid = gpGrid,
               metric = "AUC")

gpFit

#Make plot of hyper-parameters
trellis.par.set(caretTheme())
plot(gpFit)

#Get variable importance from GP
varImp(gpFit)

##########################################################################################################
#Fit a Random Forest model
rfGrid <- expand.grid(mtry = c(1,2,3,4))
                      
rfFit <- train(response1 ~ ., data = rev_df,
               method = "parRF",
               trControl = fitControl,
               tuneGrid = rfGrid,
               metric = "AUC"
              )
rfFit

#Make plot of hyper-parameters
trellis.par.set(caretTheme())
plot(rfFit)

#Get variable importance for RF
varImp(rfFit)

###########################################################################################################
#Fit a boosted tree
xgbGrid <-  expand.grid(max_depth = c(1,2,3,4,5), 
                        nrounds = c(100,500,1000),
                        eta = c(0.1,0.2,0.3),
                        colsample_bytree = c(0.1,0.3,0.5,1.0),
                        gamma = c(0),
                        subsample = c(1.0),
                        min_child_weight = c(1,5,10,20)
                       )

xgbFit <- train(response1 ~ ., data = rev_df, 
                method = "xgbTree", 
                trControl = fitControl,
                metric = "AUC",
                verbose = TRUE,
                tuneGrid = xgbGrid,
                nthread = 1
                )
xgbFit

#Make plot of hyper-parameters
trellis.par.set(caretTheme())
plot(xgbFit)

#Get variable importance for XGB
varImp(xgbFit)

###################################################################################################################
#Make comparison of the models
resamps <- resamples(list(GBM = gbmFit,
                          SVM = svmFit,
                          RDA = rdaFit,
                          GP = gpFit,
                          RF = rfFit,
                          XGB = xgbFit))
resamps

#Make trellis plot
bwplot(resamps, layout = c(2, 2))

#Make scatterplot based comparison
splom(resamps)

#Stop parallel processing
stopCluster(cl)

save(list=c("req_columns","gbmFit","svmFit","rdaFit","gpFit","rfFit","xgbFit"),file="Model_building/Required_Files/CV_ML_models_submission_q3_v2.Rdata")
