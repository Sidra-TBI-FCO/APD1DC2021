# this is the final prediction script
library(scales)
#load data
load ("/data1/patway_score.RData")
ES <- t(ES)
ES <- as.data.frame(ES)
clinical_data <- read.csv("/data1/clinical_data.csv")
#clinical_data$patientID == rownames(ES)

#patients without TMB
no.TMB.index <- which(is.na(clinical_data$TMB))

#make table
table <- data.frame(patientID = clinical_data$patientID,
                    TMB = clinical_data$TMB)

#scale the ssGSEA scores
table$proliferation.scaled <- round(rescale(ES$`[LM] Proliferation`,to = c(1,20)))
table$proliferation.scaled_50 <- round(rescale(ES$`[LM] Proliferation`,to = c(1,50)))
table$ICR.scaled <- round(rescale(ES$`[TBI] ICR`,to = c(1,20)))
table$TGFB.scaled <- round(rescale(ES$`[HM] TGF beta signaling`,to = c(1,20)))
table$ICR_TGFB_scaled <-round(rescale(table$ICR.scaled-table$TGFB.scaled,to = c(1,20)))
table$ICR_TGFB_scaled_50 <-round(rescale(table$ICR.scaled-table$TGFB.scaled,to = c(1,50)))
#TMB.scaled <- round(rescale(clinical_data$TMB,to = c(1,60)))
table$TMB.scaled <- round(rescale(log(table$TMB,10),to = c(1,60)))

#sum the scaled scores
table$prediction.score <- table$ICR_TGFB_scaled+table$TMB.scaled+table$proliferation.scaled
table$prediction.score_50 <- table$ICR_TGFB_scaled_50+table$proliferation.scaled_50

summary(table$prediction.score)
summary(table$prediction.score_50)

#overwrite NA with 50/50
table$prediction.score_final <- table$prediction.score
table$prediction.score_final[no.TMB.index] <- table$prediction.score_50[no.TMB.index]
summary(table$prediction.score)

#save output
clean.table = table[,c("patientID","prediction.score_final")]
colnames(clean.table) = c("patientID","prediction")
dir.create("/output",showWarnings = F)
write.csv(clean.table,file = "/output/predictions.csv",row.names = F,quote=F)
#write.csv(table,file = "/output/prediction_table.csv",row.names = F)
