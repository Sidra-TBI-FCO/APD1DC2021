
# Set working directory
setwd("C:/Users/eahmed-extern/Desktop/Dream Challenge")

# Load required packages
library(birk)
library(ggplot2)
library(ggpubr)
library(corrplot)

# Load required data
load("./nsclc/NSCLC TMB.Rdata")
load("./LUAD/TMB_LUAD_ssGSEA.Rdata")
load("./LUSC/TMB_LUAD_ssGSEA.Rdata")

# Combine LUAD and LUSC 
NSCLC_TMB = rbind(TMB_LUAD, TMB_LUSC)
NSCLC_TMB$patientID = substring(NSCLC_TMB$patientID, 1, 12)

# Log2 transformation
NSCLC$nonsilent.mutation.rate = log(NSCLC$nonsilent.mutation.rate + 1,2)
NSCLC = NSCLC[!is.na(NSCLC$nonsilent.mutation.rate),]

# Add TMB_Proliferation
NSCLC$TMB_Proliferation = NSCLC_TMB$TMB_Proliferation[match(NSCLC$bcr_patient_barcode, NSCLC_TMB$patientID)]
NSCLC = NSCLC[!is.na(NSCLC$TMB_Proliferation),]

# New object to introduce NAs 
NSCLC_NA = NSCLC

# Introduce random NAs (50%)
NSCLC_NA$bcr_patient_barcode <- as.character(as.vector(NSCLC_NA$bcr_patient_barcode))
set.seed(420)
ind <- which(NSCLC_NA$nonsilent.mutation.rate %in% sample(NSCLC_NA$nonsilent.mutation.rate, 200))  # patients with removed nas
NSCLC_NA$nonsilent.mutation.rate[ind]<-NA

# Create a numeric vector
vect = replace(NSCLC_NA$TMB_Proliferation,which(is.na(NSCLC_NA$nonsilent.mutation.rate)),NA)
# Replace NAs 
TMB_imputed = as.data.frame(transform(NSCLC_NA,TMB_imputed=nonsilent.mutation.rate[sapply(TMB_Proliferation,which.closest,vec=vect)]))

#ids = which(is.na(TMB_imputed$nonsilent.mutation.rate))

#correlation values
cor(NSCLC$nonsilent.mutation.rate[ind], TMB_imputed$TMB_imputed[ind])  # cor with NAs 
cor(NSCLC$nonsilent.mutation.rate, TMB_imputed$TMB_imputed) # cor without NAs

# new object to replace NAs with median 
NSCLC_median = NSCLC_NA

# replace missing NAs with median 
NSCLC_median$nonsilent.mutation.rate[is.na(NSCLC_median$nonsilent.mutation.rate)] <- median(NSCLC$nonsilent.mutation.rate, na.rm=TRUE)
cor(NSCLC$nonsilent.mutation.rate, NSCLC_median$nonsilent.mutation.rate)   

#df$value[is.na(df$value)] <- median(df$value, na.rm=TRUE)

# Medians
median(NSCLC$nonsilent.mutation.rate)
median(TMB_imputed$TMB_imputed)
median(NSCLC$nonsilent.mutation.rate[ind])
median(TMB_imputed$TMB_imputed[ind])
median(NSCLC_median$nonsilent.mutation.rate)

save(NSCLC, file = "./nsclc/TMB_Proliferation_NSCLC.Rdata")

################################################################################################################

# matrix with patients with imputed NAs only 
cor.na = cbind(NSCLC$nonsilent.mutation.rate[ind], TMB_imputed$TMB_imputed[ind])
cor.na = as.data.frame(cor.na)

ggscatter(cor.na, x = "V1", y = "V2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TMB", ylab = "TMB_imputed")


# matrix with all patients 
cor.im = cbind(NSCLC$nonsilent.mutation.rate, TMB_imputed$TMB_imputed)
cor.im = as.data.frame(cor.im)

ggscatter(cor.im, x = "V1", y = "V2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TMB", ylab = "TMB_imputed")

