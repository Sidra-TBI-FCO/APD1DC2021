rm(list = ls())

load("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/synthetic_table_cluster.RData")
annotaion_file = read.csv("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/clinical_data.csv")

annotaion_file$ICRscore = clustering$ICRscore[match(annotaion_file$patientID,rownames(clustering))]
load("~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/Analysis/004.ssGSEASelected.pathways.Rdata")
#ES = ES[which(rownames(ES) == "[LM] Proliferation"), , drop = FALSE]
annotaion_file$Proliferation_score = ES["[LM] Proliferation",][match(annotaion_file$patientID , colnames(ES))]
annotaion_file$TGF_B_score = ES["[HM] TGF beta signaling",][match(annotaion_file$patientID , colnames(ES))]

write.csv(annotaion_file,file = "~/OneDrive - Hamad bin Khalifa University/Draem challenge/Synthetic_data/Variables.file",quote = FALSE,row.names = FALSE)
