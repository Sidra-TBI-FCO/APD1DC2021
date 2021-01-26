# library(devtools)
# install_github("miccec/yaGST")
# install_github("tolgaturan-github/Miracle")


setwd("Anti-PD1 Dream Challenge Team/")
library(yaGST)
library(Miracle)
library(GSVA)
# ?Calculate_Miracle

# the dataset are in the google drive
expr0 <- readRDS("Response Datasets Michele/ICI response datasets/GSE126044_Normalized_Expression.rds")
meta0 <- readRDS("Response Datasets Michele/ICI response datasets/GSE126044_Meta_data.rds")

expr <- readRDS("Response Datasets Michele/ICI response datasets/Normalized_expression_Melanoma_Response (1).rds")
meta <- readRDS("Response Datasets Michele/ICI response datasets/Meta_data_Melanoma_Response (1).rds")

#### Pathways
# load("GM_AntiPD1/Data/Selected.pathways.3.4.RData")
# Selected.pathways[[55]] <- c("NUF2", "NEK2", "TPX2", "KIF2C", "MCM10")
# names(Selected.pathways)[[55]] <- "TMB_Proliferation"
# library("biomaRt")
# SelPath_entrex <- list()
# SelPath_ensembl <- list()
# SelPath_affy <- list()
# SelPath_ilm <- list()
# for(i in 1:length(Selected.pathways)){
#   ensembl_human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#   geneID_human = getBM(attributes = c('ensembl_gene_id','hgnc_symbol','entrezgene_id', 'affy_hg_u133_plus_2','illumina_humanht_12_v4'),
#                        filters =  'hgnc_symbol',
#                        values=  Selected.pathways[[i]],
#                        mart = ensembl_human)
#   entr <- unique(geneID_human$entrezgene_id)
#   SelPath_entrex <- c(SelPath_entrex, list(entr[entr != ""] ))
#   ens <- unique(geneID_human$ensembl_gene_id)
#   SelPath_ensembl <- c(SelPath_ensembl, list(ens[ens != ""]))
#   aff <- unique(geneID_human$affy_hg_u133_plus_2)
#   SelPath_affy <- c(SelPath_affy, list(aff[aff != ""]))
#   ill <- unique(geneID_human$illumina_humanht_12_v4)
#   SelPath_ilm <- c(SelPath_ilm, list(ill[ill != ""]))
# }
# names(SelPath_entrex) <- names(Selected.pathways)
# names(SelPath_ensembl) <- names(Selected.pathways)
# names(SelPath_affy) <- names(Selected.pathways)
# names(SelPath_ilm) <- names(Selected.pathways)
# SelPath_Symb <- Selected.pathways
# save(SelPath_Symb, SelPath_entrex, SelPath_ensembl, SelPath_affy, SelPath_ilm,
#      file="GM_AntiPD1/Data/SelPath.RData")

load("GM_AntiPD1/Data/SelPath.RData")

### GSE126044
Mir0 <- Calculate_Miracle(expr0, platform = "gene")
# ES = t(gsva(expr0,SelPath_Symb,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(expr0)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(expr0[,i], x)$nes )))
}
Cho_GSE126044 <- cbind(Mir0, meta0[ rownames(Mir0),], resMWW)
Cho_GSE126044$Response <- factor(Cho_GSE126044$Response, levels = c("Nonresponse", "Response"))


### Riaz GSE91061
Riaz <- Calculate_Miracle(expr$Riaz)
RiazMeta <- cbind(meta$Riaz_Melanoma_All, Cohort=NA,  Treat=NA)
RiazMeta[rownames(meta$Riaz_Melanoma_NAIVE), "Cohort"] <- "NAIVE"
RiazMeta[rownames(meta$Riaz_Melanoma_PROG), "Cohort"] <- "PROG"
RiazMeta[rownames(meta$Riaz_Melanoma_On), "Treat"] <- "On"
RiazMeta[rownames(meta$Riaz_Melanoma_Pre), "Treat"] <- "Pre"

# ES = t(gsva(expr$Riaz,SelPath_ensembl,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(expr$Riaz)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_ensembl,function(x) mwwGST(expr$Riaz[,i], x)$nes )))
}
Riaz_GSE91061 <- cbind(Riaz, RiazMeta[rownames(Riaz), ], resMWW)

#### PRJEB23709 Gide
Gide <- Calculate_Miracle(expr$Gide)
GideMeta <- cbind(meta$Gide_Melanoma_All, Therapy=NA,  Treat=NA)
GideMeta[rownames(meta$Gide_Melanoma_Combo_OnTreatment), "Therapy"] <- "Combo"
GideMeta[rownames(meta$Gide_Melanoma_Combo_Pooled), "Therapy"] <- "Combo"
GideMeta[rownames(meta$Gide_Melanoma_Combo_Pretreatment), "Therapy"] <- "Combo"
GideMeta[rownames(meta$Gide_Melanoma_PD1_OnTreatment), "Therapy"] <- "PD1"
GideMeta[rownames(meta$Gide_Melanoma_PD1_Pretreatment), "Therapy"] <- "PD1"
GideMeta[rownames(meta$Gide_Melanoma_PD1_Pooled), "Therapy"] <- "PD1"
GideMeta[rownames(meta$Gide_Melanoma_Combo_Pretreatment), "Treat"] <- "Pre"
GideMeta[rownames(meta$Gide_Melanoma_PD1_Pretreatment), "Treat"] <- "Pre"
GideMeta[rownames(meta$Gide_Melanoma_Combo_OnTreatment), "Treat"] <- "On"
GideMeta[rownames(meta$Gide_Melanoma_PD1_OnTreatment), "Treat"] <- "On"

# ES = t(gsva(expr$Gide,SelPath_ensembl,method="ssgsea"))

resMWW <- c()
for(i in 1:ncol(expr$Gide)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_ensembl,function(x) mwwGST(expr$Gide[,i], x)$nes )))
}

Gide_PRJEB23709 <- cbind(Gide, GideMeta[rownames(Gide), ], resMWW)


#### Van Allen
vanAllen_res <- Calculate_Miracle(expr$vanAllen, platform = "entrez")
vanMeta <- meta$`VanAllenMelanoma_anti-CTLA4` 

# ES = t(gsva(expr$vanAllen,SelPath_entrex,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(expr$vanAllen)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_entrex,function(x) mwwGST(expr$vanAllen[,i], x)$nes )))
}

vanAllen <- cbind(vanAllen_res, vanMeta[rownames(vanAllen_res), ], resMWW)

# ggplot(VanAllen_Miracle_Response, aes(x=Status, y=Miracle, colour=Status))+geom_boxplot(outlier.shape=NA)+geom_jitter(shape=16, position=position_jitter(0.2))



#### Ribas GSE78220
Ribas <- Calculate_Miracle(expr$Ribas)
RibasMeta <- cbind(meta$GSE78220_Melanoma_All, mapki =NA)
RibasMeta[rownames(meta$GSE78220_Melanoma_no_previous_mapki), "mapki"] <- "No"
RibasMeta[rownames(meta$GSE78220_Melanoma_previous_mapki), "mapki"] <- "Yes"

# ES = t(gsva(expr$Ribas,SelPath_ensembl,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(expr$Ribas)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_ensembl,function(x) mwwGST(expr$Ribas[,i], x)$nes )))
}

Ribas_GSE78220 <- cbind(Ribas, RibasMeta[rownames(Ribas), ], resMWW)


#### Dizier GSE35640

Dizier <- Calculate_Miracle(expr$Dizier, platform= "u133p2")
rownames(Dizier) <- substr(rownames(Dizier), 1,9)
DizierMeta <- meta$GSE35640_Melanoma_MAGEA3_IO

# ES = t(gsva(expr$Dizier,SelPath_affy,method="ssgsea"))
# rownames(ES) <- substr(rownames(ES), 1,9)

resMWW <- c()
for(i in 1:ncol(expr$Dizier)){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_affy,function(x) mwwGST(expr$Dizier[,i], x)$nes )))
}

UlloaMontoya_GSE35640 <- cbind(Dizier, DizierMeta[rownames(Dizier), ], resMWW)



#### Chen
# chenExpr <- read.csv("Training Data/chen et at/gene counts_symbol.csv", row.names = 1)
# chenExpr$symbol <- NULL 
# Chen_res <- Calculate_Miracle(as.matrix(chenExpr), platform= "entrez") 
# Note that I also tried Miracle with Gene Symbols, but there are even less genes available

# ChenMeta <- read.delim("Training Data/chen et at/clinical_data.txt", row.names = 1)

# ES = t(gsva(as.matrix(chenExpr),SelPath_entrex,method="ssgsea"))
# resMWW <- c()
# for(i in 1:ncol(as.matrix(chenExpr))){
#   resMWW <- rbind(resMWW, unlist(lapply(SelPath_entrex,function(x) mwwGST(as.matrix(chenExpr)[,i], x)$nes )))
# }
# Chen <- cbind(Chen_res, Status= ChenMeta[rownames(Chen_res), ], ES)


ChenExpr <- read.delim("Training Data/chen et at/dataChen.txt", row.names = 1,header = T, check.names = F)
colnames(ChenExpr) <- paste0(colnames(ChenExpr),1:54)
Chen_res <- Calculate_Miracle(as.matrix(ChenExpr), platform = "gene") 

ChenMeta <- read.delim("Training Data/chen et at/metadataChen.txt",  header = T, check.names = F)
ChenMeta$ResponseOverall <- ChenMeta$`anti-PD-1 response`
ChenMeta$ResponseOverall[which(ChenMeta$ResponseOverall == "N/A")] <- "R"
resMWW <- c()
for(i in 1:ncol(as.matrix(ChenExpr))){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(as.matrix(ChenExpr)[,i], x)$nes )))
}
Chen <- cbind(Chen_res, ChenMeta, resMWW)



#### GSE40419 Jang - Lung

JangExpr0 <- read.delim("Training Data/GSE40419 - lung/GSE40419_LC-87_RPKM_expression.txt")
JangExpr <- JangExpr0[!duplicated(JangExpr0$gene), ]
rownames(JangExpr) <-JangExpr$gene
JangExpr <- JangExpr[,grep("LC_", colnames(JangExpr))]
Jang_res <- Calculate_Miracle(as.matrix(JangExpr), platform= "gene") 

JangMeta <- read.csv("Training Data/GSE40419 - lung/clinical_data.csv", row.names = 2)

# ES = t(gsva(as.matrix(JangExpr),SelPath_Symb,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(as.matrix(JangExpr))){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(as.matrix(JangExpr)[,i], x)$nes )))
}
Jang_GSE40419 <- cbind(Jang_res, JangMeta[rownames(Jang_res), ], resMWW)

### GSE52562  NOT DONE YET

# there is no response to treatment, but there is survival


### GSE67501 Ascierto
AsciertoExpr0 <- read.csv("Training Data/GSE67501_Ascierto/GSE67501_non_normalized_read.counts.csv")
# AsciertoExpr <- AsciertoExpr0[,grep("RCC", colnames(AsciertoExpr0))]
# rownames(AsciertoExpr) <- AsciertoExpr0$ID_REF
AsciertoExpr <- AsciertoExpr0[!duplicated(AsciertoExpr0$gene.sympol), ]
AsciertoExpr <- AsciertoExpr[,grep("RCC", colnames(AsciertoExpr))]
rownames(AsciertoExpr) <- AsciertoExpr0$gene.sympol[!duplicated(AsciertoExpr0$gene.sympol)]


Ascierto_res <- Calculate_Miracle(as.matrix(AsciertoExpr), platform= "gene") 

AsciertoMeta <- read.csv("Training Data/GSE67501_Ascierto/clinical_data.csv", row.names = 2)
AsciertoMeta$Sample_description <- gsub("-", ".",AsciertoMeta$Sample_description)
rownames(AsciertoMeta) <- AsciertoMeta$Sample_description

# ES = t(gsva(as.matrix(AsciertoExpr),SelPath_Symb,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(as.matrix(AsciertoExpr))){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(as.matrix(AsciertoExpr)[,i], x)$nes )))
}

Ascierto_GSE67501 <- cbind(Ascierto_res, AsciertoMeta[rownames(Ascierto_res), ], resMWW)



### GSE79691 Ascierto
AsciertoExpr0 <- read.csv("Training Data/GSE79691_Ascierto/GSE79691_non_normalized_read.counts.csv")
AsciertoExpr <- AsciertoExpr0[!duplicated(AsciertoExpr0$gene.sympol), ]
AsciertoExpr <- AsciertoExpr[,grep("^M", colnames(AsciertoExpr))]
rownames(AsciertoExpr) <- AsciertoExpr0$gene.sympol[!duplicated(AsciertoExpr0$gene.sympol)]

Ascierto_res <- Calculate_Miracle(as.matrix(AsciertoExpr), platform= "gene") 

AsciertoMeta <- read.csv("Training Data/GSE79691_Ascierto/clinical_data.csv")
rownames(AsciertoMeta) <- AsciertoMeta$Sample_description

# ES = t(gsva(as.matrix(AsciertoExpr),SelPath_Symb,method="ssgsea"))

resMWW <- c()
for(i in 1:ncol(as.matrix(AsciertoExpr))){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(as.matrix(AsciertoExpr)[,i], x)$nes )))
}

Ascierto_GSE79691 <- cbind(Ascierto_res, AsciertoMeta[rownames(Ascierto_res), ], resMWW)


### GSE115821 Auslander

Auslander0 <- read.csv("Training Data/GSE115821_Auslander/GSE115821_MGH_counts.csv", check.names = F)
AuslanderExpr <- Auslander0[!duplicated(Auslander0$Geneid), ]
rownames(AuslanderExpr) <- AuslanderExpr$Geneid
AuslanderExpr <- AuslanderExpr[,grep("^M|^[0-9]", colnames(AuslanderExpr))]
colnames(AuslanderExpr) <- gsub(".bam", "", colnames(AuslanderExpr))

Auslander_res <- Calculate_Miracle(as.matrix(AuslanderExpr), platform= "gene") 

AuslanderMeta <- read.csv("Training Data/GSE115821_Auslander/clinical_data.csv")
AuslanderMeta$Sample_title <- gsub(".bam", "", AuslanderMeta$Sample_title)
rownames(AuslanderMeta) <- AuslanderMeta$Sample_title

# ES = t(gsva(as.matrix(AuslanderExpr),SelPath_Symb,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(as.matrix(AuslanderExpr))){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_Symb,function(x) mwwGST(as.matrix(AuslanderExpr)[,i], x)$nes )))
}
Auslander_GSE115821 <- cbind(Auslander_res, AuslanderMeta[rownames(Auslander_res), ], resMWW)

#### PRJEB25780 Kim2018
Kim0 <- read.csv("Training Data/PRJEB25780_Kim/gene.counts_symbol.csv", check.names = F)
rownames(Kim0) <- Kim0$Entrez
Kim0$symbol <- NULL
Kim0$Entrez <- NULL
KimExpr <- Kim0

Kim_res <- Calculate_Miracle(as.matrix(KimExpr), platform= "entrez") 

KimMeta <- read.csv("Training Data/PRJEB25780_Kim/clinical_data _R_NR.csv")
rownames(KimMeta) <- gsub("-", ".",KimMeta$patient)

# ES = t(gsva(as.matrix(KimExpr),SelPath_entrex,method="ssgsea"))
resMWW <- c()
for(i in 1:ncol(as.matrix(KimExpr))){
  resMWW <- rbind(resMWW, unlist(lapply(SelPath_entrex,function(x) mwwGST(as.matrix(KimExpr)[,i], x)$nes )))
}
Kim_PRJEB25780 <- cbind(Kim_res, KimMeta[rownames(Kim_res), ], resMWW)



save(Cho_GSE126044, Riaz_GSE91061, Gide_PRJEB23709, vanAllen,
     Ribas_GSE78220, UlloaMontoya_GSE35640, Chen, Jang_GSE40419, Ascierto_GSE67501,
     Ascierto_GSE79691, Auslander_GSE115821, Kim_PRJEB25780, 
     file= "Master_Datasets.RData")
load("GM_AntiPD1/Master_Datasets.RData")


save(Cho_GSE126044, Riaz_GSE91061, Gide_PRJEB23709, vanAllen,
     UlloaMontoya_GSE35640, Chen, Kim_PRJEB25780, 
     file= "GM_AntiPD1/Master_Datasets_Selected.RData")


#### GSE121810 Cloughesy
### Only few genes here, so I had to exclude it
# Hwang0 <- read.delim("Training Data/GSE136961_Hwang_Lung/GSE136961_raw_count.tsv", check.names = F)
# HwangExpr <- Hwang0[!duplicated(Hwang0$Gene), ]
# rownames(HwangExpr) <- HwangExpr$Gene
# HwangExpr <- HwangExpr[,grep("^D[0-9]|^N[0-9]", colnames(HwangExpr))]
# 
# Hwang_res <- Calculate_Miracle(as.matrix(HwangExpr), platform= "gene") 
# 
# HwangMeta <- read.csv("Training Data/GSE136961_Hwang_Lung/clinical_data.csv")
# rownames(HwangMeta) <- HwangMeta$Sample_title
# 
# Cloughesy_GSE121810 <- cbind(Hwang_res, HwangMeta[rownames(Hwang_res), ])






