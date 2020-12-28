library(data.table)
library(GSVA)

setwd("/export/cse02/rmall/AntiPD1_Challenge/APD1DC2021/Model_Building/")

#Run GSVA for all the pathways on Riaz datasets to get each immune related pathways score

#Load Riaz data
load("Processed data/Riaz etal/002.Riaz_Data_normalized.gene.counts.Rdata")
riaz_clinical_data <- read.table("Processed data/Riaz etal/clinical_data - all cohort.csv",header=TRUE,sep=",")
riaz_clinical_data$Sample_title <- as.vector(as.character(riaz_clinical_data$Sample_title))
riaz_clinical_data$Response1 <- as.vector(as.character(riaz_clinical_data$Response1))
riaz_clinical_data$Response2 <- as.vector(as.character(riaz_clinical_data$Response2))
rownames(riaz_clinical_data) <- riaz_clinical_data$Sample_title

#Load all pathways 
load("Required Files/Selected.pathways.3.4.RData")

#Log2 the expression data
Expression.data = log(dataNorm +1, 2)
available_genes = rownames(Expression.data)
Gene.list = Selected.pathways
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Expression.data))]

##Perform ssGSEA for each pathway
ES = gsva(Expression.data,Gene.list,method="ssgsea")   # ES is the enrichment score

#Make plot of all the pathways and compare the GSVA score for responders vs non-responders
plot_df <- NULL
for (i in 1:nrow(ES))
{
  samples <- colnames(ES)
  pathway_i <- rep(rownames(ES)[i],length(samples))
  response_i <- riaz_clinical_data[samples,]$Response1
  pathway_scores <- as.vector(ES[i,])
  temp <- cbind(samples, pathway_i, pathway_scores, response_i)
  plot_df <- rbind(plot_df, temp)
}
plot_df <- as.data.frame(plot_df)
colnames(plot_df) <- c("Samples","Pathways","Scores","Response")
plot_df$Samples <- as.character(as.vector(plot_df$Samples))
plot_df$Pathways <- as.character(as.vector(plot_df$Pathways))
plot_df$Scores <- as.numeric(as.vector(plot_df$Scores))
plot_df$Response <- as.character(as.vector(plot_df$Response))

plot_df <- plot_df[plot_df$Response!=" UNK",]
g <- ggplot(data=plot_df, aes(x=Response, y=Scores, fill = Response)) + geom_boxplot() + facet_wrap( ~ Pathways, nrow=6, ncol=9)

#Compute the p-values using wilcox test between pathways scores for responders vs non-responders
graphLabels <- NULL
for (i in 1:nrow(ES))
{
  pathway <- rownames(ES)[i]
  pval <- wilcox.test(plot_df[plot_df$Pathways==pathway & plot_df$Response=="NR",]$Scores,plot_df[plot_df$Pathways==pathway & plot_df$Response=="R",]$Scores,exact=F)
  graphLabels <- rbind(graphLabels,cbind(pathway,signif(pval$p.value,4)))
}
graphLabels <- as.data.frame(graphLabels)
colnames(graphLabels) <- c("Pathways","Pval")
graphLabels$Pathways <- as.character(as.vector(graphLabels$Pathways))
graphLabels$Pval <- as.character(as.vector(graphLabels$Pval))
graphLabels$Pval <- paste0("p = ",graphLabels$Pval)

g1 <- g + geom_text(data = graphLabels, aes(x = 1.5, y = 0.5, label = Pval), inherit.aes = FALSE) +
        theme(strip.text.x = element_text(size = 8, colour = "darkgreen", angle = 0))
ggsave(filename="Results/Riaz_Pathway_Score_Comparison.pdf",plot=g1, device=pdf(),width=20, height=14, units="in")
dev.off()
