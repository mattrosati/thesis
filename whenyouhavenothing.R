library(dplyr)
library(hgug4112a.db)
library(annotate)
library(biomaRt)
library(genefilter)
library(limma)

ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
a <- listAttributes(ensembl)
f <- listFilters(ensembl)

active <- read.delim("~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/a549/active.txt")
control <- read.delim("~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/a549/control.txt")
colnames(control)[colnames(control)=="VALUE"] <- "control"
colnames(active)[colnames(active)=="VALUE"] <- "active"

qc <- quantile(control$control, probs = 0.75)
qa <- quantile(active$active, probs = 0.75)
active$active <- active$active/qa
control$control <- control$control/qc

controlsc <- grep("\\(", control$Reporter.Identifier)
controlsa <- grep("\\(", active$Reporter.Identifier)
datac <- control[-c(controlsc), ]
dataa <- active[-c(controlsa),]
probesc <- grep("A_", datac$Reporter.Identifier)
probesa <- grep("A_", dataa$Reporter.Identifier)
dataa <- dataa[c(probesa),]
datac <- datac[c(probesc),]
data <- full_join(datac, dataa)

data$active <- log2(data$active)
data$control <- log2(data$control)
data$logFC <- data$active - data$control
# data.logs <- data[, c("Reporter.Identifier", "logFC")]
data.logs.sort <- data[order(abs(data$logFC), decreasing = TRUE),]
data.logs.sort <- data.logs.sort[!duplicated(data.logs.sort$Reporter.Identifier),]

# rownames(data.logs.sort) <- data.logs.sort$Reporter.Identifier

annotationset <- getBM(attributes=c("agilent_wholegenome_4x44k_v1", "ensembl_gene_id"),
                       filters="agilent_wholegenome_4x44k_v1", values=data.logs.sort$Reporter.Identifier, mart=ensembl)

colnames(annotationset)[colnames(annotationset)=="agilent_wholegenome_4x44k_v1"] <- "Reporter.Identifier"
annotated <- full_join(annotationset, data.logs.sort)
annotated <- annotated[!is.na(annotated$ensembl_gene_id),]
annotated.sort <- annotated[order(abs(annotated$logFC), decreasing = TRUE),]
annotated.sort <- annotated.sort[!duplicated(annotated.sort$ensembl_gene_id),]

annotationset <- getBM(attributes=c("external_gene_name", "entrezgene", "ensembl_gene_id", "description"), 
                       filters="ensembl_gene_id", values=annotated.sort$ensembl_gene_id, 
                       mart=ensembl)

final <- full_join(annotationset, annotated.sort[,c(2:5)])
write.csv(final, file = "~/Desktop/MacMicking Lab/IFN_microarray/Processed Data/a549f2.csv")




