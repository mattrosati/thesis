library(lumi)
library(affyPLM)
library(simpleaffy)
library(dplyr)
library(biomaRt)
library(genefilter)
library(hgu133a.db)
library(limma)
library(annotate)

ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#import and manipulate
raw <- read.delim("~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/microglia/microglia.txt")
probes <- raw[-1,1]
datacols <- grep("24hrs.2", colnames(raw))
data.raw <- as.data.frame(raw[-1, c(datacols)])
rownames(data.raw) <- probes
colnames(data.raw) <- c("active1", "control2", "active2", "active3", "active4",
                        "control1", "control4", "control3")

#build expression set
datamatrix <- as.matrix(data.raw)
datamatrix <- apply(datamatrix, 2, as.numeric)
rownames(datamatrix) <- rownames(data.raw)
experimentData <- new("MIAME", name = "Pierre Fermat", 
                      lab = "Francis Galton Lab", 
                      contact = "pfermat@lab.not.exist", 
                      title = "Smoking-Cancer Experiment", 
                      abstract = "An example ExpressionSet", 
                      url = "www.lab.not.exist", 
                      other = list(notes = "Created from text files"))
pd <- data.frame(population = factor(c("active","control", "active", 
                                       "active","active", "control", 
                                       "control", "control")), 
                 replicate = c(1, 2, 2, 3, 4, 1, 4, 3))
rownames(pd) <- colnames(datamatrix)
vl <- data.frame(labelDescription = c("control is control, active is treatment", "replicate number"))
phenoData <- new("AnnotatedDataFrame", data = pd, varMetadata = vl)
data <- new("ExpressionSet", exprs = datamatrix, phenoData = phenoData, 
            experimentData = experimentData, annotation = "hgu133a")

names <- rownames(data)
data.b <- lumiB(data, method="bgAdjust.affy")
data.l <- lumiT(data.b, method="log2")
data.n <- lumiN(data.l, method="quantile")

#p-values
design <- model.matrix(~pd$population)
fit <- lmFit(data.n, design)
ebayes <- eBayes(fit)

tab <- tab[!is.na(tab$external_gene_name),]
annotated.sort <- annotated[order(-abs(annotated$logFC), annotated$P.Value),]
annotated.sort <- annotated.sort[!duplicated(annotated.sort$ensembl_gene_id),]
final <- annotated.sort
