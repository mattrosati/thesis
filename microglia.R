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
datamatrix <- as.matrix(data.raw)
datamatrix <- apply(datamatrix, 2, as.numeric)
rownames(datamatrix) <- rownames(data.raw)

#build expression set so that I can filter
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

# normalize and get p-values, not needed because data is already normalized
# names <- rownames(data)
# data.b <- lumiB(data, method="bgAdjust.affy")
# data.l <- lumiT(data.b, method="log2")
# data.n <- lumiN(data.l, method="quantile")
# 
# design <- model.matrix(~pd$population)
# fit <- lmFit(data.n, design)
# ebayes <- eBayes(fit)

#filter to remove control probes
data.filter <- nsFilter(data, require.entrez=TRUE, require.GOBP=FALSE, 
                        require.GOCC=FALSE, require.GOMF=FALSE, 
                        require.CytoBand=FALSE, remove.dupEntrez=FALSE, var.filter=FALSE, 
                        feature.exclude="^AFFX")$eset

#create data frame with fold change, sort by logFC, and remove duplicates
gene_expr <- data.frame(data.filter@featureData@data, data.filter@assayData[["exprs"]])
genenames <- rownames(gene_expr)
gene_expr$external_gene_name <- as.character(getSYMBOL(genenames, "hgu133a.db"))
gene_expr$avec <- c(rowMeans(gene_expr[, c(2, 6, 7, 8)]))
gene_expr$avea <- c(rowMeans(gene_expr[, c(1, 3, 4, 5)]))
gene_expr$logFC <- c(log2(with(gene_expr, avea/avec)))

gene_expr.sort <- gene_expr[order(-abs(gene_expr$logFC)),]
gene_expr.sort <- gene_expr.sort[!duplicated(gene_expr.sort$external_gene_name),]

#annotate
annotationset <- getBM(attributes=c("external_gene_name", "entrezgene", "ensembl_gene_id", "description"), 
                       filters="external_gene_name", values=gene_expr.sort$external_gene_name, mart=ensembl)

#make final and save
final <- full_join(annotationset, gene_expr.sort[,9:12])
final.sort <- final[order(-abs(final$logFC)),]
final.sort <- final.sort[!duplicated(final.sort$ensembl_gene_id),]
write.csv(final.sort, file = "~/Desktop/MacMicking Lab/IFN_microarray/Processed Data/microglialf.csv")

