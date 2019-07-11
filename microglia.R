source("packages.R")

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

#transform to log2
datamatrix <- log2(datamatrix)

#build expression set so that I can filter
experimentData <- new("MIAME", name = "Matteo Rosati",
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

# normalize, not needed because data is already normalized
# names <- rownames(data)
# data.b <- lumiB(data, method="bgAdjust.affy")
# data.l <- lumiT(data.b, method="log2")
# data.n <- lumiN(data.l, method="quantile")

# get p-values
pd$population <- relevel(pd$population, "control")    #so that in design matrix, active comes up as 1
design <- model.matrix(~pd$population)
fit <- lmFit(data, design)
ebayes <- eBayes(fit)

#filter to remove control probes
data.filter <- nsFilter(data, require.entrez=TRUE, require.GOBP=FALSE, 
                        require.GOCC=FALSE, require.GOMF=FALSE, 
                        require.CytoBand=FALSE, remove.dupEntrez=FALSE, var.filter=FALSE, 
                        feature.exclude="^AFFX")$eset

#toptable
tab <- topTable(ebayes, coef = 2, adjust.method = "BH", n = 100000)

#annotate gene expression data
gene_expr <- data.frame(data.filter@featureData@data, data.filter@assayData[["exprs"]])
gene_expr$PROBEID <- rownames(gene_expr)
annotate <- AnnotationDbi::select(hgu133a.db, keys = gene_expr$PROBEID, 
                   columns = c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID"), keytype = "PROBEID")
annotate <- annotate[!is.na(annotate$SYMBOL),]
gene_expr <- left_join(gene_expr, annotate)

#merge with tab, sort by logFC, remove duplicates
tab$PROBEID <- rownames(tab)
tab2 <- left_join(gene_expr, tab)
tab2.sort <- tab2[order(-abs(tab2$logFC), tab2$P.Value),]
tab2.sort <- tab2.sort[!duplicated(tab2.sort$SYMBOL),]

#make final and save
final <- tab2.sort[,c(10, 11, 12, 13, 6, 2, 8, 7, 1, 3, 4, 5, 14, 15, 16, 17, 18, 19)]
write.csv(final, file = "~/Desktop/MacMicking Lab/IFN_microarray/Processed Data/microglialf.csv")

