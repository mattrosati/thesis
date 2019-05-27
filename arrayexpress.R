library(ArrayExpress)
library(dplyr)
library(hugene10sttranscriptcluster.db)
library(annotate)
library(biomaRt)
library(genefilter)
library(limma)


data_exp <- getAE("E-GEOD-36331", type = "processed", 
              path= "~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/arpe")
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#import and process
data_list <- list.files("~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/fibroblasts",full.names = T)
control1 <- as.data.frame(read.table(file=data_list[1], header = T, sep = "\t"))
control2 <- as.data.frame(read.table(file=data_list[2], header = T, sep = "\t"))
control3 <- as.data.frame(read.table(file=data_list[3], header = T, sep = "\t"))
active1 <- as.data.frame(read.table(file=data_list[4], header = T, sep = "\t"))
active2 <- as.data.frame(read.table(file=data_list[5], header = T, sep = "\t"))
active3 <- as.data.frame(read.table(file=data_list[6], header = T, sep = "\t"))
colnames(control1)[colnames(control1)=="VALUE"] <- "control1"
colnames(control2)[colnames(control2)=="VALUE"] <- "control2"
colnames(control3)[colnames(control3)=="VALUE"] <- "control3"
colnames(active1)[colnames(active1)=="VALUE"] <- "active1"
colnames(active2)[colnames(active2)=="VALUE"] <- "active2"
colnames(active3)[colnames(active3)=="VALUE"] <- "active3"
control1 <- control1[,c(1:2)]
control2 <- control2[,c(1:2)]
control3 <- control3[,c(1:2)]
active1 <- active1[,c(1:2)]
active2 <- active2[,c(1:2)]
active3 <- active3[,c(1:2)]
datac <- full_join(control1, full_join(control2, control3))
dataa <- full_join(active1, full_join(active2, active3))
df <- full_join(datac, dataa)
rownames(df) <- df$ID_REF
datamatrix <- as.matrix(df[, 2:6])
experimentData <- new("MIAME", name = "Pierre Fermat", 
                      lab = "Francis Galton Lab", 
                      contact = "pfermat@lab.not.exist", 
                      title = "Smoking-Cancer Experiment", 
                      abstract = "An example ExpressionSet", 
                      url = "www.lab.not.exist", 
                      other = list(notes = "Created from text files"))

#build pData
pd <- data.frame(population = factor(c("control","control", "control", "active","active")), replicate = c(1, 2, 3, 1, 2))
rownames(pd) <- colnames(datamatrix)
vl <- data.frame(labelDescription = c("control is control, active is treatment", "replicate number"))
phenoData <- new("AnnotatedDataFrame", data = pd, varMetadata = vl)
data <- new("ExpressionSet", exprs = datamatrix, phenoData = phenoData, 
                  experimentData = experimentData, annotation = "pd.hugene.1.0.st.v1")

#p-values
design <- model.matrix(~pd$population)
fit <- lmFit(data, design)
ebayes <- eBayes(fit)

#top
tab <- topTable(ebayes, coef = 2, adjust.method = "BH", n = 100000)

#renames
genenames <- rownames(tab)
tab$external_gene_name <- as.character(getSYMBOL(genenames, "hugene10sttranscriptcluster.db"))
annotationset <- getBM(attributes=c("external_gene_name", "entrezgene", "ensembl_gene_id", "description"), filters="external_gene_name", values=tab$external_gene_name, mart=ensembl)

#finish
final <- full_join(annotationset, tab)
final <- final[, c(1:10)]
write.csv(final, file = "~/Desktop/MacMicking Lab/IFN_microarray/Processed Data/arpef.csv")


