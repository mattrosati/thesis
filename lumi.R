library(lumi)
library(biomaRt)
library(lumiHumanIDMapping)
library(affy)
library(illuminaHumanv4.db)
library(annotate)

ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

data_list <- list.files("~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/keratinocytes",full.names = T)
data <- lumiR.batch(data_list, lib.mapping = 'lumiHumanIDMapping')

data <- data[, c("5688887011_A", "5688887011_G", "5688887011_I", "5688879024_A", "5688879024_E", "5688879024_G")]

data.b <- lumiB(data, method="forcePositive")
data.t <- lumiT(data.b, method="log2")
data.n <- lumiN(data.t, method="quantile")

data.filter <- data.n

genenames <- tab$PROBE_ID
tab$external_gene_name <- as.character(getSYMBOL(genenames, "illuminaHumanv4.db"))
tab <- tab[!is.na(tab$external_gene_name),]



