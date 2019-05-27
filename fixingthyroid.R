library(biomaRt)
library(dplyr)

ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

activated <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/activatedmatrix.csv")
suppressed <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/suppressedmatrix.csv")
cross <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/thyroid_genestocheck.csv", header = FALSE)

cross_ensembl <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters="hgnc_symbol", values=cross$V1, mart=ensembl)

colnames(cross) <- c("hgnc_symbol")
finalcross <- inner_join(cross_ensembl, cross)

finalcross <- unique(finalcross)

thyroid_a <- as.character(activated$gene[activated$thyroid == 1])
thyroid_s <- as.character(suppressed$gene[suppressed$thyroid == 1])

missing <- ifelse(c(finalcross$ensembl_gene_id) %in% thyroid_a | c(finalcross$ensembl_gene_id) %in% thyroid_s, 'yes', c(finalcross$ensembl_gene_id))
