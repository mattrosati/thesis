library("biomaRt")
library("dplyr")
library("org.Hs.eg.db")
library("topGO")
library(reshape2)
library(ggplot2)
library(httr)
library(jsonlite)
library(GO.db)

options(stringsAsFactors = FALSE)
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
goterms <- Term(GOTERM)

suppressed.m <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/suppressedmatrix.csv")
activated.m <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/activatedmatrix.csv")
entrez.a.m <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", 
                  values = activated.m$gene, mart = ensembl)
colnames(entrez.a.m)[1] <- "gene"
activated.entrez <- full_join(activated.m, entrez.a.m)[2:11]
activated.entrez <- activated.entrez[!duplicated(activated.entrez$entrezgene),]


tf.a <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/tfactivated.csv")

tf.a.unique <- tf.a[rowSums(tf.a[3:10]) == 1, -1]
entrez.a <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", 
                     values = tf.a.unique$hgnc_symbol, mart = ensembl)
tf.a.unique.entrez <- full_join(entrez.a, tf.a.unique)

#fetch interaction data
url <- "http://bioplex.hms.harvard.edu"
paths <- c(paste("bioplexDisplay/api/api.php?pintLow=0.9&geneQuery=", tf.a.unique.entrez$entrezgene, sep = ""))
list.interactions <- sapply(paths, function(x) fromJSON(rawToChar(GET(url = url, path = x)$content)))
interactions <- do.call(what = rbind, args = list.interactions)
interaction.names <- interactions[, c(1,2)]
rownames(interaction.names) <- NULL

#uncharacterized
length <- lapply(list.interactions, length)
uncharacterized <- tf.a.unique.entrez$hgnc_symbol[length == 0]
cell <- c()
for (i in c(1:length(uncharacterized))) {
  genename <- uncharacterized[i]
  celltype <- colnames(tf.a.unique.entrez)[which(tf.a.unique.entrez[which(tf.a.unique.entrez$hgnc_symbol == genename),] == 1)]
  cell <- c(cell, celltype)
}
uncharacterized.df <- data.frame(gene = uncharacterized, cell = cell)

write.csv(uncharacterized.df, file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/tfa_unique_uncharacterized.csv")

#fetch same state interactions df
one <- c()
two <- c()
three <- c()
for (i in c(1:nrow(interaction.names))) {
  main <- interaction.names[i, which(is.element(interaction.names[i,], tf.a.unique.entrez$entrezgene))[1]]
  interact <- interaction.names[i, which(interaction.names[i,] != main)]
  if (interaction.names[i, 1] == interaction.names[i, 2]) {
    interact <- main
  }
  celltype.main <- colnames(activated.entrez)[which(activated.entrez[which(activated.entrez$entrezgene == main),] == 1)]
  state.interact <- activated.entrez[which(activated.entrez$entrezgene == interact), celltype.main]
  if (length(state.interact) > 0) {
    one <- c(one, main)
    two <- c(two, interact)
    three <- c(three, celltype.main)
  }
}

activated.tf.interactions <- data.frame(gene1 = one, gene2 = two, cell = three)
activated.tf.interactions$symbol1 <- getSYMBOL(activated.tf.interactions$gene1, data = "org.Hs.eg.db")
activated.tf.interactions$symbol2 <- getSYMBOL(activated.tf.interactions$gene2, data = "org.Hs.eg.db")

write.csv(activated.tf.interactions[,c(3:5)], 
          file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/a.tf.interactions.csv")


