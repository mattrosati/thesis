library("biomaRt")
library(tidyverse)
library("org.Hs.eg.db")
library("topGO")
library(reshape2)
library(GO.db)
library(annotate)

options(stringsAsFactors = FALSE)
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
goterms <- Term(GOTERM)

#import BioPLex data
bioplex <- read.delim("~/Desktop/MacMicking Lab/IFN_microarray/BioPlex/BioPlex_interactionList_v4a.tsv")

#read/manipulate matrix for entrez gene querying
suppressed.m <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/suppressedmatrix.csv")
activated.m <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/activatedmatrix.csv")
entrez.s.m <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", 
                  values = suppressed.m$gene, mart = ensembl)
entrez.a.m <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", 
                    values = activated.m$gene, mart = ensembl)
colnames(entrez.s.m)[1] <- "gene"
colnames(entrez.a.m)[1] <- "gene"
suppressed.entrez <- full_join(suppressed.m, entrez.s.m)[2:11]
suppressed.entrez <- suppressed.entrez[!duplicated(suppressed.entrez$entrezgene),]
activated.entrez <- full_join(activated.m, entrez.a.m)[2:11]
activated.entrez <- activated.entrez[!duplicated(activated.entrez$entrezgene),]

write.csv(activated.entrez, 
          file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/a_entrez.csv")
write.csv(suppressed.entrez, 
          file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/s_entrez.csv")

# loop through data to get interactions in Python, using R would take ages, 
# pipeline file is interaction_loop.py


# tf.a <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/tfactivated.csv")
# 
# tf.a.unique <- tf.a[rowSums(tf.a[3:10]) == 1, -1]
# entrez.a <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", 
#                      values = tf.a.unique$hgnc_symbol, mart = ensembl)
# tf.a.unique.entrez <- full_join(entrez.a, tf.a.unique)

list.interactions <- sapply(paths, function(x) fromJSON(rawToChar(GET(url = url, path = x)$content)))
length <- lapply(list.interactions, length)
interaction.v <- list.interactions[!length==0]
interaction.v <- interaction.v[c(1:2284)]
interactions <- do.call(what = rbind, args = interaction.v)
interaction.names <- interactions[, c(1,2)]
rownames(interaction.names) <- NULL




uncharacterized <- suppressed.unique$entrezgene[length == 0]
cell <- c()
for (i in c(1:length(uncharacterized))) {
  entrez <- uncharacterized[i]
  celltype <- colnames(suppressed.unique)[which(suppressed.unique[which(suppressed.unique$entrezgene == entrez),] == 1)]
  cell <- c(cell, celltype)
}
uncharacterized.df <- data.frame(entrez = uncharacterized, 
                                 symbol = getSYMBOL(as.character(uncharacterized), data = "org.Hs.eg.db"), 
                                 cell = cell)
uncharacterized.df <- uncharacterized.df[,c(2:3)]

write.csv(uncharacterized.df, file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/s_unique_uncharacterized.csv")

one <- c()
two <- c()
three <- c()
for (i in c(1:nrow(interaction.names))) {
  main <- interaction.names[i, which(is.element(interaction.names[i,], suppressed.unique$entrezgene))[1]]
  interact <- interaction.names[i, which(interaction.names[i,] != main)]
  if (interaction.names[i, 1] == interaction.names[i, 2]) {
    interact <- main
  }
  celltype.main <- colnames(suppressed.entrez)[which(suppressed.entrez[which(suppressed.entrez$entrezgene == main),] == 1)]
  state.interact <- suppressed.entrez[which(suppressed.entrez$entrezgene == interact), celltype.main]
  if (length(state.interact) > 0) {
    one <- c(one, main)
    two <- c(two, interact)
    three <- c(three, celltype.main)
  }
}

suppressed.unique.interactions <- data.frame(gene1 = one, gene2 = two, cell = three)
suppressed.unique.interactions$symbol1 <- getSYMBOL(suppressed.unique.interactions$gene1, data = "org.Hs.eg.db")
suppressed.unique.interactions$symbol2 <- getSYMBOL(suppressed.unique.interactions$gene2, data = "org.Hs.eg.db")

write.csv(suppressed.unique.interactions[,c(3:5)], 
          file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/s.u.interactions.csv")

