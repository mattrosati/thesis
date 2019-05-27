library("biomaRt")
library("dplyr")
library("org.Hs.eg.db")
library("topGO")
library(reshape2)
library(ggplot2)
library('plyr')
library("openxlsx")

ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

activated <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/activatedmatrix.csv")
suppressed <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/suppressedmatrix.csv")
terms <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/goterms.txt", sep = "")

cell_types <- c(colnames(activated)[3:10])
names <- activated$gene

result_s <- vector("list", 8)
result_a <- vector("list", 8)

for (c in cell_types) {
  v_a = factor(activated[,c])
  v_s = factor(suppressed[,c])
  names(v_a) <- names
  names(v_s) <- names
  AGOdata <- new( "topGOdata", ontology = "BP", allGenes = v_a, nodeSize = 10, annot = annFUN.org, 
                  mapping = "org.Hs.eg.db", ID="ensembl" )
  SGOdata <- new( "topGOdata", ontology = "BP", allGenes = v_s, nodeSize = 10, annot = annFUN.org, 
                  mapping = "org.Hs.eg.db", ID="ensembl" )
  goTestA <- runTest( AGOdata, algorithm = "elim", statistic = "fisher" )
  goTestS <- runTest( SGOdata, algorithm = "elim", statistic = "fisher" )
  result_a[[match(c, cell_types)]] <- GenTable( AGOdata, goTestA)
  result_s[[match(c, cell_types)]] <- GenTable( SGOdata, goTestS)
}

go_list_s <- c()
go_list_a <- c()

#making GO descriptions dataframe
goterms_a <- data.frame(a549=result_a[[1]]$GO.ID, aaleb=result_a[[2]]$GO.ID, arpe=result_a[[3]]$GO.ID, fibroblast=result_a[[4]]$GO.ID)
goterms_a <- data.frame(goterms_a, hepatocyte=result_a[[5]]$GO.ID, keratinocyte=result_a[[6]]$GO.ID, microglial=result_a[[7]]$GO.ID, thyroid=result_a[[8]]$GO.ID)
goterms_s <- data.frame(a549=result_s[[1]]$GO.ID, aaleb=result_s[[2]]$GO.ID, arpe=result_s[[3]]$GO.ID, fibroblast=result_s[[4]]$GO.ID)
goterms_s <- data.frame(goterms_s, hepatocyte=result_s[[5]]$GO.ID, keratinocyte=result_s[[6]]$GO.ID, microglial=result_s[[7]]$GO.ID, thyroid=result_s[[8]]$GO.ID)
map = setNames(terms$x, rownames(terms))
goterms_a[] <- map[unlist(goterms_a)]
goterms_s[] <- map[unlist(goterms_s)]

wb <- createWorkbook("godesc")
addWorksheet(wb, "activated")
addWorksheet(wb, "suppressed")
writeData(wb, "activated", goterms_a, startCol = 2, startRow = 3)
writeData(wb, "suppressed", goterms_s, startCol = 2, startRow = 3)

saveWorkbook(wb, "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/godesc.xlsx", overwrite = TRUE)

#continue with rest of GO analysis
for (i in c(1:8)) {
  go_list_s <- c(go_list_s, result_s[[i]]$GO.ID)
  go_list_a <- c(go_list_a, result_a[[i]]$GO.ID)
}
go_list_s <- unique(go_list_s)
go_list_a <- unique(go_list_a)
df_s <- data.frame(go_list_s)
df_a <- data.frame(go_list_a)

for (i in c(1:8)) {
  factor_s <- factor(go_list_s %in% result_s[[i]]$GO.ID, labels = c(0, 1))
  factor_a <- factor(go_list_a %in% result_a[[i]]$GO.ID, labels = c(0, 1))
  num_s <- as.numeric(levels(factor_s)[factor_s])
  num_a <- as.numeric(levels(factor_a)[factor_a])
  df_s <- cbind(df_s, num_s)
  df_a <- cbind(df_a, num_a)
  colnames(df_s)[colnames(df_s) == "num_s"] <- cell_types[i]
  colnames(df_a)[colnames(df_a) == "num_a"] <- cell_types[i]
}
colnames(df_s)[colnames(df_s) == "go_list_s"] <- "go_term"
colnames(df_a)[colnames(df_a) == "go_list_a"] <- "go_term"

df_s <- df_s[df_s$go_term != "GO:0010469",]
df_a <- df_a[df_a$go_term != "GO:0010469",]

df_s$go <- goterms[df_s$go_term]
df_a$go <- goterms[df_a$go_term]

df_s <- df_s[order(df_s$a549, df_s$aaleb, df_s$arpe, df_s$fibroblast, df_s$hepatocyte, df_s$keratinocyte,
                   df_s$microglial, df_s$thyroid),]
df_a <- df_a[order(df_a$a549, df_a$aaleb, df_a$arpe, df_a$fibroblast, df_a$hepatocyte, df_a$keratinocyte,
                   df_a$microglial, df_a$thyroid),]

df_s_melt <- melt(df_s)
df_a_melt <- melt(df_a)

ggplot(df_s_melt, aes(variable, go_term)) + 
  geom_tile(aes(fill = factor(value)), colour = "white") + 
  scale_fill_manual(values=c("lightblue1", "dodgerblue3")) + 
  theme(axis.text.x  = element_text(angle=90, hjust = 1, size=10), legend.position = "none") +
  labs(title = "Suppressed Genes", x = "Cell Type", y="GO Term")
ggsave(filename = "suppressed_heatmapv2.jpg", width = 5, height = 10.92,
path="~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/")
ggplot(df_a_melt, aes(variable, go_term)) + 
  geom_tile(aes(fill = factor(value)), colour = "white") + 
  scale_fill_manual(values=c("lightblue1", "dodgerblue3")) + 
  theme(axis.text.x  = element_text(angle=90, hjust = 1, size=10), legend.position = "none") +
  labs(title = "Activated Genes", x = "Cell Type", y="GO Term")
ggsave(filename = "activated_heatmapv2.jpg", width = 5, 
       path="~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/")

#anova has to be done on fractional table! correct all of this!!!
s.sums <- data.frame(df_s$go_term, sums=rowSums(df_s[,2:9]))
a.sums <- data.frame(df_a$go_term,sums=rowSums(df_a[,2:9]))

s.sums$sums <- s.sums$sums/sum(s.sums$sums)
a.sums$sums <- a.sums$sums/sum(a.sums$sums)


anova.data <- data.frame(sum=c(s.sums$sums, a.sums$sums), 
                         state = factor(rep(c("Suppressed", "Activated"), 
                                            times = c(length(s.sums$sums), length(a.sums$sums)))))
fm1 <- aov(sum~state, data=anova.data)
anova(fm1)

#make bar graph
rs_a <- data.frame(Sum = rowSums(df_a[,c(2:9)]))
rs_s <- data.frame(Sum = rowSums(df_s[,c(2:9)]))

a_t <- as.data.frame(table(rs_a))
a_t$Freq <- a_t$Freq/sum(c(a_t$Freq))*100
s_t <- as.data.frame(table(rs_s))
s_t$Freq <- s_t$Freq/sum(c(s_t$Freq))*100
colnames(s_t) <- c('cell_number', 'suppressed')
colnames(a_t) <- c('cell_number', 'activated')

t <- left_join(a_t, s_t, by = 'cell_number')
t[is.na(t)] <- 0
melted <- melt(t)

barplot <- ggplot(melted, aes(fill=variable, y=value, x=cell_number)) + geom_bar(position="dodge", stat="identity")
barplot + ggtitle('Fraction of Significant GO Terms by Cell Type Occurrence') + 
  xlab('Number of occurrences across cell types') + ylab('Percent of GO Terms') + 
  scale_fill_discrete(name = 'State', labels = c('Activated', 'Suppressed'))
ggsave(filename = "go_bar.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")

#transcription factor analysis
suppressed.v2 <- suppressed
activated.v2 <- activated
rownames(suppressed.v2) <- suppressed$gene
rownames(activated.v2) <- activated$gene
activated.v2[,c(1:2)] <- NULL
suppressed.v2[,c(1:2)] <- NULL

suppressed.goterms <- select(org.Hs.eg.db, keys = rownames(suppressed.v2), columns = c("ENSEMBL", "GO"), keytype = "ENSEMBL")
activated.goterms <- select(org.Hs.eg.db, keys = rownames(activated.v2), columns = c("ENSEMBL", "GO"), keytype = "ENSEMBL")
suppressed.goterms <- suppressed.goterms[!is.na(suppressed.goterms$GO),]
activated.goterms <- activated.goterms[!is.na(activated.goterms$GO),]

suppressed.tf <- suppressed.goterms[suppressed.goterms$GO == "GO:0003700",]
activated.tf <- activated.goterms[activated.goterms$GO == "GO:0003700",]

suppressed.tf.db <- suppressed[suppressed$gene %in% suppressed.tf$ENSEMBL,]
activated.tf.db <- activated[activated$gene %in% activated.tf$ENSEMBL,]
suppressed.tf.db <- suppressed.tf.db[rowSums(suppressed.tf.db[3:10]) != 0, -1]
activated.tf.db <- activated.tf.db[rowSums(activated.tf.db[3:10]) != 0, -1]

rowsums.s <- data.frame(gene = suppressed.tf.db$gene, sums = rowSums(suppressed.tf.db[,2:9,]))
rowsums.a <- data.frame(gene = activated.tf.db$gene, suma = rowSums(activated.tf.db[,2:9]))

#anova, need to use partial sums
rowsums.s$sums <- rowsums.s$sums/sum(rowsums.s$sums)
rowsums.a$suma <- rowsums.a$suma/sum(rowsums.a$suma)

anova.data <- data.frame(sum=c(rowsums.s$sums, rowsums.a$suma), 
                         state = factor(rep(c("Suppressed", "Activated"), 
                                            times = c(length(rowsums.s$sums), length(rowsums.a$suma)))))
fm1 <- aov(sum~state, data=anova.data)
anova(fm1)


#make bar graph
a_tf <- as.data.frame(table(rowsums.a[,2]))
a_tf$Freq <- a_tf$Freq/sum(c(a_tf$Freq))*100
s_tf <- as.data.frame(table(rowsums.s[,2]))
s_tf$Freq <- s_tf$Freq/sum(c(s_tf$Freq))*100
colnames(s_tf) <- c('cell_number', 'suppressed')
colnames(a_tf) <- c('cell_number', 'activated')

tf <- left_join(a_tf, s_tf, by = 'cell_number')
tf[is.na(tf)] <- 0
melted_tf <- melt(tf)

barplot <- ggplot(melted_tf, aes(fill=variable, y=value, x=cell_number)) + geom_bar(position="dodge", stat="identity")
barplot + ggtitle('Fraction of Transcription Factors by Cell Type Occurrence') + 
  xlab('Number of occurrences across cell types') + ylab('Percent of Transcription Factors') + 
  scale_fill_discrete(name = 'State', labels = c('Activated', 'Suppressed'))
ggsave(filename = "tf_bar.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")


hgnc.s <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                values = suppressed.tf.db$gene, mart = ensembl)
hgnc.a <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                values = activated.tf.db$gene, mart = ensembl)
colnames(hgnc.s)[1] <- "gene"
colnames(hgnc.a)[1] <- "gene"

tf.s <- full_join(hgnc.s, suppressed.tf.db)[,2:10]
tf.a <- full_join(hgnc.a, activated.tf.db)[,2:10]
tf.s <- tf.s[!duplicated(tf.s$hgnc_symbol),]
tf.a <- tf.a[!duplicated(tf.a$hgnc_symbol),]

write.csv(tf.s, file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/tfsuppressed.csv")
write.csv(tf.a, file = "~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/tfactivated.csv")

tf.s.melt <- melt(tf.s)
tf.a.melt <- melt(tf.a)

ggplot(tf.s.melt, aes(variable, hgnc_symbol)) + 
  geom_tile(aes(fill = factor(value)), colour = "white") + 
  scale_fill_manual(values=c("lightblue1", "dodgerblue3")) + 
  theme(axis.text.x  = element_text(angle=90, hjust = 1, size=10), legend.position = "none") +
  labs(title = "Suppressed Transcription Factors", x = "Gene", y="Cell Type")
ggsave(filename = "suppressed_heatmapvtf.jpg", width = 5, height = 10.92,
       path="~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/")
ggplot(tf.a.melt, aes(variable, hgnc_symbol)) + 
  geom_tile(aes(fill = factor(value)), colour = "white") + 
  scale_fill_manual(values=c("lightblue1", "dodgerblue3")) + 
  theme(axis.text.x  = element_text(angle=90, hjust = 1, size=10), 
        axis.text.y  = element_text(size=8), legend.position = "none") +
  labs(title = "Activated Transcription Factors", x = "Gene", y="Cell Type")
ggsave(filename = "activated_heatmapvtf.jpg", width = 5,
       path="~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/")

#shared go term: GO:0010469
suppressed.goterms <- select(org.Hs.eg.db, keys = rownames(suppressed.v2), columns = c("ENSEMBL", "GO"), keytype = "ENSEMBL")
activated.goterms <- select(org.Hs.eg.db, keys = rownames(activated.v2), columns = c("ENSEMBL", "GO"), keytype = "ENSEMBL")
suppressed.goterms <- suppressed.goterms[!is.na(suppressed.goterms$GO),]
activated.goterms <- activated.goterms[!is.na(activated.goterms$GO),]

suppressed.term <- suppressed.goterms[suppressed.goterms$GO == "GO:0010469",]
activated.term <- activated.goterms[activated.goterms$GO == "GO:0010469",]

suppressed.term.db <- suppressed[suppressed$gene %in% suppressed.term$ENSEMBL,]
activated.term.db <- activated[activated$gene %in% activated.term$ENSEMBL,]
suppressed.term.db <- suppressed.term.db[rowSums(suppressed.term.db[3:10]) != 0, -1]
activated.term.db <- activated.term.db[rowSums(activated.term.db[3:10]) != 0, -1]

hgnc.term.s <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                values = suppressed.term.db$gene, mart = ensembl)
hgnc.term.a <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                values = activated.term.db$gene, mart = ensembl)
colnames(hgnc.term.s)[1] <- "gene"
colnames(hgnc.term.s)[1] <- "gene"

term.s <- full_join(hgnc.term.s, suppressed.term.db)[,2:10]
term.a <- full_join(hgnc.term.s, activated.term.db)[,2:10]
term.s <- term.s[!duplicated(term.s$hgnc_symbol),]
term.a <- term.a[!duplicated(term.a$hgnc_symbol),]



