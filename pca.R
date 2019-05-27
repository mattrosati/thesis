library(dplyr)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(reshape2)

logfold <- read.csv(file = "~/Desktop/MacMicking Lab/IFN_microarray/full_logchange_p.csv",
                       header = T)
logfold$X <- NULL
genenames <- logfold[,1]
celltype <- colnames(logfold)[-1]
logfold.nogene <- logfold
rownames(logfold.nogene) <- logfold.nogene$gene
logfold.nogene$gene <- NULL

logfold.pca <- prcomp(logfold.nogene, scale. = T)
ggbiplot(logfold.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '')

genes.s <- c()
genes.a <- c()
for (i in c(2:9)) {
  temps <- logfold[logfold[,i] <= -1, c(1, i)]
  tempa <- logfold[logfold[,i] >= 1, c(1, i)]
  genes.s <- c(genes.s, as.character(temps[,1]))
  genes.a <- c(genes.a, as.character(tempa[,1]))
}

genes.a <- genes.a[!duplicated(genes.a)]
genes.s <- genes.s[!duplicated(genes.s)]

logfold.s <- logfold[match(genes.s, logfold$gene), ]
logfold.a <- logfold[match(genes.a, logfold$gene), ]

logfold.s.pca <- prcomp(logfold.s[2:9], scale. = F)
maximumchange.celltype.s <- factor(colnames(logfold.s[2:9])[apply(logfold.s[2:9],1,which.max)])
logfold.a.pca <- prcomp(logfold.a[2:9], scale. = T)
maximumchange.celltype.a <- factor(colnames(logfold.a[2:9])[apply(logfold.a[2:9],1,which.max)])

ggbiplot(logfold.s.pca, obs.scale = 1, var.scale = 1, groups = maximumchange.celltype.s, 
         ellipse = TRUE, circle = TRUE) + scale_color_discrete(name = '')
ggbiplot(logfold.a.pca, obs.scale = 1, var.scale = 1, groups = maximumchange.celltype.a, 
         ellipse = TRUE, circle = TRUE) + scale_color_discrete(name = '')






mat.s <- as.matrix(logfold.s[2:9])
mat.a <- as.matrix(logfold.a[2:9])
heatmap(mat.s)
heatmap(mat.a)


