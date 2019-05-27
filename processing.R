library(affy)
library(affyPLM)
library(simpleaffy)
library(dplyr)
library(biomaRt)
library(genefilter)
library(hgu133a.db)
library(limma)
library(annotate)

ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#import
data_list <- list.files("~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/thyroid",full.names = T)
data <- ReadAffy(filenames = data_list) 

#make phenodata
sampleNames(data.filter) <- c("control1", "control2", "control3", "active1", "active2", "active3")
pd <- data.frame(population = factor(c("control", "control", "control", "active", "active", "active")), replicate = c(1, 2, 3, 1, 2, 3))
rownames(pd) <- sampleNames(data)
vl <- data.frame(labelDescription = c("control is control, active is treatment", "replicate number"))
phenoData(data.filter) <- new("AnnotatedDataFrame", data = pd, varMetadata = vl)


#quality check
dataPLM <- fitPLM(data)
Mbox(dataPLM, main="RLE", ylim = c(-0.3, 0.3), outline=F, col="mistyrose", las=3,
     whiskylty=0, staplelty=0)
boxplot(dataPLM, main="NUSE", ylim = c(0.95, 1.22), outline=F, col="lightblue", las=3, 
        whisklty=0, staplelty=0)

#rma
data.rma <- affy::rma(data)
indexc <- which(data.rma$population == "control")
indexa <- which(data.rma$population == "active")

#filtering
data.filter <- nsFilter(data.rma, require.entrez=TRUE, require.GOBP=FALSE, 
                             require.GOCC=FALSE, require.GOMF=FALSE, 
                             require.CytoBand=FALSE, remove.dupEntrez=F, var.cutoff=0.5, 
                             feature.exclude="^AFFX")$eset

#p-values/log2change
# data6 <- as.data.frame(Biobase::exprs(data.filter[, c(1,2)]))
# data7<- as.data.frame(Biobase::exprs(data.filter[, c(3,4)]))
# data6$log2change <- data6$huh6a - data6$huh6c
# data$log2change <- data7$huh7a - data7$huh7c
design <- model.matrix(~data.filter$population)
fit <- lmFit(data.filter, design)
ebayes <- eBayes(fit)

#toptable
tab <- topTable(ebayes, coef = 2, adjust.method = "BH", n = 100000)

#annotate both tab and data.filter
gene_expr <- data.frame(data.filter@featureData@data, data.filter@assayData[["exprs"]])
genenames <- rownames(tab)
tab$external_gene_name <- as.character(getSYMBOL(genenames, "hgu133plus2.db"))
tab2 <- merge(gene_expr, tab, by=0, all.x = F, all.y = T)
tab2.sort <- tab2[order(-abs(tab2$logFC), tab2$P.Value),]
tab2.sort <- tab2.sort[!duplicated(tab2.sort$external_gene_name),]
annotationset <- getBM(attributes=c("external_gene_name", "entrezgene", "ensembl_gene_id", "description"), filters="external_gene_name", values=tab2.sort$external_gene_name, mart=ensembl)
colnames(tab2.sort)[c(2, 3, 4, 5, 6, 7)] <- c("control1", "active1", "control2", "active2", "control3", "active3")

#make final and save
final <- full_join(annotationset, tab2.sort[,2:14])
final.sort <- final[order(-abs(final$logFC), final$P.Value),]
final.sort <- final.sort[!duplicated(final.sort$ensembl_gene_id),]
write.csv(final.sort, file = "~/Desktop/MacMicking Lab/IFN_microarray/Processed Data/thyroidf2.csv")




# final.control <- data.frame(data.logs.control$means, row.names = row.names(data.logs.control))
# final.active <- data.frame(data.logs.active$means, row.names = row.names(data.logs.active))
# colnames(final.control) <- "log2rma.mean"
# colnames(final.active) <- "log2rma.mean"
# final.control$probe_id <- rownames(final.control)
# final.active$probe_id <- rownames(final.active)
# 
# annotate.control <- getBM(attributes=c("affy_hg_u133_plus_2", "ensembl_gene_id", "description"), filters="affy_hg_u133_plus_2", values=final.control$probe_id, mart=ensembl)
# annotate.active <- getBM(attributes=c("affy_hg_u133_plus_2", "ensembl_gene_id", "description"), filters="affy_hg_u133_plus_2", values=final.active$probe_id, mart=ensembl)
# 
# colnames(annotate.active)[1] <- "probe_id"
# colnames(annotate.control)[1] <- "probe_id"
# final.control <- inner_join(annotate.control, final.control)
# final.active <- inner_join(annotate.active, final.active)

# data_list <- list.files("~/Desktop/MacMicking Lab/IFN_microarray/Raw Data/aaleb",full.names = T)
# data <- read.celfiles(data_list)
# 
# ph = data@phenoData
# ph@data[ ,1] = c("control1","control2", "active1","active2")
# 
# for (i in 1:4)
# {
#   image(data[,i],main=ph@data$sample[i])
# }
# 
# boxplot(data,which='pm',col='red',names=ph@data$index)