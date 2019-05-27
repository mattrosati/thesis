library(dplyr)

data_list <- list.files("~/Desktop/MacMicking Lab/IFN_microarray/Processed Data",full.names = T)

#read files
a549 <- read.csv(file = data_list[1], header = T)
aaleb <- read.csv(file = data_list[2], header = T)
arpe <- read.csv(file = data_list[3], header = T)
fibroblast <- read.csv(file = data_list[4], header = T)
hepatocyte <- read.csv(file = data_list[5], header = T)
# intestine <- read.csv(file = data_list[8], header = T)
keratinocyte <- read.csv(file = data_list[9], header = T)
microglial <- read.csv(file = data_list[10], header = T)
thyroid <- read.csv(file = data_list[11], header = T)

#take out gene id and log fc
# a549p <- a549[, c("ensembl_gene_id", "logFC")]
# aalebp <- aaleb[, c("ensembl_gene_id", "logFC")]
# arpep <- arpe[, c("ensembl_gene_id", "logFC")]
# fibroblastp <- fibroblast[, c("ensembl_gene_id", "logFC")]
# hepatocytep <- hepatocyte[, c("ensembl_gene_id", "logFC")]
# # intestinep <- intestine[intestine$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
# keratinocytep <- keratinocyte[, c("ensembl_gene_id", "logFC")]
# microglialp <- microglial[, c("ensembl_gene_id", "logFC")]
# thyroidp <- thyroid[, c("ensembl_gene_id", "logFC")]
a549p <- a549[, c("ensembl_gene_id", "logFC")]
aalebp <- aaleb[aaleb$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
arpep <- arpe[arpe$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
fibroblastp <- fibroblast[fibroblast$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
hepatocytep <- hepatocyte[hepatocyte$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
# intestinep <- intestine[intestine$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
keratinocytep <- keratinocyte[keratinocyte$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
microglialp <- microglial[microglial$P.Value <= 0.1 , c("ensembl_gene_id", "logFC")]
thyroidp <- thyroid[thyroid$P.Value <= 0.1, c("ensembl_gene_id", "logFC")]

#eliminate missing genes
a549p <- a549p[!is.na(a549p$ensembl_gene_id),]
aalebp <- aalebp[!is.na(aalebp$ensembl_gene_id),]
arpep <- arpep[!is.na(arpep$ensembl_gene_id),]
fibroblasp <- fibroblastp[!is.na(fibroblastp$ensembl_gene_id),]
hepatocytep <- hepatocytep[!is.na(hepatocytep$ensembl_gene_id),]
# intestinep <- intestinep[!is.na(intestinep$ensembl_gene_id),]
keratinocytep <- keratinocytep[!is.na(keratinocytep$ensembl_gene_id),]
microglialp <- microglialp[!is.na(microglialp$ensembl_gene_id),]
thyroidp <- thyroidp[!is.na(thyroidp$ensembl_gene_id),]

#change title
colnames(a549p) <- c("gene", "a549")
colnames(aalebp) <- c("gene", "aaleb")
colnames(arpep) <- c("gene", "arpe")
colnames(fibroblastp) <- c("gene", "fibroblast")
colnames(hepatocytep) <- c("gene", "hepatocyte")
colnames(keratinocytep) <- c("gene", "keratinocyte")
colnames(microglialp) <- c("gene", "microglial")
colnames(thyroidp) <- c("gene", "thyroid")

#join all
temps1 <- full_join(a549p, full_join(aalebp, arpep))
temps2 <- full_join(fibroblastp, full_join(hepatocytep, keratinocytep))
temps3 <- full_join(microglialp, thyroidp)
globals <- full_join(temps1, full_join(temps2, temps3))

#fix duplicates and nas and negatvie numbers
globals <- globals[!is.na(globals$gene),]
globals.nodup <- globals[!duplicated(globals$gene),]
globals.nodup[is.na(globals.nodup)] <- 0
globals.nodup.neg <- cbind(globals.nodup[,c(1,2)], -1*globals.nodup[,c(3:9)])

write.csv(globals.nodup.neg, file = "~/Desktop/MacMicking Lab/IFN_microarray/full_logchange_p.csv")
