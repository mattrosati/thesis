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

#convert to 0 or 1 and separate into a and s
a549fa <- transmute(a549p, gene = ensembl_gene_id, 
                    a549 = ifelse(a549p$logFC >= 1, 1, 0))
a549fs <- transmute(a549p, gene = ensembl_gene_id, 
                    a549 = ifelse(a549p$logFC <= -1, 1, 0))

aalebfs <- transmute(aalebp, gene = ensembl_gene_id, 
                    aaleb = ifelse(aalebp$logFC >= 1, 1, 0))
aalebfa <- transmute(aalebp, gene = ensembl_gene_id, 
                    aaleb = ifelse(aalebp$logFC <= -1, 1, 0))

arpefs <- transmute(arpep, gene = ensembl_gene_id, 
                     arpe = ifelse(arpep$logFC >= 1, 1, 0))
arpefa <- transmute(arpep, gene = ensembl_gene_id, 
                     arpe = ifelse(arpep$logFC <= -1, 1, 0))

fibroblastfs <- transmute(fibroblastp, gene = ensembl_gene_id, 
                    fibroblast = ifelse(fibroblastp$logFC >= 1, 1, 0))
fibroblastfa <- transmute(fibroblastp, gene = ensembl_gene_id, 
                    fibroblast = ifelse(fibroblastp$logFC <= -1, 1, 0))

hepatocytefs <- transmute(hepatocytep, gene = ensembl_gene_id, 
                    hepatocyte = ifelse(hepatocytep$logFC >= 1, 1, 0))
hepatocytefa <- transmute(hepatocytep, gene = ensembl_gene_id, 
                    hepatocyte = ifelse(hepatocytep$logFC <= -1, 1, 0))

keratinocytefs <- transmute(keratinocytep, gene = ensembl_gene_id, 
                    keratinocyte = ifelse(keratinocytep$logFC >= 1, 1, 0))
keratinocytefa <- transmute(keratinocytep, gene = ensembl_gene_id, 
                    keratinocyte = ifelse(keratinocytep$logFC <= -1, 1, 0))

microglialfs <- transmute(microglialp, gene = ensembl_gene_id, 
                    microglial = ifelse(microglialp$logFC >= 1, 1, 0))
microglialfa <- transmute(microglialp, gene = ensembl_gene_id, 
                    microglial = ifelse(microglialp$logFC <= -1, 1, 0))

thyroidfs <- transmute(thyroidp, gene = ensembl_gene_id, 
                    thyroid = ifelse(thyroidp$logFC >= 1, 1, 0))
thyroidfa <- transmute(thyroidp, gene = ensembl_gene_id, 
                    thyroid = ifelse(thyroidp$logFC <= -1, 1, 0))

#unite into global a and s
temps1 <- full_join(a549fs, full_join(aalebfs, arpefs))
temps2 <- full_join(fibroblastfs, full_join(hepatocytefs, keratinocytefs))
temps3 <- full_join(microglialfs, thyroidfs)
tempa1 <- full_join(a549fa, full_join(aalebfa, arpefa))
tempa2 <- full_join(fibroblastfa, full_join(hepatocytefa, keratinocytefa))
tempa3 <- full_join(microglialfa, thyroidfa)
globals <- full_join(temps1, full_join(temps2, temps3))
globala <- full_join(tempa1, full_join(tempa2, tempa3))

#fix duplicates and nas
globals <- globals[!is.na(globals$gene),]
globala <- globala[!is.na(globala$gene),]
globals.nodup <- globals[!duplicated(globals$gene),]
globala.nodup <- globala[!duplicated(globala$gene),]
globals.nodup[is.na(globals.nodup)] <- 0
globala.nodup[is.na(globala.nodup)] <- 0

#save
write.csv(globals.nodup, file = "~/Desktop/MacMicking Lab/IFN_microarray/suppressedmatrix.csv")
write.csv(globala.nodup, file = "~/Desktop/MacMicking Lab/IFN_microarray/activatedmatrix.csv")
