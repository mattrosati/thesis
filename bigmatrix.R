library(plyr)
library(dplyr)


for(i in c(1:length(data_list))) {
  assign(names[i], read_csv(data_list[i]))
}

data_list <- list.files("~/Desktop/MacMicking Lab/IFN_microarray/Processed Data",full.names = T, pattern="\\.csv")
names <- c("gene", "arpe", "fibroblast", "hepatocyte", "hsg", "keratinocyte", "migroglial")
state <- c("global", "activated", "suppressed")

#read files
data <- lapply(data_list, read_csv)

#eliminate missing genes
data.nona <- lapply(data, function(x) x[!is.na(x$Gene.symbol),])

#converts all non-significant logFC to zero, keeping significant logFC (fdr of 0.1)
#NOTE: coercing symbol to factor
data.converted <- lapply(data.nona, function(x) data.frame(gene = as.character(x$Gene.symbol), 
                                                      logFC = ifelse(x$adj.P.Val < 0.1, x$logFC, 0)))

#separate into a and s
data.activated <- lapply(data.converted, function(x) data.frame(gene = x$gene, 
                                                                logFC = ifelse(x$logFC > 1, x$logFC, 0)))
data.suppressed <- lapply(data.converted, function(x) data.frame(gene = x$gene, 
                                                                logFC = ifelse(x$logFC < -1, x$logFC, 0)))

#unite into a global list of activated, suppressed, and general
matrices <- list(data.converted, data.activated, data.suppressed)
all_global <- lapply(matrices, function(x) join_all(x, by="gene", type = "inner"))

#change column names
all_global <- lapply(all_global, setNames, names)


#save
for(i in c(1:length(state))) {
  path <- paste("~/Desktop/MacMicking Lab/IFN_microarray/", state[i], ".csv", sep = "")
  entry <- all_global[i]
  write.csv(entry, file = path)
}
