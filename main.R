#import relevant functions and paackages

source("packages.R")
source("differential_expression.R")

#all sources dataframe
# I like: Microglial, ARPE, Fibroblast (right skewed intensities), 
#         Hepatocyte (right skewed intensities), Keratinocyte, ~HSG
# AALEB is a problem because of extremely skewed intensities, 
# Thyroid just doesn't have any matches?

all_series <- data.frame(Cell = c("Thyroid", "Microglial", "AALEB", "ARPE", "Fibroblast", "Hepatocyte", "Keratinocyte", "HSG"),
                         Series = c("GSE5054", "GSE1432", "GSE77154", "GSE36331", "GSE67737", "GSE38147", "GSE12109", "GSE41291"),
                         Platform = c("GPL96", "GPL96", "GPL570", "GPL6244", "GPL10558", "GPL6244", "GPL571", "GPL6947"),
                         Selections = c("01XX01XX01XX", "XXXXXXXX1111XXXXXXXX0000", "0011XXXX", "00011XXXXXXXXXXXX", "000XXXXXXXXXXXXXXX111XXX", "00XXXXXX11", "1XX01XX01XX0",
                                        "XXXXXXXXXXXXXXXXXX000XXX111XXXXXXXXX"))
i <- "Microglial"
#loop over all cell types and perform differential expression
#for (i in all_series$Cell) {

#load into environment
gset <- load(i)

#need: quality control 
#      probe filtering if multiple mappings
#      probe filtering if no mappings

#check low intensity probes
tran <- log2(exprs(gset))
medians <- rowMedians(Biobase::exprs(gset))
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#determines design matrix with blocking
design <- blocking(gset)

#returns toptable
tT <- p.vals(gset, design)

#some interesting stats are printed
of_interest <- as.integer(nrow(tT[tT$adj.P.Val <= 0.1, ]))
supp <- as.integer(nrow(tT[tT$adj.P.Val <= 0.1 & tT$logFC < 0, ]))
act <- as.integer(nrow(tT[tT$adj.P.Val <= 0.1 & tT$logFC > 0, ]))

print("Number of genes with significant adjusted p values: ", of_interest)
print("Number of activated genes: ", act)
print("Number of suppressed genes: ", supp)
hist(tT$P.Value)

#}

