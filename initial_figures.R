library("biomaRt")
library("dplyr")
library("org.Hs.eg.db")
library("topGO")
library(reshape2)
library(ggplot2)

ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

activated <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/activatedmatrix.csv")[,c(2:10)]
suppressed <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/suppressedmatrix.csv")[,c(2:10)]

row.names(activated) <- activated[,1]
row.names(suppressed) <- suppressed[,1]
activated <- activated[,-1]
suppressed <- suppressed[,-1]

rs_a <- data.frame(rowSums(activated), row.names = row.names(activated))
rs_s <- data.frame(rowSums(suppressed), row.names = row.names(suppressed))

a_t <- as.data.frame(table(rs_a))[-1,]
a_t$Freq <- a_t$Freq/sum(c(a_t$Freq))*100
s_t <- as.data.frame(table(rs_s))[-1,]
s_t$Freq <- s_t$Freq/sum(c(s_t$Freq))*100
colnames(s_t) <- c('cell_number', 'suppressed')
colnames(a_t) <- c('cell_number', 'activated')

t <- left_join(a_t, s_t, by = 'cell_number')
t[is.na(t)] <- 0

melted <- melt(t)

barplot <- ggplot(melted, aes(fill=variable, y=value, x=cell_number)) + geom_bar(position="dodge", stat="identity")
barplot + ggtitle('Fraction of Stimulated Genes by Cell Type Occurrence') + 
  xlab('Number of occurrences across cell types') + ylab('Percent of stimulated genes') + 
  scale_fill_discrete(name = 'State', labels = c('Activated', 'Suppressed'))
ggsave(filename = "initial.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")

#boxplot
rs_a_b <- data.frame(Sum=rowSums(activated), State='activated')
rs_s_b <- data.frame(Sum=rowSums(suppressed), State='suppressed')
rs_b <- rbind(rs_a_b, rs_s_b)
rs_b <- rs_b[rs_b$Sum != 0,]

melted$cell_number <- as.numeric(melted$cell_number)
melted$box <- (melted$value*melted$cell_number)/sum(melted$value*melted$cell_number)*100

bplot <- ggplot(rs_b, aes(fill = State, x=State, y=Sum)) + geom_boxplot()
bplot + labs(title = 'Gene Cell Type Occurrence by Stimulation State', x = 'State', 
             y = 'Gene Cell Type Frequency') + guides(fill=FALSE)
ggsave(filename = "initial_boxplot.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")

stat_dist <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/Processed Data/subgraph_quant_data_wund.csv")
colnames(stat_dist) <- c('Method 1', 'Method 2')
melted_stat <- melt(stat_dist)

ggplot(melted_stat,aes(x=value)) + geom_histogram(alpha=1, bins=50) + facet_grid(~variable) +
  ggtitle('Distributions of Statistic Values according to Method') + 
  xlab('Value') + ylab('Count')
ggsave(filename = "stat_hist.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")

ggplot(melted_stat,aes(x=variable, y=value, fill=variable)) + geom_boxplot() + 
  ggtitle('Box Plot of Cluster Stimulation Measurements') + 
  xlab('Method') + ylab('Values') + guides(fill=FALSE)
ggsave(filename = "initial_boxplot.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")


full_data <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/bigdata.csv")
full_data <- full_data[full_data$stimul != 0,]
state_data <- read.csv("~/Desktop/MacMicking Lab/IFN_microarray/Processed Data/boxplots_data.csv")

c_cen_s <- state_data[state_data$label == 'suppressed',]
c_cen_a <- state_data[state_data$label == 'activated',]
c_cen_s <- c_cen_s[,5]
c_cen_a <- c_cen_a[,5]

#boxplots of centrality measures according to state does not work
boxplot_data <- data.frame(closeness=c(c_cen_s, c_cen_a), 
                           state = factor(rep(c("Suppressed", "Activated"), 
                                              times = c(length(c_cen_s), length(c_cen_a)))))

bplot <- ggplot(boxplot_data, aes(fill = state, x=state, y=closeness)) + geom_boxplot()
bplot + labs(title = 'Degree Centrality by Stimulation State', x = 'State', 
             y = 'Degree Centrality of nodes') + guides(fill=FALSE)
ggsave(filename = "d_cen_boxplot.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")

#scatter plots of stimulation and centralities
boxPlot <- ggplot(full_data, aes(x=full_data$mISG, y=full_data$d_cen, group=full_data$mISG)) + geom_boxplot()
boxPlot + labs(title = 'Degree Centrality according to Stimulation Measure', x = 'Stimulation', 
               y = 'Degree Centrality of nodes') + guides(fill=FALSE)
ggsave(filename = "d_cen_misg.jpg",
       path="~/Desktop/MacMicking Lab/IFN_microarray/Figures/")



