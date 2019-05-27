library("biomaRt")
library("dplyr")
library("org.Hs.eg.db")
library("topGO")
library(reshape2)
library(ggplot2)

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
