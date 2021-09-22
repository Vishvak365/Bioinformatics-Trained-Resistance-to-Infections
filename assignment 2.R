#loading data
library("readxl")

#getting orig RNA counts
RNA_counts_orig <- read_xlsx("RNA_counts_orig.xlsx")


RNA_counts_trans <- as.data.frame(t(RNA_counts_orig))
plot(density(RNA_counts_orig[1,2: ]))
RNA_counts_orig_dataframe <- as.data.frame(RNA_counts_orig)
summed <- SummarizedExperiment(RNA_counts_orig_dataframe[2:17977,2:27 ])
plotPCA(DESeqTransform(summed))
colData(summed)

BiocManager::install("DESeq2")
install.packages("BiocManager")
library(M3C)

