#loading data
library("readxl")

#getting orig RNA counts and forcing data types
RNA_counts_orig <- read_xlsx("RNA_counts_orig.xlsx", col_types=c("text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                                 "numeric","numeric","numeric","numeric","numeric", "numeric"))


#making it a dataframe and some renaming
RNA_counts_orig_dataframe <- as.data.frame(RNA_counts_orig)
names(RNA_counts_orig_dataframe)[1] <- "condition"

#summarizing data, each col needs a colData which is different from just the name of the column, its like a hidden label for each col
summed <- SummarizedExperiment(RNA_counts_orig_dataframe[2:17977,2:27 ], colData= c("infected", "infected","infected", "infected","infected", "infected",
                                                                                    "infected", "infected","infected", "infected","infected", "infected",
                                                                                    "infected", "control", "control","control", "control",
                                                                                    "control", "control","control", "control","control", "control",
                                                                                    "control", "control", "control"))




#compute gene expression ranges and plot density of ranges
library(matrixStats)
numeric_matrix <- as.matrix(sapply(RNA_counts_orig[1:17977,2:27 ], as.numeric))
gene_minmax <- rowRanges(numeric_matrix)
gene_minmax
gene_ranges <- gene_minmax[,2] - gene_minmax[,1]
length(gene_ranges)
gene_ranges


#density plots
plot(density(gene_ranges))
plot(density(log2(gene_ranges)))

#plotPCA of summed
plotPCA(DESeqTransform(summed), intgroup = c("X"))

#Tsne plot
library(M3C)
tsne(RNA_counts_orig[,2:27], perplex=1, labels = c("infected", "infected","infected", "infected","infected", "infected",
                                                   "infected", "infected","infected", "infected","infected", "infected",
                                                   "infected", "control", "control","control", "control",
                                                   "control", "control","control", "control","control", "control",
                                                   "control", "control", "control"))

#Volcano plot
