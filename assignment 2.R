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

#plotPCA of summed
plotPCA(DESeqTransform(summed), intgroup = c("X"))


library(M3C)
#tsne? 2 diff graphs and compare? why do they change every time???
tsne(data.matrix(infected),perplex = 3, labels = c("infected"))
tsne(data.matrix(control),perplex = 3, labels = c("control"))

#differential analysis...
#prep data...