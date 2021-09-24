if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
install.packages("readxl")
install.packages("tibble")
install.packages("readr")

#loading data
library("readxl")
library("DESeq2")
#getting orig RNA counts and forcing data types
RNA_counts_orig <- read_xlsx("RNA_counts_orig.xlsx", rowNames=TRUE)
read.cs
?read_xlsx


#making it a dataframe and some renaming
RNA_counts_orig_dataframe <- as.data.frame(RNA_counts_orig)
names(RNA_counts_orig_dataframe)[1] <- "condition"
RNA_counts_orig_dataframe
columnData <- c("text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                "numeric","numeric","numeric","numeric","numeric", "numeric")

columnData
inf_data <- c("infected", "infected","infected", "infected","infected", "infected",
              "infected", "infected","infected", "infected","infected", "infected",
              "infected", "control", "control","control", "control",
              "control", "control","control", "control","control", "control",
              "control", "control", "control")
infected_data <- data.frame("infected_data"=inf_data)
table(infected_data)
dim(RNA_counts_orig_dataframe)
ncol(RNA_counts_orig_dataframe)
ddset <- DESeqDataSetFromMatrix(countData = RNA_counts_orig_dataframe, colData = infected_data, design = infected_data )
?DESeqDataSetFromMatrix
