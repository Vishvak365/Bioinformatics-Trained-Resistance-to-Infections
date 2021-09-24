if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
install.packages("readxl")
install.packages("tibble")
install.packages("readr")
install.packages("dplyr")
#loading data
library("readxl")
library("DESeq2")
library("tibble")
library("dplyr")

expression_df <- readr::read_tsv("RNA_counts_orig.tsv") %>%
  tibble::column_to_rownames("Gene")
expression_df

filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(countData = gene_matrix, colData = infected_data, design = ~infected_data )
deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
head(deseq_results)



#-----OLD-------#



#getting orig RNA counts and forcing data types
RNA_counts_orig <- read_xlsx("RNA_counts_orig.xlsx")
samp2 <- RNA_counts_orig[,-1]
rownames(samp2) <- RNA_counts_orig[,1]
#RNA_counts_orig %>% remove_rownames %>% column_to_rownames(var="Gene")
RNA_counts_orig


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
ddset <- DESeqDataSetFromMatrix(countData = expression_df, colData = infected_data, design = infected_data )
?DESeqDataSetFromMatrix
