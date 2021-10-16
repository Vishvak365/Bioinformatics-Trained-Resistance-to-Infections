if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("EnhancedVolcano")

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
library("apeglm")
?read_tsv
expression_df <- readr::read_tsv("RNA_counts_orig.tsv",col_types = c("c","d","d","d","d","d","d","d","d","d","d",
                                                                     "d","d","d","d","d","d","d","d","d","d",
                                                                     "d","d","d","d","d", "d")) %>%  tibble::column_to_rownames("Gene")
#?column_to_rownames
#rownames(expression_df) = make.names(expression_df[1], unique=TRUE)
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

gene_matrix <- round(filtered_expression_df)


inf_data <- c("infected", "infected","infected", "infected","infected", "infected",
              "infected", "infected","infected", "infected","infected", "infected",
              "infected", "control", "control","control", "control",
              "control", "control","control", "control","control", "control",
              "control", "control", "control")
infected_data <- data.frame("infected_data"=inf_data)

ddset <- DESeqDataSetFromMatrix(countData = gene_matrix, colData = infected_data, design = ~infected_data )
deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
deseq_results



deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))


#order by pvalue to get most variable genes, not right though, how to get most variable genes?
deseq_df[order(-deseq_df$baseMean),]
deseq_df <- deseq_df[order(-deseq_df$baseMean),]
#get first 5000 and select only gene name and pval
dataset <- deseq_df[1:5000,c (1, 2)]

dataset

ran <- sample(1:nrow(dataset), .9* nrow(dataset))
nor <- function(x) {(x-min(x))/(max(x)-min(x))}

datasetnorm <- as.data.frame(lapply(dataset[2], nor))
summary(datasetnorm)
datasetnorm

dataset_train <- datasetnorm[ran,1]
dataset_test <- datasetnorm[-ran,1]

dataset_target_cat <- dataset[ran,1]
dataset_test_cat <- dataset[-ran, 1]

dataset_test_cat
library("class")
as.factor(dataset_train)
pr <- knn(data.frame(dataset_train), data.frame(dataset_test), cl = as.factor(dataset_target_cat), k = 70)

tab <- table(pr, dataset_test_cat)

accuracy <- function(x) {sum(diag(x))/(sum(rowSums(x))) * 100}
accuracy(tab)

