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

deseq_df

#order by pvalue to get most variable genes, not right though, how to get most variable genes?
deseq_df <- deseq_df[order(abs(-deseq_df$log2FoldChange), decreasing = TRUE),]
#get first 5000 and select only gene name and pval
dataset <- deseq_df[c(1,3)]
dataset
dataset <- dataset %>% 
  mutate(condition = 
           case_when(log2FoldChange < 0 ~ "inf",
                     log2FoldChange > 0 ~ "control"))
dataset
dataset[1:100, c(2)]

get_dataset <- function(index=5000) {
  toret <- dataset[1:index, c(2)]
  
  return (toret)
}

gene10 <- get_dataset(10)
gene100 <- get_dataset(100)
gene1000 <- get_dataset(1000)
gene10000 <- get_dataset(10000)
gene5000 <- get_dataset(5000)

#different k values
clusters5000 <-  kmeans(gene5000, 5, iter.max = 25, nstart = 1)
clusters5000$size
clusters5000 <-  kmeans(gene5000, 10, iter.max = 25, nstart = 1)
clusters5000$size
clusters5000 <-  kmeans(gene5000, 15, iter.max = 25, nstart = 1)
clusters5000$size
clusters5000 <-  kmeans(gene5000, 3, iter.max = 25, nstart = 1)
clusters5000$size

#different gene numbers, 10, 100, 1000, 10000, replace with your clustering algorithm
clusters10 <- kmeans(gene10, 5, iter.max = 25, nstart = 1)
clusters100 <- kmeans(gene100, 5, iter.max = 25, nstart = 1)
clusters1000 <- kmeans(gene1000, 5, iter.max = 25, nstart = 1)
clusters10000 <- kmeans(gene10000, 5, iter.max = 25, nstart = 1)

clusters10$size
clusters100$size
clusters1000$size
clusters10000$size

#install.packages("ggalluvial")
library("ggalluvial")

#creating sanky plot
all10 <- c( clusters10$size)
all100 <- c( clusters100$size)
all1000 <- c( clusters1000$size)
all10000 <- c( clusters10000$size)
clusternum <- c( "1", "2", "3", "4", "5")
type <- c("gene10", "gene10","gene10","gene10","gene10",
          "gene100","gene100","gene100","gene100","gene100",
          "gene1000","gene1000","gene1000","gene1000","gene1000",
          "gene10000","gene10000","gene10000","gene10000","gene10000")
scalar1 <- function(x) {x / sqrt(sum(x^2))}
scalar1(clusters10$size)
normalize(clusters10$size)
freq <- c(clusters10$size/10, clusters100$size/100,clusters1000$size/1000, clusters10000$size/10000)

#k = 5
group <- c("1", "2", "3", "4", "5",
           "1", "2", "3", "4", "5",
           "1", "2", "3", "4", "5",
           "1", "2", "3", "4", "5")

#plotting
toplot <- data.frame(type, freq, group)

toplot

p <- ggplot(data = toplot,mapping= aes(x = type,stratum = group, alluvium = group, y = freq, fill = group))+
  geom_flow(stat = "alluvium") +
  geom_stratum(alpha = .5) +
  scale_fill_manual(values = c("grey", "green", "red", "blue", "black"))  +
  ggtitle("Changes in group memberships for different number of genes")
show(p)


