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
dataset <- deseq_df[3]
dataset
dataset <- dataset %>% 
  mutate(condition = 
           case_when(log2FoldChange < 0 ~ "inf",
                     log2FoldChange > 0 ~ "control"))
dataset
dataset[1:100, c(1)]

get_dataset <- function(index=5000) {
  toret <- dataset[1:index, c(1)]
  
  return (toret)
}

gene10 <- get_dataset(10)
gene100 <- get_dataset(100)
gene1000 <- get_dataset(1000)
gene10000 <- get_dataset(10000)

clusters10 <- kmeans(gene10, 5, iter.max = 25, nstart = 1)
clusters100 <- kmeans(gene100, 5, iter.max = 25, nstart = 1)
clusters1000 <- kmeans(gene1000, 5, iter.max = 25, nstart = 1)
clusters10000 <- kmeans(gene10000, 5, iter.max = 25, nstart = 1)

#install.packages("ggalluvial")
library("ggalluvial")


clusters10

ggplot(data = as.data.frame(clusters),
       aes(axis1 = gene_Num, y = accuracy/3)) +
  scale_x_discrete(limits = c("gene_Num"), expand = c(.02, .05)) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  ggtitle("Accuracy KNN on different number of genes")

clusters$size

#top ten significant genes
heatmap_genes = deseq_df$log2FoldChange[1:10]
heatmap_genes %in% rownames(ddset@assays@data@listData[["counts"]])
#subset the counts that match the top ten from full data set
heatmap_expression_data = ddset@assays@data@listData[["counts"]][rownames(ddset@assays@data@listData[["counts"]]) %in% heatmap_genes,]

install.packages("ComplexHeatmap")

library(ComplexHeatmap)
#set colors of two groups- infected and not infected
col_fun = colorRamp2(c(0, 1), c("black", "white"))
#generate heatmap and save to png
png(filename = "gene_expr_heatmap.png",
    width = 10, height = 8)
Heatmap(heatmap_expression_data, bottom_annotation = HeatmapAnnotation(sample_type = ifelse(inf_data=="infected",
                                                                                            0,
                                                                                            1),
                                                                       col = list(sample_type = col_fun),
                                                                       show_annotation_name = FALSE),
        name = "Gene\nExpression")
