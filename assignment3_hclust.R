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
deseq_df <- deseq_df[order(abs(-deseq_df$padj), decreasing = FALSE),]
#get first 5000 and select only gene name and pval
dataset <- deseq_df[3]
dataset
dataset <- dataset %>% 
  mutate(condition = 
           case_when(log2FoldChange < 0 ~ "inf",
                     log2FoldChange > 0 ~ "control"))
dataset
dataset[1:100, c(1)]

get_gene_data <- function(genes=5000) {
  return( gene_matrix[ deseq_df[1:genes, "Gene" ], ] )
}

gene10 <- get_gene_data(10)
gene100 <- get_gene_data(100)
gene1000 <- get_gene_data(1000)
gene10000 <- get_gene_data(10000)


###Using euclidean algorithm for determined distance
dist_gene10 <- dist(scale(t(gene10)), method = "euclidean")
dist_gene100 <- dist(scale(t(gene100)), method = "euclidean")
dist_gene1000 <- dist(scale(t(gene1000)), method = "euclidean")
dist_gene10000 <- dist(scale(t(gene10000)), method = "euclidean")


#ward method
hclust_gene10 <- hclust(dist_gene10, method = "ward.D")
plot(hclust_gene10)

hclust_gene100 <- hclust(dist_gene100, method = "ward.D")
plot(hclust_gene100)

hclust_gene1000 <- hclust(dist_gene1000, method = "ward.D")
plot(hclust_gene1000)

hclust_gene10000 <- hclust(dist_gene10000, method = "ward.D")
plot(hclust_gene10000)


#complete method
hclust_gene10 <- hclust(dist_gene10, method = "complete")
plot(hclust_gene10)

hclust_gene100 <- hclust(dist_gene100, method = "complete")
plot(hclust_gene100)

hclust_gene1000 <- hclust(dist_gene1000, method = "complete")
plot(hclust_gene1000)

hclust_gene10000 <- hclust(dist_gene10000, method = "complete")
plot(hclust_gene10000)



#install.packages("ggalluvial")
library("ggalluvial")

#creating sankey plot
all10 <- c( hclust_gene10$size)
all100 <- c( hclust_gene100$size)
all1000 <- c( hclust_gene1000$size)
all10000 <- c( hclust_gene1000$size)
clusternum <- c( "1", "2", "3", "4", "5")
type <- c("gene10", "gene10","gene10","gene10","gene10",
          "gene100","gene100","gene100","gene100","gene100",
          "gene1000","gene1000","gene1000","gene1000","gene1000",
          "gene10000","gene10000","gene10000","gene10000","gene10000")
scalar1 <- function(x) {x / sqrt(sum(x^2))}
scalar1(hclust_gene10$size)
normalize(hclust_gene10$size)
freq <- c(hclust_gene10$size/10, hclust_gene100$size/100,hclust_gene1000$size/1000, hclust_gene10000$size/10000)

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



##heatmap
heatmap(as.matrix(gene100), scale="row")


##chi-squared test
#can test different values of k
cutree(hclust_gene1000, k=6)
data.frame(clust_assign = cutree(hclust_gene1000, k=6) )
cluster_info <- data.frame(clust_assign = cutree(hclust_gene1000, k=6) )
cluster_info$status <- "inf"
cluster_info$status[14:-1] <- "ctrl"
cluster_info
chisq.test( table(cluster_info$clust_assign, cluster_info$status) )



