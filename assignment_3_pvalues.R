# inf_data[1:13]
# inf_data[14:26]
#infected_data <- data.frame("infected_data"=inf_data[14:26])
#ddset <- DESeqDataSetFromMatrix(countData = gene_matrix[14:26], colData = infected_data, design = ~1)
library("MASS")
infected_data <- data.frame("infected_data"=inf_data[1:13])
ddset <- DESeqDataSetFromMatrix(countData = gene_matrix[1:13], colData = infected_data, design = ~1)
deseq_object <- DESeq(ddset)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

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
deseq_df <- deseq_df[order(abs(-deseq_df$log2FoldChange), decreasing = TRUE),]
#get first 5000 and select only gene name and pval
dataset <- deseq_df[c(1,3)]
dataset <- dataset %>% 
  mutate(condition = 
           case_when(log2FoldChange < 0 ~ "inf",
                     log2FoldChange > 0 ~ "control"))
dataset

get_dataset_compare <- function(subset=5000){
  #Below is the data partitioning steps
  dataset <- deseq_df[1:subset,c (1, 2)]
  nor <- function(x) {(x-min(x))/(max(x)-min(x))}
  datasetnorm <- as.data.frame(lapply(dataset[2], nor))
  return(datasetnorm)
}
data <- get_dataset_compare()
classes_gmm <- predict(fit,data)$classification
as.data.frame(table(classes_gmm)) 

control_gmm_cluster <- c(3720,1041,239)
infected_gmm_cluster <- c(3737,1028,235)

gmm_cluster_all <- c(3097, 1594, 309)
kmeans_cluster_all <- c(4999,1,0,0)
consensus_cluster_all <- c(286,286,260,234,208,182,156,130,104,78,52,26)

chisq.test(control_gmm_cluster,infected_gmm_cluster)
chisq.test(kmeans_cluster_all,gmm_cluster_all)
chisq.test(kmeans_cluster_all[1:4],consensus_cluster_all[1:4])
chisq.test(consensus_cluster_all[1:3],gmm_cluster_all[1:3])

p_values <- c(.2231,.1911,.2231)
p.adjust(p_values, method = "holm", n = length(p_values))
