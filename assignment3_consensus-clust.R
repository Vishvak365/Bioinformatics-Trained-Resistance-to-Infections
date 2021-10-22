#LOAD AND SET DIFFERENTIALLY EXPRESSED GENES
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


#extract 5000 most variable based on log2fold change
#order genes by variability
most_variable_genes <- deseq_df[order(-deseq_df$log2FoldChange),]
#take top 5000
var_genes <- head(most_variable_genes, 5004)


#RUN CONSENSUS CLUSTERING
#BiocManager::install("ConsensusClusterPlus")
library("ConsensusClusterPlus")

#getting orig RNA counts and forcing data types
RNA_counts_orig <- read_xlsx("RNA_counts_orig.xlsx", col_types=c("text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                                 "numeric","numeric","numeric","numeric","numeric", "numeric"))
length(var_genes$Gene)

#transpose so that we can index columns by gene name
counts_transposed <- data.frame(t(data.frame(RNA_counts_orig)[-1]))
colnames(counts_transposed) <- t(RNA_counts_orig[,1])

colnames(counts_transposed)

#remove undefined names
undefined_ct <- 0
for (name in var_genes$Gene){
  if(!(name %in% colnames(counts_transposed))){
    print(name)
    undefined_ct <- undefined_ct + 1
  }
}
undefined_ct


proc <- var_genes$Gene[!var_genes$Gene %in% c("2-Mar","5-Mar", "7-Sep", "6-Sep")]
length(proc)

#get columns of variable genes from original data
clust_mat <- counts_transposed[,proc]
#transpose to conform with concensus clustering expectations
#columns are samples, rows are features (genes)
clust_t <- t(data.frame(clust_mat))

#median center normalization
clust_t = sweep(clust_t,1, apply(clust_t,1,median,na.rm=T))

results = ConsensusClusterPlus(clust_t,maxK=10,reps=50,pItem=0.8,pFeature=1,title=tempdir(),clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot=NULL)
icl = calcICL(results,title=tempdir(),plot=NULL)
#icl[["itemConsensus"]]
#5000
as.data.frame(table(icl$itemConsensus['cluster']))


(most_variable_genes[1:2000,])

head(most_variable_genes$Gene,20)

run_clustering <- function(num_genes){
  vars <- head(most_variable_genes$Gene, num_genes+4)
  proc <- vars[!vars %in% c("2-Mar","5-Mar", "7-Sep", "6-Sep")]
  clust_mat <- counts_transposed[,proc]
  #transpose to conform with concensus clustering expectations
  #columns are samples, rows are features (genes)
  clust_t <- t(data.frame(clust_mat))
  
  #median center normalization
  clust_t = sweep(clust_t,1, apply(clust_t,1,median,na.rm=T))
  
  results = ConsensusClusterPlus(clust_t,maxK=5,reps=50,pItem=0.8,pFeature=1,title=tempdir(),clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot=NULL)
  icl = calcICL(results,title=tempdir(),plot=NULL)

  return (as.data.frame(table(icl$itemConsensus['cluster'])))
}

cl10 <- run_clustering(10)
cl100 <- run_clustering(100)
cl1000 <- run_clustering(1000)
cl10000 <- run_clustering(10000)

scale <- function(x){
  return (as.double((x-min(x)))/as.double((max(x)-min(x))))
}

#freq <- c(scale(cl10$Freq), scale(cl100$Freq),scale(cl1000$Freq), scale(cl10000$Freq))
group <- c("1", "2", "3", "4", "5",
           "1", "2", "3", "4", "5",
           "1", "2", "3", "4", "5",
           "1", "2", "3", "4", "5")
type <- c("gene10", "gene10","gene10","gene10","gene10",
          "gene100","gene100","gene100","gene100","gene100",
          "gene1000","gene1000","gene1000","gene1000","gene1000",
          "gene10000","gene10000","gene10000","gene10000","gene10000")
toplot <- data.frame(type, freq, group)

toplot

library("ggplot2")
library("ggalluvial")
p <- ggplot(data = toplot,mapping= aes(x = type,stratum = group, alluvium = group, y = freq, fill = group))+
  geom_flow(stat = "alluvium") +
  geom_stratum(alpha = .5) +
  scale_fill_manual(values = c("grey", "green", "red", "blue", "black"))  +
  ggtitle("Changes in group memberships for different number of genes - Consensus Clustering")
show(p)

freq

