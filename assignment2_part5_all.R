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

########Part 5i ClusterProfiler (Gene Ontology) ####################
BiocManager::install("org.Rn.eg.db")
BiocManager::install("clusterProfiler")
library("org.Rn.eg.db")
library("clusterProfiler")

#sorting deseq_df's p values in decending order (unnecessary)
deseq_results_sorted <- deseq_df[order(-deseq_df$pvalue),]
deseq_results_sorted
#making a new dataframe with just gene name and sorted p values (unnecessary, just need gene names)
df <- deseq_results_sorted[,c(1,5)]


#library 
BiocManager::install("topGO")
library("topGO")


#getting entrez ids from gene name list, not sure if correct
entrez <- bitr(df$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Rn.eg.db")
entrez
#needs entrez ids for gene

enriched <- enrichGO(entrez$ENTREZID, 
         'org.Rn.eg.db',
         minGSSize = 1,
         ont = "CC",
         pvalueCutoff = .01)

enriched
head(enriched)
###################################################################


######################part 5ii #####################################
#Gprofiler2 gene ontology
install.packages("gprofiler2")
library("gprofiler2")

gostres <- gost(entrez$SYMBOL, organism="mmusculus", user_threshold=0.05)

#result dataframe
gostres$result

######################################################################



########Part 5iii clusterProfiler for Disease Ontology####################
BiocManager::install("DOSE")
library("DOSE")


#getting entrez ids and enriching with disease ontology
enriched_DO <- enrichDO(entrez$ENTREZID, 
                        ont = "DO",
                        minGSSize = 1,
                        pAdjustMethod = "bonferroni",
                        pvalueCutoff = .01)
#can adjust pvalue, but higher pvalue still produces null result

enriched_DO


####################################################################

########Part 5iv topGO for Gene Ontology####################

BiocManager::install("GenomicSuperSignature")
deseq_df
browseVignettes("GenomicSuperSignature")
BiocManager::install("topGO")
BiocManager::install("ALL")
BiocManager::install("affyLib")
BiocManager::install("hgu95av2.db")


library(topGO)
library(ALL)
data(ALL)
data(geneList)
geneList
type(deseq_df['pvalue'])

affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)
sum(topDiffGenes(deseq_df['pvalue']))

GOdata <- new("topGOdata", description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)
GOdata
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultKS
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes
pValue.classic <- score(resultKS)
pValue.classic


