if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("EnhancedVolcano")
BiocManager::install("ComplexHeatmap")

install.packages("readxl")
install.packages("tibble")
install.packages("readr")
install.packages("dplyr")
BiocManager::install("ComplexHeatmap", force=TRUE)
#loading data
library("readxl")
library("DESeq2")
library("tibble")
library("dplyr")
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

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
volcano_plot
ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "volcano_plot.png")
)



#-----Part 4--------

#only 1 significant gene, uncomment for actual threshold
#heatmap_genes = deseq_df$Gene[deseq_df$threshold==TRUE] 

#top ten significant genes
heatmap_genes = deseq_df$Gene[1:10]
heatmap_genes %in% rownames(ddset@assays@data@listData[["counts"]])
#subset the counts that match the top ten from full data set
heatmap_expression_data = ddset@assays@data@listData[["counts"]][rownames(ddset@assays@data@listData[["counts"]]) %in% heatmap_genes,]

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
dev.off()
