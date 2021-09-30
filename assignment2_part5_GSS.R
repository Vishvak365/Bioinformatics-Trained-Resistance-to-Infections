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
