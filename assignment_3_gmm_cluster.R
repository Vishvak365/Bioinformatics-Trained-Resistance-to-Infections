# !Run the contents of the 'assignment 3.R' file first in R Studio before running the below lines
library(mclust, quietly=TRUE)

get_dataset <- function(subset=5000){
  #Below is the data partitioning steps
  dataset <- deseq_df[1:subset,c (1, 2)]
  nor <- function(x) {(x-min(x))/(max(x)-min(x))}
  datasetnorm <- as.data.frame(lapply(dataset[2], nor))
  return(datasetnorm)
}
sub_5000 <- get_dataset(subset = 5000)
sub_10 <- get_dataset(subset = 10)
sub_100 <- get_dataset(subset = 100)
sub_1000 <- get_dataset(subset = 1000)
sub_10000 <- get_dataset(subset = 10000)

fit = Mclust(sub_10000, G=3)
plot(fit, what="density", main="",xlim=c(0,.02))
summary(fit)
clust_5000 <- c(3097, 1594, 309)
clust_10 <- c(6, 2, 2)
clust_100 <- c(48, 44, 8)
clust_1000 <- c(612, 336, 52) 
clust_10000 <- c(6384, 3130, 486)

all10 <- c( clust10)
all100 <- c( clust100)
all1000 <- c( clust1000)
all10000 <- c( clust10000)
clusternum <- c( "1", "2", "3")
type <- c("gene10", "gene10","gene10",
          "gene100","gene100","gene100",
          "gene1000","gene1000","gene1000",
          "gene10000","gene10000","gene10000")
scalar1 <- function(x) {x / sqrt(sum(x^2))}
scalar1(clusters10$size)
normalize(clusters10$size)
freq <- c(clusters10$size/10, clusters100$size/100,clusters1000$size/1000, clusters10000$size/10000)

#k = 3
group <- c("1", "2", "3",
           "1", "2", "3",
           "1", "2", "3",
           "1", "2", "3",)

#plotting
toplot <- data.frame(type, freq, group)

toplot

p <- ggplot(data = toplot,mapping= aes(x = type,stratum = group, alluvium = group, y = freq, fill = group))+
  geom_flow(stat = "alluvium") +
  geom_stratum(alpha = .5) +
  scale_fill_manual(values = c("grey", "green", "red"))  +
  ggtitle("Changes in group memberships for different number of genes")
show(p)