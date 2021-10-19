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

