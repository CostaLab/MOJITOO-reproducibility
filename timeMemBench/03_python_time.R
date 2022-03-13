mute <- suppressPackageStartupMessages
mute(library(MOFA2))
mute(library(MOJITOO))
mute(library(Signac))
mute(library(rliger))
mute(library(mclust))
mute(library(glue))
mute(library(FNN))
mute(library(RANN)) ## nn2
mute(library(dplyr))
mute(library(reshape2))
mute(library(cli))
mute(library(Rcpp))
mute(library(RcppArmadillo))
mute(library(cluster))
mute(library(parallelDist))
mute(library(S4Vectors))
mute(library(Matrix))
mute(library(assertthat))
mute(library(SeuratDisk))
mute(library(knn.covertree))
mute(library(Seurat))
mute(library(pls))
mute(library(optparse))
mute(library(CCA))
mute(library(Signac))
mute(library(MOJITOO))
mute(library(dplyr))
mute(library(mclust))
mute(library(data.table))
mute(library(bench))
mute(library(MOFA2))
mute(library(rbenchmark))



source("../util/util.R")


ASSAY_1 = "RNA"
ASSAY_2 = "ATAC"


nums <- (1:10)*3000 

differ_list <- list()

for(num in nums){
  t1 <- Sys.time()  
  system(glue("python 00_schema_prepare.py -d {num}"))
  t2 <- Sys.time()
  differ_list[[as.character(num)]] <-(difftime(t2, t1, units = "secs")[[1]]) 
}
saveRDS(differ_list, "save/schema_list.Rds")


