mute <- suppressPackageStartupMessages
mute(library(S4Vectors))
mute(library(Seurat))
mute(library(glue))
mute(library(cli))
mute(library(optparse))

source("../util/util_tri.R")

AllOptions <- function(){
  parser <- OptionParser()
  parser <- add_option(parser, c("-d", "--dir"), type="character", default="ASAP_BM",
                       help="directory [default %default]",
                       metavar="character")
  return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)


dataset <- pa$dir
cf <- getSettings(dataset = dataset)

reduction_dims <- SimpleList(
                             "A1"            = cf$dims_1,
                             "A2"            = cf$dims_2,
                             "schema"        = 1:30,
                             "MOJITOO"       = NULL,
                             "WNN"           = NULL,
                             "MOFARed"       = NULL,
                             "DIABLORed"     = NULL
)



nms <- names(reduction_dims)
nms[1] <- cf$Redu_1
nms[2] <- cf$Redu_2
names(reduction_dims) <- nms


metrics_vec <- c(
                 "structure"         =  1,
                 "silhouette"        =  1,
                 "rangeclustering"   =  1
)


CALC_DIM12 = TRUE
if(CALC_DIM12){
  calc_dims <- 1:length(reduction_dims)
}else{
  calc_dims <- 3:length(reduction_dims)
}


multiome <- load_object(file=file.path(dataset, "/save/multiome.Rds"))


## set dimensions
for(redu in c("MOJITOO", "DIABLORed",  "MOFA", "MOFARed")){
  reduction_dims[[redu]] <- 1:ncol(multiome[[redu]])
}


##This avoid clustering using wrong graph
rednames <- names(multiome@reductions)
sapply(rednames, function(x) {multiome@reductions[[x]]@assay.used <- cf$ASSAY_1})


if(metrics_vec["structure"]  == 1 ){
  message(date(), " preparing structure")
  for(nm in names(reduction_dims)[3:length(reduction_dims)]){
    message("        ", date(), " ", nm)
    multiome <- getStructure(multiome,
                             nredu_1 = cf$Redu_1,
                             nredu_2 = cf$Redu_2,
                             nredu_3 = cf$Redu_3,
                             dims_1 = cf$dims_1,
                             dims_2 = cf$dims_2,
                             dims_3 = cf$dims_3,
                             nredu=nm, 
                             dims=reduction_dims[[nm]],
                             cor_method = "pearson")
  }
}

if(metrics_vec["silhouette"]  == 1 ){
  ## get Silhouette
  message(date(), " preparing silhouette")
  for(nm in names(reduction_dims)[calc_dims]){
    message("        ", date(), " ", nm)
    multiome <- getSilhouette(multiome, nredu=nm, dims=reduction_dims[[nm]], method="euclidean")

  }
}


if(metrics_vec["rangeclustering"]  == 1 ){
  message(date(), " preparing range clusters")
  for(nm in names(reduction_dims)[calc_dims]){
    message("        ", date(), " ", nm)
    multiome <- getRangeClusters(multiome, ASSAY= cf$ASSAY_1, reduction=nm,
                                 dims=reduction_dims[[nm]],
                                 #vres=seq(0.1, 3, by=0.1))
                                 vres=seq(2, 3, by=0.1))
  }
}
system(glue("cp {dataset}/save/multiome.Rds {dataset}/save/multiome.Rds.bak"))
save_object(multiome, file=file.path(dataset, "save/multiome.Rds"))
message("done.")
