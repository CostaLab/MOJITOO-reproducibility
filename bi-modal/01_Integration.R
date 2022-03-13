mute <- suppressPackageStartupMessages
mute(library(Seurat))
mute(library(Signac))
mute(library(MOJITOO))
mute(library(dplyr))
mute(library(mclust))
mute(library(data.table))
mute(library(bench))
mute(library(optparse))
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html

source("../util/util.R")
source("../util/load.R")
source("../util/peaks_tile.R")


AllOptions <- function(){
  parser <- OptionParser()
  parser <- add_option(parser, c("-d", "--dir"), type="character", default="kidney",
                       help="directory [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-i", "--init"), type="character", default="FALSE",
                       help="directory [default %default]",
                       metavar="character")

  return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)


dataset <- pa$dir
INIT <- pa$init

cf <- getSettings(dataset = dataset)

if(INIT == "TRUE"){
  multiome <- loadmulti(
      dataset = dataset,
      ASSAY_1=cf$ASSAY_1,
      ASSAY_2=cf$ASSAY_2,
      workdir="./"
  )
  bench_list <- list()
}else if(INIT == "FALSE"){
  multiome <- load_object(file.path(dataset, "save/multiome.Rds"))
  bench_list <- load_object(file.path(dataset, "save/bench_list.Rds"))
}
#########################################################################

UPDATE_VEC <- c(
    "Signac"             =   1,
    "MOJITOO"            =   1,
    "MOFA"               =   1,
    "MOFARed"            =   1,
    "DIABLORed"          =   1,
    "symphonyRed"        =   1,
    "WNN"                =   1,
    "schema"             =   1,
    "liger"              =   1,
    "scAI"               =   1
)

if(INIT == "TRUE" ){
  UPDATE_VEC["Signac"] = 1
}

message("dataset: ", dataset)
message("run: ", paste(names(UPDATE_VEC[UPDATE_VEC == 1]), sep=" "))

if(UPDATE_VEC["Signac"] == 1){
  message(">>>>>>>>>>>>>>>>>>>>>>>Signac ", date())
  multiome <- SignacProcess(object=multiome,
                           ASSAY_1=cf$ASSAY_1,
                           ASSAY_2=cf$ASSAY_2,
                           nredu_1=cf$Redu_1,
                           nredu_2=cf$Redu_2,
                           dims_1=cf$dims_1,
                           dims_2=cf$dims_2,
                           name = dataset
  )
}



if(UPDATE_VEC["liger"] == 1){
  message(">>>>>>>>>>>>>>>>>>>>>>>liger ", date())
  a_bench <- bench::mark(
  multiome <- liger_run(object=multiome,
                 ASSAY_1=cf$ASSAY_1,
                 ASSAY_2=cf$ASSAY_2,
                 nredu_1=cf$Redu_1,
                 nredu_2=cf$Redu_2,
                 dims_1=cf$dims_1,
                 dims_2=cf$dims_2)
  )
  a_bench$result <- NULL
  bench_list[["liger"]] <- a_bench
}

if(UPDATE_VEC["MOJITOO"] == 1){
message(">>>>>>>>>>>>>>>>>>>>>>>MOJITOO ", date())
  a_bench <- bench::mark(
  multiome <- MOJITOO_run(
       object=multiome,
       ASSAY_1=cf$ASSAY_1,
       ASSAY_2=cf$ASSAY_2,
       nredu_1=cf$Redu_1,
       nredu_2=cf$Redu_2,
       dims_1=cf$dims_1,
       dims_2=cf$dims_2
  )
  )
  a_bench$result <- NULL
  bench_list[["MOJITOO"]] <- a_bench
}



if(UPDATE_VEC["DIABLORed"] == 1){
message(">>>>>>>>>>>>>>>>>>>>>>>DIABLORed ", date())
  a_bench <- bench::mark(
  multiome <- DIABLORed_run(
       object=multiome,
       ASSAY_1=cf$ASSAY_1,
       ASSAY_2=cf$ASSAY_2,
       nredu_1=cf$Redu_1,
       nredu_2=cf$Redu_2,
       dims_1=cf$dims_1,
       dims_2=cf$dims_2
  )
  )
  a_bench$result <- NULL
  bench_list[["DIABLORed"]] <- a_bench
}

if(UPDATE_VEC["symphonyRed"] == 1){
message(">>>>>>>>>>>>>>>>>>>>>>>symphonyRed ", date())
  a_bench <- bench::mark(
  multiome <- symphonyRed_run(
       object=multiome,
       ASSAY_1=cf$ASSAY_1,
       ASSAY_2=cf$ASSAY_2,
       nredu_1=cf$Redu_1,
       nredu_2=cf$Redu_2,
       dims_1=cf$dims_1,
       dims_2=cf$dims_2,
       dataset=dataset
  )
  )
  a_bench$result <- NULL
  bench_list[["symphonyRed"]] <- a_bench
}


if(UPDATE_VEC["MOFARed"] == 1){
message(">>>>>>>>>>>>>>>>>>>>>>>MOFARed ", date())
  a_bench <- bench::mark(
  multiome <- MOFARed_run(
       object=multiome,
       ASSAY_1=cf$ASSAY_1,
       ASSAY_2=cf$ASSAY_2,
       nredu_1=cf$Redu_1,
       nredu_2=cf$Redu_2,
       dims_1=cf$dims_1,
       dims_2=cf$dims_2
  )
  )
  a_bench$result <- NULL
  bench_list[["MOFARed"]] <- a_bench
}

if(UPDATE_VEC["schema"] == 1){
schema_dims=1:50
if(startsWith(dataset, "cite")){
  schema_dims=1:10
}
message(">>>>>>>>>>>>>>>>>>>>>>>schema ", date())
  a_bench <- bench::mark(
  multiome <- schema_run(
       object=multiome,
       ASSAY_1=cf$ASSAY_1,
       ASSAY_2=cf$ASSAY_2,
       dims = schema_dims
  )
  )
  a_bench$result <- NULL
  bench_list[["schema"]] <- a_bench
}

if(UPDATE_VEC["MOFA"] == 1){
  message(">>>>>>>>>>>>>>>>>>>>>>>MOFA ", date())
  a_bench <- bench::mark(
  multiome <- MOFA_run(object=multiome,
                 ASSAY_1=cf$ASSAY_1,
                 ASSAY_2=cf$ASSAY_2)
  )
  a_bench$result <- NULL
  bench_list[["MOFA"]] <- a_bench
}


if(UPDATE_VEC["scAI"] == 1){ ## don't benchmark ate too much memory
  message(">>>>>>>>>>>>>>>>>>>>>>>scAI ", date())
  multiome <- scAI_run(object=multiome,
                 ASSAY_1=cf$ASSAY_1,
                 ASSAY_2=cf$ASSAY_2,
                 nredu_1=cf$Redu_1,
                 nredu_2=cf$Redu_2,
                 dims_1=cf$dims_1,
                 dims_2=cf$dims_2,
                 rerun=FALSE,
                 #rerun=TRUE,
                 dataset=dataset)
}

if(UPDATE_VEC["WNN"] == 1){
  message(">>>>>>>>>>>>>>>>>>>>>>>WNN ", date())
  a_bench <- bench::mark(
  multiome <- FindMultiModalNeighbors(multiome, reduction.list = list(cf$Redu_1, cf$Redu_2), dims.list = list(cf$dims_1, cf$dims_2)) %>% RunUMAP(nn.name = "weighted.nn", reduction.name = "WNN_UMAP", reduction.key = "wnnUMAP_", verbose=F)
  )
  a_bench$result <- NULL
  bench_list[["WNN"]] <- a_bench
}

message("end ", date())
dir.create(glue("{dataset}/save"))
multiome@tools$dataset=dataset
if(class(multiome) == "Seurat"){
  system(glue("cp {dataset}/save/multiome.Rds {dataset}/save/multiome.Rds.bak"))
  save_object(multiome, file=file.path(dataset, "save/multiome.Rds" ))
  save_object(bench_list, file=file.path(dataset,"save/bench_list.Rds"))
}else{
  stop("wrong object class")
}

