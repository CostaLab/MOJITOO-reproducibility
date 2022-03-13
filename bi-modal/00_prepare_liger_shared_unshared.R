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
  parser <- add_option(parser, c("-i", "--INIT"), type="character", default="unshared",
                       help="shared or unshared [default %default]",
                       metavar="character")

  return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)

dataset <- pa$dir
which_mtx <- pa$INIT

is_sort = F
if(dataset == "pbmc"){
  organ = "hg38"
  fragment_file = "data/multiome_PBMC/output/pbmc_granulocyte_sorted_10k_atac_fragments_ordered.tsv"

}
if(dataset == "skin"){
  organ = "mm10"
  fragment_file = "data/mouse_skin/raw/GSM4156597_skin.late.anagen.atac.fragments.bed.gz"
  issort = T
}

if(dataset == "TEA"){
  organ = "hg38"
  fragment_file = "data/TEA/data/merged/fragments_merge.tsv"
  issort = T
}



if(which_mtx == "shared"){
  shared_mtx <- shared_counts(fragment_file, dataset=dataset, organ=organ, sort_fragment=issort)
  save_object(shared_mtx, glue("{dataset}/save/shared_mtx.Rds"))
}

if(which_mtx == "unshared"){
  object <- load_object(glue("{dataset}/save/multiome.Rds"))
  if(dataset %in% c("TEA", "ASAP_BM")){
    peaks <- GetAssayData(object, assay="atac", slot="counts")
  }else{
    peaks <- GetAssayData(object, assay="Peaks", slot="counts")
  }
  rownames(peaks) <- stringr::str_replace(rownames(peaks), "_", "-")
  rownames(peaks) <- stringr::str_replace(rownames(peaks), ":", "-")
  unshared_mtx <- unshared_bin_mtx(peaks, organ=organ, cpus=30)
  save_object(unshared_mtx, glue("{dataset}/save/unshared_mtx.Rds"))
}



