library(Seurat)
library(dplyr)
library(SeuratDisk)

source("../util/util.R")

nums <- (1:10)*3000

ASSAY_1 = "RNA"
ASSAY_2 = "ATAC"

for(i in nums){
  message(date(), " ", i)
  fname <- paste0("save/", i, ".Rds")
  object <- load_object(fname)

  RNA  <- GetAssayData(object, slot="counts", assay=ASSAY_1)
  ATAC <- GetAssayData(object, slot="counts", assay=ASSAY_2)
  meta <- object@meta.data
  
  rna <- CreateSeuratObject(counts=RNA, meta.data=meta, assay=ASSAY_1)
  SaveH5Seurat(rna, filename = glue("save/{i}RNA.h5Seurat"))
  Convert(glue("save/{i}RNA.h5Seurat"), dest = "h5ad")
  
  atac <- CreateSeuratObject(counts=ATAC, meta.data=meta, assay=ASSAY_2)
  SaveH5Seurat(atac, filename = glue("save/{i}ATAC.h5Seurat"))
  Convert(glue("save/{i}ATAC.h5Seurat"), dest = "h5ad")
}

