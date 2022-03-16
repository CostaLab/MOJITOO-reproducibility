library(Seurat)
library(dplyr)

source("../util/util.R")

workdir <- "../"
mex_dir_path <- file.path(workdir, "data/mouse_skin/ATAC/")

mtx_path <- paste(mex_dir_path, "GSM4156597_skin.late.anagen.counts.txt.gz", sep = '/')
feature_path <- paste(mex_dir_path, "GSM4156597_skin.late.anagen.peaks.bed.gz", sep = '/')
barcode_path <- paste(mex_dir_path, "GSM4156597_skin.late.anagen.barcodes.txt.gz", sep = '/')
features <-  data.table::fread(cmd = paste0("zcat < ", feature_path), sep="\t", header=FALSE) %>% 
                    data.frame() %>% tidyr::unite(feature)
barcodes <- data.table::fread(cmd = paste0("zcat < ", barcode_path), sep="\t", header=FALSE) %>% 
                    data.frame() %>% tidyr::unite(barcode)

peak_mtx <- Matrix::readMM(mtx_path) %>% 
                     as("dgCMatrix") %>% 
                     magrittr::set_rownames(features$feature) %>% 
                     magrittr::set_colnames(barcodes$barcode)
peaks <- rownames(peak_mtx)
peak_mtx<- peak_mtx[which(startsWith(peaks, "chr")), ]
rna_mtx <- Read10X(file.path(workdir, "data/mouse_skin/RNA/outs/filtered_feature_bc_matrix/"))
meta <- read.csv(file.path(workdir, "data/mouse_skin/ATAC/GSM4156597_skin_celltype.txt"), sep="\t")


RNA <- rna_mtx[, meta$rna.bc]
ATAC <- peak_mtx[, meta$rna.bc]
rownames(meta) <- meta$rna.bc

ASSAY_1 <- "RNA"
ASSAY_2 <- "ATAC"


object <- CreateSeuratObject(counts=RNA, meta.data=meta, assay=ASSAY_1)
object[[ASSAY_2]] <- CreateAssayObject(counts=ATAC, assay=ASSAY_2)


cell_num <- (1:10)*3000

set.seed(10)

new_order <- sample(colnames(object))
for(a_num in cell_num){
  cells <-new_order[1:a_num]
  a_sub <- subset(object, cells=cells)
  save_object(a_sub, glue("save/{a_num}.Rds"))
}

dir.create("save")
save_object(object, glue("save/all.Rds"))
