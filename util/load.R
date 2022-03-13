library(SeuratData)

## need source(util.R)
loadmulti <- function(dataset, ASSAY_1, ASSAY_2, ASSAY_3=NULL, workdir="./"){
  if(dataset == "kidney"){
    mtx_path <- file.path(workdir, "data/MouseKidney/scATAC/matrix.mtx.gz")
    peak_path <- file.path(workdir, "data/MouseKidney/scATAC/peaks.bed.gz")
    barcode_path <- file.path(workdir, "data/MouseKidney/scATAC/barcodes.tsv.gz")
    features <-  data.table::fread(cmd = paste0("zcat < ", peak_path), sep="\t", header=FALSE) %>%
                        data.frame() %>% tidyr::unite(feature)
    barcodes <- data.table::fread(cmd = paste0("zcat < ", barcode_path), sep="\t", header=FALSE) %>%
                        data.frame() %>% tidyr::unite(barcode)



    peak_mtx <- Matrix::readMM(mtx_path) %>%
                         as("dgCMatrix") %>%
                         magrittr::set_rownames(features$feature) %>%
                         magrittr::set_colnames(barcodes$barcode)


    rna_mtx <- Read10X(file.path(workdir, "data/MouseKidney/scRNA/"))

    shared_cells <- intersect(colnames(peak_mtx), colnames(rna_mtx))



    meta <- read.csv(file.path(workdir, "data/MouseKidney/scRNA/GSM3271044_RNA_mouse_kidney_cell.txt"))
    rownames(meta) <- meta$sample
    meta <- meta[shared_cells,]
    RNA  <- rna_mtx[, shared_cells]
    ATAC <- peak_mtx[, shared_cells]
    peaks <- rownames(ATAC)
    ATAC <- ATAC[which(startsWith(peaks, "chr")), ]
    object <- CreateSeuratObject(counts=RNA, meta.data=meta, assay=ASSAY_1)
    object[[ASSAY_2]] <- CreateAssayObject(counts=ATAC, assay=ASSAY_2)
    object$celltype <- object$cell_name#####################################################
  } else if(dataset == "skin"){
   # peak-bc matrix
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
    object <- CreateSeuratObject(counts=RNA, meta.data=meta, assay=ASSAY_1)
    object[[ASSAY_2]] <- CreateAssayObject(counts=ATAC, assay=ASSAY_2)

  }else if(dataset == "skin_hairfollicle"){
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
    object <- CreateSeuratObject(counts=RNA, meta.data=meta, assay=ASSAY_1)
    object[[ASSAY_2]] <- CreateAssayObject(counts=ATAC, assay=ASSAY_2)
    Idents(object) <- "celltype"
    object <- subset(object, idents=c("TAC-1",
                                      "TAC-2",
                                      "Medulla",
                                      "IRS",
                                      "Hair Shaft-cuticle.cortex"))
  }else if(dataset == "skin_deMix"){
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
    object <- CreateSeuratObject(counts=RNA, meta.data=meta, assay=ASSAY_1)
    object[[ASSAY_2]] <- CreateAssayObject(counts=ATAC, assay=ASSAY_2)
    cells <- colnames(object)[!(object$celltype %in% "Mix")]
    Idents(object) <- "celltype"
    object <- subset(object, cells=cells)

  }else if(dataset == "pbmc"){
    mtxs <- Read10X_h5(file.path(workdir, "data/multiome_PBMC/output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"))
    RNA <- mtxs[["Gene Expression"]]
    ATAC <- mtxs[["Peaks"]]
    peaks <- rownames(ATAC)
    ATAC <- ATAC[which(startsWith(peaks, "chr")), ]
 
    meta <- read.csv(file.path(workdir, dataset, "save/annotation.tsv"), sep="\t")



    RNA <- RNA[, meta$barcode]
    ATAC <- ATAC[, meta$barcode]
    rownames(meta) <- meta$barcode
    object <- CreateSeuratObject(counts=RNA, meta.data=meta, assay=ASSAY_1)
    object[[ASSAY_2]] <- CreateAssayObject(counts=ATAC, assay=ASSAY_2)
    object$celltype <- object$annotation #####################################################
  }else if(dataset == "cite_30k"){
    object <- LoadData(ds = "bmcite")
    object$celltype = object$celltype.l2
    DefaultAssay(object) <- 'RNA'
    object <- NormalizeData(object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
    DefaultAssay(object) <- 'ADT'
    VariableFeatures(object) <- rownames(object[["ADT"]])
    object <- NormalizeData(object, normalization.method = 'CLR', margin = 2) %>%
                          ScaleData() %>% RunPCA(reduction.name = 'apca')
 }else if(dataset == "cite_elife"){
    object <- readRDS(file=file.path(workdir, "data/elife/swh:1:dir:c28aab565256405ee61ad2019c0746e83e3fe84e/5P-CITE-seq_Titration.rds"))
    object$celltype = as.character(object$fineCluster)
    object[["ADT"]] <- object[["ADT.kallisto"]]
    object[["RNA"]] <- object[["RNA.kallisto"]]
    DefaultAssay(object) <- 'RNA'
    object <- NormalizeData(object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
    DefaultAssay(object) <- 'ADT'
    VariableFeatures(object) <- rownames(object[["ADT"]])
    object <- NormalizeData(object, normalization.method = 'CLR', margin = 2) %>%
                          ScaleData() %>% RunPCA(reduction.name = 'apca')
    object$dataset <- "cite_elife"

  }else if(dataset == "ASAP_BM"){
    tmp <- load_object(file=file.path(workdir, "data/ASAP_BM/pbmc_LLL_processed.rds"))
    object <- CreateSeuratObject(counts=GetAssayData(tmp, slot="counts", assay="RNA"), assay="RNA", meta.data=tmp@meta.data)
    object$celltype <-object$predicted.celltype.l2
    object[["atac"]] <- CreateAssayObject(counts=GetAssayData(tmp, slot="counts", assay="peaks"))
    object[["adt"]] <- CreateAssayObject(counts=GetAssayData(tmp, slot="counts", assay="ADT"))
  }else{
    stop("wrong dataset name")
  }
  return(object)
}


