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
mute(library(dplyr))
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
mute(library(scAI))
mute(library(MOFA2))
mute(library(rbenchmark))
mute(library(rliger))
mute(library(peakRAM))

set.seed(10)


AllOptions <- function(){
  parser <- OptionParser()
  parser <- add_option(parser, c("-t", "--tool"), type="character", default="WNN",
                       help="tool default[default %default]",
                       metavar="character")

  return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)


tool <- pa$tool



source("../util/util.R")


ASSAY_1 = "RNA"
ASSAY_2 = "ATAC"



preprocess <- function(object){

    DefaultAssay(object) <- ASSAY_1
    object <-  NormalizeData(object)
    object  <- FindVariableFeatures(object, nfeatures=3000, verbose=F)
    object  <- ScaleData(object, nfeatures=3000, verbose=F)
    object <- RunPCA(object, npcs=50, reduction.name="RNA_PCA", verbose=F)

    DefaultAssay(object) <- ASSAY_2

    object <- RunTFIDF(object, verbose=F)
    object <- FindTopFeatures(object, min.cutoff = 'q0', verbose=F)
    object <- RunSVD(object, verbose=F)
    return(object)
}



symphony_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="ATAC",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      reduction.name = "symphony"
){


  library(mixOmics)
  DefaultAssay(object) <- ASSAY_1
  object <-  NormalizeData(object)
  object  <- FindVariableFeatures(object, nfeatures=3000, verbose=F)
  object  <- ScaleData(object, nfeatures=3000, verbose=F)

  DefaultAssay(object) <- ASSAY_2

  object <- RunTFIDF(object, verbose=F)
  object <- FindVariableFeatures(object, nfeatures=5000)
  object <- ScaleData(object, features=VariableFeatures(object))

  DefaultAssay(object) <- ASSAY_1

  Y = 1:ncol(object) ## set all Y to be
  X1 = GetAssayData(object, assay=ASSAY_1, slot="scale.data") ## use normalized data
  X2 = GetAssayData(object, assay=ASSAY_2, slot="scale.data") ## use normalized data

  cc_out <- cc(t(X1), t(X2))
  dimension <- cc_out$scores$xscores[, 1:20]

  dimension <- cbind(diablo$variates[[ASSAY_1]], diablo$variates[[ASSAY_2]])
  dims <- ncol(dimension)
  object[[reduction.name]] <- CreateDimReducObject(dimension, assay = "RNA", key=reduction.name)
  object <- RunUMAP(object, reduction=reduction.name, dims=1:dims, reduction.name =glue("{reduction.name}_UMAP"), verbose=F)

  return(object)
}
symphonyRed_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="ATAC",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      reduction.name = "symphony"
){


  library(mixOmics)
  object <- preprocess(object)

  DefaultAssay(object) <- ASSAY_1

  X1 = Embeddings(object[[nredu_1]])[, dims_1] ## use normalized data
  X2 = Embeddings(object[[nredu_2]])[, dims_2] ## use normalized data

  cc_out <- cc((X1), (X2))
  dimension <- cc_out$scores$xscores[, 1:20]

  dims <- ncol(dimension)
  object[[reduction.name]] <- CreateDimReducObject(dimension, assay = "RNA", key=reduction.name)
  object <- RunUMAP(object, reduction=reduction.name, dims=1:dims, reduction.name =glue("{reduction.name}_UMAP"), verbose=F)

  return(object)
}

MOFA_run <- function(object,
                     ASSAY_1="RNA",
                     ASSAY_2="ATAC"
){


  object <- preprocess(object)
  DefaultAssay(object) <- ASSAY_2
  object <- FindTopFeatures(object, min.cutoff = 2000, verbose=F)
  object <- ScaleData(object, verbose=F)
  mofa <- create_mofa(object, assays = c(ASSAY_1,ASSAY_2))
  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- 15
  mofa <- prepare_mofa(mofa,
    model_options = model_opts
  )
  mofa <- run_mofa(mofa)
}

schema_run <- function(dataset=3000
){
  system(glue("python 00_schema_prepare.py -d {dataset}"))
}

RAW_run <- function(object){
  preprocess(object)
}


MOJITOO_run <- function(object,
                        ASSAY_1="RNA",
                        ASSAY_2="ATAC",
                        nredu_1="RNA_PCA",
                        nredu_2="lsi",
                        dims_1 = 1:50,
                        dims_2 = 2:50,
                        reduction.name = "MOJITOO"

){
  object <- preprocess(object)
  mojitoo.Seurat(object,
         reduction.list = list(nredu_1, nredu_2),
         dims.list = list(dims_1, dims_2),
         reduction.name="MOJITOO"
  )

}

MOFARed_run <- function(object,
                    ASSAY_1="RNA",
                    ASSAY_2="ATAC",
                    nredu_1="RNA_PCA",
                    nredu_2="lsi",
                    dims_1 = 1:50,
                    dims_2 = 1:50,
                    name = "skin"
){

  mute(library(MOFA2))
  object <- preprocess(object)
  mofa <- create_mofa_from_matrix(
    list(t(Embeddings(object[[nredu_1]])[, dims_1]),
         t(Embeddings(object[[nredu_2]])[, dims_2])
    )
  )
  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- 15
  mofa <- prepare_mofa(mofa,
    model_options = model_opts
  )
  mofa <- run_mofa(mofa)

  plot_factor_cor(mofa)
  factors <- 1:get_dimensions(mofa)[["K"]]
  Z <- get_factors(mofa, factors = factors, groups = "all")[[1]]
  Z <- Z[colnames(object), ]
  object[["MOFARed"]] <- CreateDimReducObject(embeddings=as.matrix(Z), key="factor", assay=ASSAY_1)

  return(object)

}

DIABLORed_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="ATAC",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      reduction.name = "DIABLORed"
){


  library(mixOmics)
  object <- preprocess(object)
  Y = 1:ncol(object) ## set all Y to be
  X1 = Embeddings(object[[nredu_1]])[, dims_1] ## use normalized data
  X2 = Embeddings(object[[nredu_2]])[, dims_2] ## use normalized data

  X <- list(X1, X2)
  names(X) <- c(ASSAY_1, ASSAY_2)
  design =  matrix(c(0,1,1,0),ncol=2, nrow=2, byrow=TRUE,
                   dimnames = list(names(X), names(X)))
  diablo <- block.plsda(X,Y, design=design, ncomp=30)

  dimension <- cbind(diablo$variates[[ASSAY_1]], diablo$variates[[ASSAY_2]])
  dims <- ncol(dimension)
  object[[reduction.name]] <- CreateDimReducObject(dimension, assay = "RNA", key=reduction.name)
  object
}


WNN_run <- function(object){
  object <- preprocess(object)
  object <- FindMultiModalNeighbors(object,
                                    reduction.list = list("RNA_PCA", "lsi"),
                                    dims.list = list(1:50, 2:50))
}

scAI_run <- function(object){
  start_time <- Sys.time()
  object <- preprocess(object)
  DefaultAssay(object) <- ASSAY_2
  object <- FindVariableFeatures(object, nfeatures=5000)
  object <- ScaleData(object, features=VariableFeatures(object))

  mtx_list <- list(rna = (GetAssayData(object, slot="counts", assay=ASSAY_1)[object@assays[[ASSAY_1]]@var.features, ]),
                   atac = (GetAssayData(object, slot="counts", assay=ASSAY_2)[object@assays[[ASSAY_2]]@var.features, ]))

  X <- list(RNA=mtx_list[[1]], ATAC=mtx_list[[2]])
  labels <- data.frame(labels=object$celltype) # the prior labels of cells
  scAI_outs <- create_scAIobject(raw.data = X, do.sparse = F)
  scAI_outs <- preprocessing(scAI_outs, assay = NULL)
  scAI_outs <- addpData(scAI_outs, pdata = labels, pdata.name = "Cell types")

  scAI_outs <- run_scAI(scAI_outs, K = 20, nrun = 1, do.fast = T)
  end_time <- Sys.time()
  print(end_time - start_time)
}

liger_run <- function(object){

  start_time <- Sys.time()

  liger <- rliger::createLiger(list(peaks = unshared_atac))
  liger <- rliger::normalize(liger)
  norm <- liger@norm.data$peaks

  se = CreateSeuratObject(norm)
  vars_2000 <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
  top2000 <- head(VariableFeatures(vars_2000),2000)
  top2000_feats <-  norm[top2000,]
  liger <- selectGenes(liger)
  liger@var.genes <- top2000
  liger <- scaleNotCenter(liger)
  unshared_feats = liger@scale.data$peaks


  liger <- rliger::createLiger(list(rna = rna, atac = shared_atac))
  liger <- rliger::normalize(liger)
  liger <- selectGenes(liger, var.thresh = 0.1, datasets.use =1 , unshared = TRUE,  unshared.datasets = list(2), unshared.thresh= 0.2)
  liger <- scaleNotCenter(liger)
  peak_names <- rownames(unshared_feats)
  liger@var.unshared.features[[2]] = peak_names
  liger@scale.unshared.data[[2]] = t(unshared_feats)
  liger <- optimizeALS(liger, k=30, use.unshared = TRUE, max_iters =30,thresh=1e-10)
  liger <- quantile_norm(liger)
  liger <- louvainCluster(liger)
  liger <- runUMAP(liger) ## by default on H.norm
  umap_plots <-plotByDatasetAndCluster(liger, axis.labels = c("UMAP1","UMAP2"), return.plots = TRUE)
  umap_plots[[2]]

  reduction <-  liger@H.norm[, 1:ncol(liger@H.norm)]
  rna_idx <-  which(startsWith(rownames(reduction), "rna_") )
  atac_idx <-  which(!startsWith(rownames(reduction), "rna_") )
  reduction_rna <- reduction[rna_idx, ]
  reduction_atac <- reduction[atac_idx, ]
  rownames(reduction_rna) <-  stringr::str_replace(rownames(reduction_rna),"rna_", "")
  reduction_rna <- reduction_rna[rownames(reduction_atac), ]

  bind_reduction <- cbind(reduction_rna, reduction_atac)

  cells = rownames(bind_reduction)
  nd_time <- Sys.time()
  print(end_time - start_time)
}


nums <- (1:10)*3000

run_tool <- function(tool){
  func_name = paste0(tool, "_run")
  func_call = paste0(func_name, "(object)")
  message(paste("executing", func_name))
  eval(parse(text=func_call))
}

df_list <- list()
object <- NULL

df_list <- list()


rna <<- NULL
shared_atac <<- NULL
unshared_atac <<- NULL
shared_atac_all <- load_object(file=glue("../skin/save/shared_mtx.Rds"))
unshared_atac_all <- load_object(file=glue("../skin/save/unshared_mtx.Rds"))

fn <- paste0("save/", "all.Rds")
object_all <- load_object(fn)


shared_atac_all <- shared_atac_all[rownames(shared_atac_all) != "", ]
atac_all_bc <- stringr::str_replace_all(colnames(shared_atac_all), ",", ".")
map_vec <- unname(object_all$rna.bc)
names(map_vec) <- unname(object_all$atac.bc)
rna_bc_all <-  plyr::revalue(atac_all_bc, map_vec, warn_missing=F)
colnames(shared_atac_all) <- rna_bc_all
rna_all <- GetAssayData(object_all, assay="RNA", slot="counts")
shared_bc_all <- Reduce(intersect, list(colnames(rna_all),colnames(shared_atac_all),colnames(unshared_atac_all)))
sample_shared_bc_all <- sample(shared_bc_all)


for(num in nums){
  fn <- paste0("save/", num, ".Rds")
  num <- as.character(num)
  if(tool != "liger"){
    object <- load_object(fn)
  }
  if(tool == "liger"){
    dataset <- "skin"
    shared_bc <- sample_shared_bc_all[1:min(length(shared_bc_all), as.integer(num))]
    message("actual num: ", length(shared_bc))
    rna <- rna_all[, shared_bc]
    shared_atac <- shared_atac_all[, shared_bc]
    unshared_atac  <- unshared_atac_all[, shared_bc]
    colnames(rna) <- paste0("rna_", colnames(rna))
  }


  message(date(), "  ====================num=================:  ", num)
  df <- peakRAM(run_tool(tool))
  df$tool = tool
  df$num = num
  if(tool == "liger"){
    df$num <- as.character(length(shared_bc))
  }
  df_list[[num]] <- df
}

dff <- do.call(rbind, df_list)
saveRDS(dff, glue("save/{tool}_RAM_TIME.RDS"))
saveRDS(df_list, glue("save/{tool}_RAM_TIME_list.RDS"))
