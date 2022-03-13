mute <- suppressPackageStartupMessages
mute(library(Seurat))
mute(library(rliger))
mute(library(mclust))
mute(library(glue))
mute(library(FNN))
mute(library(RANN)) ## nn2
mute(library(dplyr))
mute(library(reshape2))
mute(library(cli))
mute(library(scAI))
mute(library(Rcpp))
mute(library(RcppArmadillo))
mute(library(plyr))
mute(library(cluster))
mute(library(parallelDist))
mute(library(S4Vectors))
mute(library(Matrix))
mute(library(assertthat))
mute(library(CCA))
mute(library(SeuratDisk))
mute(library(knn.covertree))
mute(library(reticulate))

cite_elife_dic <- c("0" = "T.CD4",
                    "3" = "T.CD.EM",
                    "14" = "T.Treg",
                    "2" = "T.CD8",
                    "7" = "T.CD8.CM",
                    "10" = "T.CD8.EM",
                    "16" = "T.Prolif",
                    "9" = "NK",
                    "6" = "B.Blood",
                    "13" = "B.Lung",
                    "8" = "B.Plasma.1",
                    "20" = "B.Plasma.2",
                    "11" = "Mono.Int",
                    "15" = "mDC.Blood",
                    "17" = "pDC",
                    "12" = "Macro",
                    "1" = "Mono.1",
                    "4" = "Mono.2",
                    "5" = "Mono.3",
                    "19" = "mDC.Lung",
                    "18" = "Fibro",
                    "21" = "Epithelia"
)




getSettings <- function(dataset){
    SL <- SimpleList(
      ASSAY_1="RNA",
      ASSAY_2="Peaks",
      Redu_1 = "RNA_PCA",
      Redu_2 = "lsi",
      dims_1 = 1:50,
      dims_2 = 2:50
    )

  if(dataset=="cite_30k"){
    SL <- SimpleList(
      Redu_1 = "pca",
      Redu_2 = "apca",
      dims_1 = 1:50,
      dims_2 = 1:24,
      ASSAY_1 = "RNA",
      ASSAY_2 = "ADT"
    )
  } else if(dataset=="cite_elife"){
    SL <- SimpleList(
      Redu_1 = "pca",
      Redu_2 = "apca",
      dims_1 = 1:50,
      dims_2 = 1:30,
      ASSAY_1 = "RNA",
      ASSAY_2 = "ADT"
    )
  }

  SL$dims_add = 1:(length(SL$dims_1) + length(SL$dims_2))
  SL$dims_min = 1:(min(length(SL$dims_1) , length(SL$dims_2)))
  SL$dims_bind = 1:(2*min(length(SL$dims_1) , length(SL$dims_2)))

  return(SL)
}

#Thu Nov 25 22:20:04 2021 data pbmc
#pbmc: pval>0.05 46
#Thu Nov 25 22:20:42 2021 data skin
#skin: pval>0.05 45
#Thu Nov 25 22:21:59 2021 data kidney

save_object <- function(object, file, file_format="lz4"){

  stopifnot(file_format %in% c("zstd", "lz4", "gzip", "bzip2", "xz", "nocomp"))

  if(file_format %in% "nocomp"){
    saveRDS(object = object, file = file, compress = FALSE)
    return(invisible(NULL))
  }

  if(file_format %in% c("zstd", "lz4")){
    con <- archive::file_write(file = file, filter = file_format)
    open(con)
    saveRDS(object = object, file = con)
    close(con)
  }else{
    saveRDS(object = object, file = file, compress = file_format)
  }
}

load_object <- function(file){
  con <- archive::file_read(file = file)
  res <- readRDS(file = con)
  close(con)
  return(res)
}



## calculate the distance matrix to upper
getDistanceVector <- function(object,
                              nredu="lsi",
                              dims=1:50,
                              method = "euclidean"
){
  to.upper<-function(X) Matrix::t(X)[lower.tri(X,diag=FALSE)]
  if(nredu =="WNN"){
    d <- 1- object@graphs$wsnn
    vec <- to.upper(d)

  }else{
    m <- Embeddings(object, reduction=nredu)[, dims]
    d <- parDist(m, method = method)
    vec <- as.vector(d)
  }
  rm(d)
  gc()
  return(vec)
}

getStructure <- function(object,
                         nredu_1 = "RNA_PCA",
                         nredu_2 = "lsi",
                         dims_1 = 1:50,
                         dims_2 = 1:50,
                         nredu="lsi",
                         dims=1:50,
                         method = "euclidean",
                         cor_method = "pearson"
){
  objectbak <- object
  if(nredu== "liger"){
    object <- objectbak@tools$liger
  }

  if (!("DEmbed" %in% names(object@tools))){
    object@tools$DEmbed <- list()
  }
  if(!(nredu_1  %in%names(object@tools$DEmbed))){
    object@tools$DEmbed[[nredu_1]] <- getDistanceVector(object, nredu_1, dims_1)
  }
  if(!(nredu_2  %in%names(object@tools$DEmbed))){
    object@tools$DEmbed[[nredu_2]] <- getDistanceVector(object, nredu_2, dims_2)
  }
  for(nm in names(object@tools$DEmbed)){
      if ((nm!=nredu_1) & (nm!=nredu_2)){
        object@tools$DEmbed[[nm]] <- NULL
      }
  }
  vec <- getDistanceVector(object, nredu, dims, method=method)

  if(cor_method == "pearson"){
    Correlation <- c(cor(object@tools$DEmbed[[nredu_1]], vec),
                   cor(object@tools$DEmbed[[nredu_2]], vec)
                  )

    object@tools[[glue("structure_{method}")]][[nredu]] <- Correlation
    objectbak@tools[[glue("structure_{method}")]][[nredu]] <- Correlation
  }else if(cor_method == "spearman"){
    Correlation <- c(cor(object@tools$DEmbed[[nredu_1]], vec, method="spearman"),
                   cor(object@tools$DEmbed[[nredu_2]], vec, method="spearman")
                  )
    object@tools[[glue("structure_{method}_spearman")]][[nredu]] <- Correlation
    objectbak@tools[[glue("structure_{method}_spearman")]][[nredu]] <- Correlation
  }

  if(nredu == "liger"){
    objectbak@tools$liger <- object
  }else{
    objectbak <- object
  }
  rm(object)
  gc()

  return(objectbak)
}

getNN <- function(object,
                  nredu="lsi",
                  dims=1:50,
                  k=30
){
  if(nredu=="WNN"){
    wnn <- object@graphs$wsnn
    nn_ <- matrix(rep(-1, ncol(object)*k), nrow=ncol(object), ncol=k)
    cli_progress_bar("WNN", total = ncol(object))
    for(i in 1:ncol(object)){
      indexs <- which(wnn[i, ] > 0)
      dff = data.frame(index=indexs, value=wnn[i, indexs]) %>%
                    filter(index !=i) %>%
                    top_n(n=30, wt=value)%>%
                    dplyr::arrange(-value)

      indexs = dff$index[1:min(30, nrow(dff))]
      nn_[i, 1:min(30, length(indexs))] <- indexs
      cli_progress_update()
    }
  }else{
    redu <- Embeddings(object[[nredu]])[, dims]
    nn_ <- get.knn(redu, k=30)$nn.index
  }
  if (!("NN" %in% names(object@tools))){
    object@tools$NN <- list()
  }
  object@tools$NN[[nredu]] <- nn_
  return(object)
}


SignacProcess <- function(object,
                          ASSAY_1="RNA",
                          ASSAY_2="Peaks",
                          nredu_1="RNA_PCA",
                          nredu_2="lsi",
                          dims_1 = 1:50,
                          dims_2 = 1:50,
                          name="yyy"
){
  if((name == "cite_30k") | name=="cite_elife"){
    message(">>>>>>>>cite norm<<<<<<<<<<<")
    DefaultAssay(object) <- 'ADT'
    VariableFeatures(object) <- rownames(object[["ADT"]])
    object <- NormalizeData(object, normalization.method = 'CLR', margin = 2) %>%
                          ScaleData() %>% RunPCA(reduction.name = 'apca')

    #do nothing
  }else{
    DefaultAssay(object) <- ASSAY_1
    #added just @20211123
    object <-  NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
    object <- FindVariableFeatures(object, nfeatures=3000)
    object <- ScaleData(object, nfeatures=3000)
    object <-  RunPCA(object, npcs=50, reduction.name="RNA_PCA")


    DefaultAssay(object) <- ASSAY_2

    object <- RunTFIDF(object)
    object <- FindTopFeatures(object, min.cutoff = 'q0')
    object <- RunSVD(object)

    object <- RunUMAP(object, reduction.name=glue("{nredu_1}_UMAP"), reduction=glue("{nredu_1}"), dims=dims_1, verbose=F)
    object <- RunUMAP(object, reduction.name=glue("{nredu_2}_UMAP"), reduction=glue("{nredu_2}"), dims=dims_2, verbose=F)
  }

  DefaultAssay(object) <- ASSAY_2
  object <- FindVariableFeatures(object, nfeatures=5000)
  object <- ScaleData(object, features=VariableFeatures(object))

  return(object)
}


MOFA_run <- function(object,
                    ASSAY_1,
                    ASSAY_2,
                    name = "skin"
){

  mute(library(MOFA2))

  mofa <- create_mofa(object, assays = c(ASSAY_1,ASSAY_2))
  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- 15
  mofa <- prepare_mofa(mofa,
    model_options = model_opts
  )
  mofa <- run_mofa(mofa)
  #saveRDS(mofa, file="mofa.Rds")

  plot_factor_cor(mofa)
  factors <- 1:get_dimensions(mofa)[["K"]]

  mofa <- run_umap(mofa,
    factors = factors,
    n_neighbors = 15,
    min_dist = 0.30
  )

  mofaUMAP <- mofa@dim_red$UMAP
  rownames(mofaUMAP) <- paste0(mofaUMAP$sample)
  assertthat::assert_that(all(rownames(mofaUMAP  %in% colnames(object))))
  assertthat::assert_that(all(colnames(object) %in% rownames(mofaUMAP)))
  mofaUMAP$sample = NULL
  colnames(mofaUMAP) <- paste0("UMAP_", 1:2)
  mofaUMAP <- mofaUMAP[colnames(object), ]

  object[["MOFA_UMAP"]] <- CreateDimReducObject(embeddings=as.matrix(mofaUMAP), key="mofa", assay=ASSAY_1)

  #factors <- 1:get_dimensions(mofa)[["K"]]
  Z <- get_factors(mofa, factors = factors, groups = "all")[[1]]
  Z <- Z[colnames(object), ]
  object[["MOFA"]] <- CreateDimReducObject(embeddings=as.matrix(Z), key="factor", assay=ASSAY_1)

  return(object)

}

MOFARed_run <- function(object,
                    ASSAY_1,
                    ASSAY_2,
                    nredu_1="RNA_PCA",
                    nredu_2="lsi",
                    dims_1 = 1:50,
                    dims_2 = 1:50,
                    name = "skin"
){

  mute(library(MOFA2))

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

  mofa <- run_umap(mofa,
    factors = factors,
    n_neighbors = 15,
    min_dist = 0.30
  )

  mofaUMAP <- mofa@dim_red$UMAP
  rownames(mofaUMAP) <- paste0(mofaUMAP$sample)
  assertthat::assert_that(all(rownames(mofaUMAP  %in% colnames(object))))
  assertthat::assert_that(all(colnames(object) %in% rownames(mofaUMAP)))
  mofaUMAP$sample = NULL
  colnames(mofaUMAP) <- paste0("UMAPRed_", 1:2)
  mofaUMAP <- mofaUMAP[colnames(object), ]

  object[["MOFARed_UMAP"]] <- CreateDimReducObject(embeddings=as.matrix(mofaUMAP), key="mofaRed", assay=ASSAY_1)

  Z <- get_factors(mofa, factors = factors, groups = "all")[[1]]
  Z <- Z[colnames(object), ]
  object[["MOFARed"]] <- CreateDimReducObject(embeddings=as.matrix(Z), key="factor", assay=ASSAY_1)

  return(object)

}

schema_run <- function(object,
                       ASSAY_1,
                       ASSAY_2,
                       dims=1:50,
                       name="skin"
){

  embedd <- read.csv(glue("{dataset}/save/schema.csv"), row.names = 1)
  object[["schema"]] <- CreateDimReducObject(embeddings=as.matrix(embedd[colnames(object), ]), key="schema", assay=ASSAY_1)
  object <- RunUMAP(object, reduction="schema", reduction.name=glue("schema_UMAP"), dims=dims, verbose=F)
  return(object)
}




liger_cite <- function(object,
                        ASSAY_1="RNA",
                        ASSAY_2="Peaks",
                        nredu_1="RNA_PCA",
                        nredu_2="apca",
                        dims_1 = 1:50,
                        dims_2 = 1:50,
                        reduction.name = "liger"
){
  #http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/UINMF_vignette.html

  dataset <- object@tools$dataset
  rna <-GetAssayData(object, slot="counts", assay=ASSAY_1)
  adt <-GetAssayData(object, slot="counts", assay=ASSAY_2)
  colnames(rna) <- paste0("rna_", colnames(rna))

  liger <- rliger::createLiger(list(rna = rna, adt = adt))
  liger <- rliger::normalize(liger)
  liger <- selectGenes(liger, unshared = TRUE, unshared.datasets = list(1), unshared.thresh= 0.4)
  liger <- scaleNotCenter(liger)
  if(dataset=="cite_30k"){
    liger <- optimizeALS(liger, k=11, use.unshared = TRUE) ## maximum 11, otherwise error
  }else if(dataset=="cite_elife"){
    liger <- optimizeALS(liger, k=20, use.unshared = TRUE) ## maximum 11, otherwise error
  }
  liger <- quantile_norm(liger)
  liger <- louvainCluster(liger)

  reduction <-  liger@H.norm[, 1:ncol(liger@H.norm)]
  rna_idx <-  which(startsWith(rownames(reduction), "rna_") )
  adt_idx <-  which(!startsWith(rownames(reduction), "rna_") )
  reduction_rna <- reduction[rna_idx, ]
  reduction_adt <- reduction[adt_idx, ]
  rownames(reduction_rna) <-  stringr::str_replace(rownames(reduction_rna),"rna_", "")
  reduction_rna <- reduction_rna[rownames(reduction_adt), ]

  bind_reduction <- cbind(reduction_rna, reduction_adt)

  cells = rownames(bind_reduction)
  subobject <- CreateSeuratObject(counts=GetAssayData(object, assay=ASSAY_1, slot="counts")[1:20, cells], meta.data=object@meta.data[cells, ])
  subobject[[ASSAY_2]] <- CreateAssayObject(counts=GetAssayData(object, assay=ASSAY_2, slot="counts")[1:20, cells])
  subobject[[nredu_1]] <- CreateDimReducObject(Embeddings(object[[nredu_1]])[cells,] )
  subobject[[nredu_2]] <- CreateDimReducObject(Embeddings(object[[nredu_2]])[cells,] )
  subobject[["liger"]] <- CreateDimReducObject(bind_reduction, assay = "RNA", key=reduction.name)
  subobject <- RunUMAP(subobject, reduction=reduction.name, dims=1:ncol(bind_reduction), reduction.name =glue("{reduction.name}_UMAP"), verbose=F)
  object@tools[["liger"]] <- subobject
  return(object)

}

liger_run <-   function(object,
                        ASSAY_1="RNA",
                        ASSAY_2="Peaks",
                        nredu_1="RNA_PCA",
                        nredu_2="lsi",
                        dims_1 = 1:50,
                        dims_2 = 1:50,
                        reduction.name = "liger"
){

  library(rliger)
  dataset <- object@tools$dataset
  if(dataset=="cite_30k" | dataset == "cite_elife"){
    object <-liger_cite(object,
               ASSAY_1,
               ASSAY_2,
               nredu_1,
               nredu_2,
               dims_1 ,
               dims_2 ,
               reduction.name)
    return(object)

  }
  shared_atac <- load_object(file=glue("{dataset}/save/shared_mtx.Rds"))
  unshared_atac <- load_object(file=glue("{dataset}/save/unshared_mtx.Rds"))

  shared_atac <- shared_atac[rownames(shared_atac) != "", ]

  if(dataset == "skin"){
    atac_bc <- stringr::str_replace_all(colnames(shared_atac), ",", ".")
    map_vec <- unname(object$rna.bc)
    names(map_vec) <- unname(object$atac.bc)
    rna_bc <- revalue(atac_bc, map_vec, warn_missing=F)
    colnames(shared_atac) <- rna_bc
  }

  rna <- GetAssayData(object, assay="RNA", slot="counts")
  shared_bc <- Reduce(intersect, list(colnames(rna),colnames(shared_atac),colnames(unshared_atac)))
  rna <- rna[, shared_bc]
  shared_atac <- shared_atac[, shared_bc]
  unshared_atac  <- unshared_atac[, shared_bc]
  colnames(rna) <- paste0("rna_", colnames(rna))


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
  subobject <- CreateSeuratObject(counts=GetAssayData(object, assay=ASSAY_1, slot="counts")[1:20, cells], meta.data=object@meta.data[cells, ])
  subobject[[ASSAY_2]] <- CreateAssayObject(counts=GetAssayData(object, assay=ASSAY_2, slot="counts")[1:20, cells])
  subobject[[nredu_1]] <- CreateDimReducObject(Embeddings(object[[nredu_1]])[cells,] )
  subobject[[nredu_2]] <- CreateDimReducObject(Embeddings(object[[nredu_2]])[cells,] )
  subobject[["liger"]] <- CreateDimReducObject(bind_reduction, assay = "RNA", key=reduction.name)

  subobject <- RunUMAP(subobject, reduction=reduction.name, dims=1:ncol(bind_reduction), reduction.name =glue("{reduction.name}_UMAP"), verbose=F)

  object@tools[["liger"]] <- subobject


  return(object)
}


MOJITOO_run <- function(object,
                        ASSAY_1="RNA",
                        ASSAY_2="Peaks",
                        nredu_1="RNA_PCA",
                        nredu_2="lsi",
                        dims_1 = 1:50,
                        dims_2 = 1:50,
                        reduction.name = "MOJITOO"

){

  mute(library(MOJITOO))
  object <- mojitoo(
       object=object,
       reduction.list = list(nredu_1, nredu_2),
       dims.list = list(dims_1, dims_2),
       reduction.name='MOJITOO',
       assay=ASSAY_1
  )

  DefaultAssay(object) <- ASSAY_1
  embedd <- Embeddings(object[[reduction.name]])
  object <- RunUMAP(object, reduction="MOJITOO", reduction.name="MOJITOO_UMAP", dims=1:ncol(embedd), verbose=F)
  return(object)
}

DIABLO_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="Peaks",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      reduction.name = "DIABLO"
){


  library(mixOmics)
  Y = 1:ncol(object) ## set all Y to be distinct ids
  X1 = GetAssayData(object, assay=ASSAY_1, slot="data") ## use normalized data
  X2 = GetAssayData(object, assay=ASSAY_2, slot="data") ## use normalized data

  X <- list(t(X1), t(X2))
  names(X) <- c(ASSAY_1, ASSAY_2)
  design =  matrix(c(0,1,1,0),ncol=2, nrow=2, byrow=TRUE,
                   dimnames = list(names(X), names(X)))
  diablo <- block.plsda(X,Y, design="full", ncomp=30)
  dimension <- cbind(diablo$variates[[ASSAY_1]], diablo$variates[[ASSAY_2]])

  dims <- ncol(object)
  object[[reduction.name]] <- CreateDimReducObject(dimension, assay = "RNA", key=reduction.name)
  object <- RunUMAP(object, reduction=reduction.name, dims=1:dims, reduction.name =glue("{reduction.name}_UMAP"), verbose=F)

  return(object)
}


DIABLORed_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="Peaks",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      reduction.name = "DIABLORed"
){


  library(mixOmics)
  Y = 1:ncol(object) ## set all Y to be
  X1 = Embeddings(object[[nredu_1]])[, dims_1] ## use normalized data
  X2 = Embeddings(object[[nredu_2]])[, dims_2] ## use normalized data

  X <- list(X1, X2)
  names(X) <- c(ASSAY_1, ASSAY_2)
  design =  matrix(c(0,1,1,0),ncol=2, nrow=2, byrow=TRUE,
                   dimnames = list(names(X), names(X)))
  diablo <- block.plsda(X,Y, design=design, ncomp=min(30,  ncol(X1), ncol(X2)))

  dimension <- cbind(diablo$variates[[ASSAY_1]], diablo$variates[[ASSAY_2]])
  dims <- ncol(dimension)
  object[[reduction.name]] <- CreateDimReducObject(dimension, assay = "RNA", key=reduction.name)
  object <- RunUMAP(object, reduction=reduction.name, dims=1:dims, reduction.name =glue("{reduction.name}_UMAP"), verbose=F)

  return(object)
}

symphony_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="Peaks",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      reduction.name = "symphony"
){


  library(mixOmics)
  Y = 1:ncol(object) ## set all Y to be
  X1 = GetAssayData(object, assay=ASSAY_1, slot="scale.data") ## use normalized data
  X2 = GetAssayData(object, assay=ASSAY_2, slot="scale.data") ## use normalized data

  cc_out <- cc(t(X1), t(X2))
  dimension <- cc_out$scores$xscores[, 1:20]

  dims <- ncol(dimension)
  object[[reduction.name]] <- CreateDimReducObject(dimension, assay = "RNA", key=reduction.name)
  object <- RunUMAP(object, reduction=reduction.name, dims=1:dims, reduction.name =glue("{reduction.name}_UMAP"), verbose=F)

  return(object)
}


symphonyRed_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="Peaks",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      reduction.name = "symphonyRed",
                      dataset="cite"
){


  library(mixOmics)
  library(harmony)
  Y = 1:ncol(object) ## set all Y to be
  X1 = Embeddings(object[[nredu_1]])[, dims_1]
  X2 = Embeddings(object[[nredu_1]])[, dims_2]

  cc_out <- cc(X1, X2)
  dimension <- cc_out$scores$xscores#[, 1:20]
  dims <- ncol(dimension)


  if(dataset=="cite_30k"){
      message("30k harmony" )
      dimension <- HarmonyMatrix(dimension, meta_data=object@meta.data, vars_use= "lane",
                               plot_convergence = FALSE, nclust = 100, max.iter.harmony = 20,
                               max.iter.cluster = 20, do_pca = F, verbose = T)
  }

  if(dataset=="cite_elife"){
      message("elife harmony" )
      dimension <- HarmonyMatrix(dimension, meta_data=object@meta.data, vars_use= "group",
                               plot_convergence = FALSE, nclust = 100, max.iter.harmony = 20,
                               max.iter.cluster = 20, do_pca = F, verbose = T)
  }



  object[[reduction.name]] <- CreateDimReducObject(dimension, assay = "RNA", key=reduction.name)
  object <- RunUMAP(object, reduction=reduction.name, dims=1:dims, reduction.name =glue("{reduction.name}_UMAP"), verbose=F)

  return(object)
}

scAI_run <- function(object,
                      ASSAY_1="RNA",
                      ASSAY_2="Peaks",
                      nredu_1="RNA_PCA",
                      nredu_2="lsi",
                      dims_1 = 1:50,
                      dims_2 = 1:50,
                      rerun=FALSE,
                      dataset="kidney"

){

  if(rerun){
      mtx_list <- list(rna = (GetAssayData(object, slot="counts", assay=ASSAY_1)[object@assays[[ASSAY_1]]@var.features, ]),
                       atac = (GetAssayData(object, slot="counts", assay=ASSAY_2)[object@assays[[ASSAY_2]]@var.features, ]))

      mtx_list <- list(rna = (GetAssayData(object, slot="data", assay=ASSAY_1)[object@assays[[ASSAY_1]]@var.features, ]),
                       atac = (GetAssayData(object, slot="data", assay=ASSAY_2)[object@assays[[ASSAY_2]]@var.features, ]))


      X <- list(RNA=mtx_list[[1]], ATAC=mtx_list[[2]])
      labels <- data.frame(labels=object$celltype) # the prior labels of cells
      scAI_outs <- create_scAIobject(raw.data = X, do.sparse = F)
      ## add normalization for scAI
      library(Matrix)
      scAI_outs <- preprocessing(scAI_outs, assay = NULL, minFeatures=0)
      scAI_outs <- addpData(scAI_outs, pdata = labels, pdata.name = "Cell types")

      start_time <- Sys.time()
      scAI_outs <- run_scAI(scAI_outs, K = 20, nrun = 1, do.fast = T)
      end_time <- Sys.time()
      save_object(scAI_outs, file.path(dataset, "save/scAI_outs_5.Rds"))
      save_object(end_time - start_time, file.path(dataset, "save/scAI_time.Rds"))
  }
  scai <- load_object(file=file.path(dataset, "save/scAI_outs.Rds"))
  reduction <-t(scai@fit$H[, colnames(object)])
  object[["scAI"]] <- CreateDimReducObject(embeddings=reduction, key="factor", assay=ASSAY_1)
  object <- RunUMAP(object, reduction.name=glue("scAI_UMAP"), reduction=glue("scAI"), dims=1:20, verbose=F)
  return(object)
}

## get silhouette labels to reductions
getSilhouette <- function(object,
                          nredu="lsi",
                          dims=1:50,
                          method = "euclidean"
){

  if (!("silh" %in% names(object@tools))){
    object@tools$silh <- list()
  }

  objectbak <- object
  if(nredu== "liger"){
    object <- objectbak@tools$liger
  }

  ## remove NA from celltypes
  NA_names <- names(which(is.na(object$celltype)))
  keep_names <- setdiff(names(object$celltype), NA_names)
  celltypes <- object$celltype[keep_names]
  icelltypes <- as.integer(factor(celltypes))


  if(nredu=="WNN"){
    keep_wsnn <- object@graphs$wsnn[keep_names, keep_names]
    silh_obj <- cluster::silhouette(icelltypes, as.dist(1-keep_wsnn))
    #mean(silh_obj[, 'sil_width'])
    object@tools[[glue("silh_{method}")]][["WNN"]] <- silh_obj
  }else{
    m <- Embeddings(object, reduction=nredu)[keep_names, dims]
    if(startsWith(nredu, "MOJITOO")){
      N11 <- function(x) ((2*((x-min(x))/(max(x)-min(x)))) - 1)
      m <- apply(m, 2, N11)
    }
    d <- parDist(m, method=method)
    silh_obj <- cluster::silhouette(icelltypes, d)
    object@tools[[glue("silh_{method}")]][[nredu]] <- silh_obj
    objectbak@tools[[glue("silh_{method}")]][[nredu]] <- silh_obj
  }
  if ("silh" %in% names(object@tools)){
    object@tools$silh <- NULL
  }
  if(nredu == "liger"){
    objectbak@tools$liger <- object
  }else{
    objectbak <- object
  }
  rm(object)
  gc()

  return(objectbak)
}

getRangeClusters <- function(object,
                             ASSAY,
                             reduction,
                             dims,
                             vres=seq(0.1, 2, by=0.1)
){
  message(date(), " reduction: ", reduction)
  DefaultAssay(object) <- ASSAY
  search_step = 0.1
  print(vres)

  objectbak <- object
  if(reduction== "liger"){
    object <- objectbak@tools$liger
  }

  K = length(unique(object$celltype))
  # better keep the cluster res?
  if(reduction=="WNN"){
     #vres <- seq(0.1, 2, by=0.1)
     object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, verbose = F, resolution=vres)
     for(res in vres){
       key=glue("wsnn_res.{res}")
       cluster_name <- glue("{reduction}_clusters.{res}")
       object@meta.data[, cluster_name] <- object@meta.data[, key]
       object@meta.data[, glue("{cluster_name}_{res}")] <- NULL
       object@meta.data[, key] <- NULL
     }

  }else{
       #vres <- seq(0.1, 2, by=0.1)
       object <- FindNeighbors(object, reduction=reduction, dims=dims, verbose=F)
       object <- FindClusters(object, resolution=vres, verbose=F)
       for(res in vres){
         key=glue("{ASSAY}_snn_res.{res}")
         cluster_name <- glue("{reduction}_clusters.{res}")
         object@meta.data[, cluster_name] <- object@meta.data[, key]
         object@meta.data[, key] <- NULL
         object@meta.data[, glue("{cluster_name}_{res}")] <- NULL

       }
  }
  if(reduction == "liger"){
    objectbak@tools$liger <- object
  }
  else{
    objectbak <- object
  }
  rm(object)
  gc()

  return(objectbak)
}

getRangeARI <- function(object,
                        ASSAY,
                        reduction,
                        vres=seq(0.1, 2, by=0.1)
){


  objectbak <- object
  if(reduction== "liger"){
    object <- objectbak@tools$liger
  }

  message(glue("get {reduction}..."))
  ## remove NA from celltypes
  NA_names <- names(which(is.na(object$celltype)))
  keep_names <- setdiff(names(object$celltype), NA_names)
  celltypes <- object$celltype[keep_names]
  icelltypes <- as.integer(factor(celltypes))

  aris = sapply(vres, function(res){
            aricode::ARI(icelltypes, object@meta.data[keep_names, glue("{reduction}_clusters.{res}")])})

  return(aris)
}

## same number of cluster labels
getLabelNClusters <- function(object,
                              ASSAY,
                              reduction,
                              dims
){
  message(date(), " reduction: ", reduction)

  objectbak <- object
  if(reduction== "liger"){
    object <- objectbak@tools$liger
  }
  DefaultAssay(object) <- ASSAY
  search_step = 0.1
  cluster_name <- glue("{reduction}_clusters")

  K = length(unique(object$celltype))
  # better keep the cluster res?
  if(reduction=="WNN"){
     L <- 0.01
     R <- 5
     while(TRUE){
        res  <- ((L + R) / 2)
        object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, verbose = F, resolution=res)
        k = length(unique(object$seurat_clusters))
        key=glue("wsnn_res.{res}")
        print(glue("L{L} R{R} res{res} k{k} K{K}"))
        if( k == K | (abs(L-R) < 0.000001)){
          object@meta.data[, cluster_name] <- object$seurat_clusters
          object@meta.data[, glue("{cluster_name}_{res}")] <- object$seurat_clusters
          object@meta.data[, key] <- NULL
          message("Bingo!")
          break
        }else if(k < K){
          L <- res + search_step
        }
        else if(k > K){
          R <- res - search_step
        }
        if(L>=R ){
           R = R + search_step
           L = max(0.01, L-search_step)
           search_step <- search_step/10
           #stop("binary search failed, use smaller steps")
        }
        object@meta.data[, key] <- NULL
     }

  }else{
    L <- 0.01
    R <- 5
    while(TRUE){
       res  <- ((L + R) / 2)
       object <- FindNeighbors(object, reduction=reduction, dims=dims, verbose=F)
       object <- FindClusters(object, resolution=res, verbose=F)
       key=glue("{ASSAY}_snn_res.{res}")
       k = length(unique(object$seurat_clusters))
       print(glue("L{L} R{R} res{res} k{k} K{K}"))
       if( k == K| (abs(L-R) < 0.000001)){
         object@meta.data[, cluster_name] <- object@meta.data$seurat_clusters
         object@meta.data[, glue("{cluster_name}_{res}")] <- object@meta.data$seurat_clusters
         print(str(object$seurat_clusters))
         #object@meta.data[, cluster_name] <- object@meta.data[, key]
         #object@meta.data[, glue("{cluster_name}_{res}")] <- object@meta.data[, key]

         object@meta.data[, key] <- NULL
         message("Bingo!")
         break
       }else if(k < K){
         L <- res + search_step
       }
       else if(k > K){
         R <- res - search_step
       }
       if(L >=R ){
           R = R + search_step
           L = max(0.01, L-search_step)
           search_step <- search_step/10
          #stop("binary search failed, use smaller steps")
       }
       object@meta.data[, key] <- NULL
    }
  }
  if(reduction == "liger"){
    objectbak@tools$liger <- object
  }else{
    objectbak <- object
  }
  rm(object)
  gc()
  return(objectbak)
}

## get best nmi & ari on one specific cluster resolution
getBestCluster <- function(object,
               ASSAY,
               reduction,
               dims
){
  resolutions = seq(0.1, 2, 0.1)
  if(reduction == "WNN"){
    object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, verbose = F, resolution=reductions)
    nmis = sapply(resolutions, function(res){
              aricode::NMI(object$celtype, object@meta.data[, glue("wsnn_res.{res}")])})
    idx <- which.max(nmis)
    nmi <- nmis[idx]
    nmi_keep <- resolutions[idx]

    aris = sapply(resolutions, function(res){
              aricode::ARI(object$celtype, object@meta.data[, glue("wsnn_res.{res}")])})
    idx <- which.max(aris)
    ari <- aris[idx]
    ari_keep <- resolutions[idx]

    sapply(resolutions, function(res)object@meta.data[, glue("wsnn_res.{res}")] <- NULL)

  }else{
    object <- FindNeighbors(object, reduction=reduction, dims=dims, verbose=F)
    object <- FindClusters(object, resolution=resolutions, verbose=F)
    nmis = sapply(resolutions, function(res){
              aricode::NMI(object$celtype, object@meta.data[, glue("{ASSAY}_snn_res.{res}")])})
    idx <- which.max(nmis)
    nmi <- nmis[idx]
    nmi_keep <- resolutions[idx]

    aris = sapply(resolutions, function(res){
              aricode::ARI(object$celtype, object@meta.data[, glue("{ASSAY}_snn_res.{res}")])})
    idx <- which.max(aris)
    ari <- aris[idx]
    ari_keep <- resolutions[idx]


    sapply(resolutions, function(res){if(res != keep){
      object@meta.data[, glue("{ASSAY}_snn_res.{res}")] <- NULL
    }})
  }
  if(!("ARI" %in% names(object@tools))){
    object@tools$ARI <- list()
  }
  if(!("NMI" %in% names(object@tools))){
    object@tools$NMI <- list()
  }
  object@tools$NMI[reduction] <- data.frame("res"=nmi_keep, nmi=nmi)
  object@tools$ARI[reduction] <- data.frame("res"=ari_keep, ari=ari)
  return(object)
}



Seurat2Scanpy <- function(file=NULL){
  message(date(), " loading SeuratObject...")
  object <- load_object(file=file)
  midfile=paste0(tools::file_path_sans_ext(file), ".h5Seurat")
  message(date(), " To h5seurat...")
  SaveH5Seurat(object, filename = midfile)
  message(date(), " To h5ad...")
  Convert(midfile, dest = "h5ad")
  message(date(), " Done.")
}

