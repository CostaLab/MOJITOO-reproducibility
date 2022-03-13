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
mute(library(Rcpp))
mute(library(RcppArmadillo))
mute(library(cluster))
mute(library(parallelDist))
mute(library(S4Vectors))
mute(library(Matrix))
mute(library(assertthat))
mute(library(CCA))
mute(library(SeuratDisk))
mute(library(knn.covertree))
mute(library(harmony))
mute(library(reticulate))

getSettings <- function(dataset){
    SL <- SimpleList(
      ASSAY_1="RNA",
      ASSAY_2="Peaks",
      Redu_1 = "RNA_PCA",
      Redu_2 = "lsi",
      dims_1 = 1:50,
      dims_2 = 2:50
    )

if(dataset == "TEA" | dataset == "ASAP_BM"){
    SL <- SimpleList(
      ASSAY_1="RNA",
      ASSAY_2="atac",
      ASSAY_3="adt",
      Redu_1 = "rpca",
      Redu_2 = "lsi",
      Redu_3 = "apca",
      dims_1 = 1:50,
      dims_2 = 2:50,
      dims_3 = 1:30
      )
  }


  SL$dims_add = 1:(length(SL$dims_1) + length(SL$dims_2) + length(SL$dims_3))
  SL$dims_min = 1:(min(length(SL$dims_1) , length(SL$dims_2), length(SL$dims_3) ))
  SL$dims_bind = 1:(3*min(length(SL$dims_1) , length(SL$dims_2), length(SL$dims_3)))

  return(SL)
}


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
                         nredu_1 = "rpca",
                         nredu_2 = "lsi",
                         nredu_3 = "apca",
                         dims_1 = 1:50,
                         dims_2 = 1:50,
                         dims_3 = 1:30,
                         nredu="lsi",
                         dims=1:50,
                         method = "euclidean",
                         cor_method = "pearson"
){

  if (!("DEmbed" %in% names(object@tools))){
    object@tools$DEmbed <- list()
  }
  if(!(nredu_1  %in%names(object@tools$DEmbed))){
    object@tools$DEmbed[[nredu_1]] <- getDistanceVector(object, nredu_1, dims_1, method)
  }
  if(!(nredu_2  %in%names(object@tools$DEmbed))){
    object@tools$DEmbed[[nredu_2]] <- getDistanceVector(object, nredu_2, dims_2, method)
  }
  if(!(nredu_3  %in%names(object@tools$DEmbed))){
    object@tools$DEmbed[[nredu_3]] <- getDistanceVector(object, nredu_3, dims_3, method)
  }

  vec <- getDistanceVector(object, nredu, dims, method=method)

  if(cor_method == "pearson"){
    Correlation <- c(cor(object@tools$DEmbed[[nredu_1]], vec),
                     cor(object@tools$DEmbed[[nredu_2]], vec),
                     cor(object@tools$DEmbed[[nredu_3]], vec)
                  )

    object@tools[[glue("structure_{method}")]][[nredu]] <- Correlation
  }else if(cor_method == "spearman"){
    Correlation <- c(cor(object@tools$DEmbed[[nredu_1]], vec, method="spearman"),
                     cor(object@tools$DEmbed[[nredu_2]], vec, method="spearman"),
                     cor(object@tools$DEmbed[[nredu_3]], vec, method="spearman")
                  )
    object@tools[[glue("structure_{method}_spearman")]][[nredu]] <- Correlation
  }

  return(object)
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
                          ASSAY_list=list("RNA", "atac", "adt"),
                          nredu_list=list("rpca", "lsi", "apca"),
                          dims_list=list(1:50, 1:50, 1:30),
                          name="yyy"
){

  integration = FALSE
  if(name == "TEA"){
    integration = FALSE
  }else if(name =="ASAP_BM"){
    integration = TRUE
  }

  if((name == "cite_30k")){
    #do nothing
  }else{

    ## rna
    DefaultAssay(object) <- ASSAY_list[[1]]
    object <-  NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
    object  <- FindVariableFeatures(object, nfeatures=3000)
    object  <- ScaleData(object, nfeatures=3000)
    if(integration){
      object <- RunPCA(object, npcs=50, reduction.name="rpca.raw", verbose=F)
      object <- RunHarmony(object, group.by.vars = 'stim',
                           reduction = 'rpca.raw', assay.use = 'RNA',
                           project.dim = FALSE,  reduction.save = "rpca")
    }else{
      object <- RunPCA(object, npcs=50, reduction.name="rpca", verbose=F)
    }

    ## atac
    DefaultAssay(object) <- ASSAY_list[[2]]
    object <- RunTFIDF(object)
    object <- FindTopFeatures(object, min.cutoff = 'q0')
    if(integration){
      object <- RunSVD(object, reduction.name = "lsi.raw")
      object <- RunHarmony(object, group.by.vars = 'stim',
                           reduction = 'lsi.raw', assay.use = 'atac',
                           project.dim = FALSE,  reduction.save = "lsi")
    }else{
      object <- RunSVD(object, reduction.name = "lsi")
    }


    ## apca
    DefaultAssay(object) <- ASSAY_list[[3]]
    VariableFeatures(object) <- rownames(object[["adt"]])
    object <- NormalizeData(object, normalization.method = 'CLR', margin = 2) %>%  ScaleData()
    if(integration){
      object <- RunPCA(object, reduction.name = 'apca.raw', verbose=F, npcs=30)
      object <- RunHarmony(object, group.by.vars = 'stim',
                           reduction = 'apca.raw', assay.use = 'adt',
                           project.dim = FALSE,  reduction.save = "apca")
    }else{
      object <-  RunPCA(object, reduction.name = 'apca', verbose=F, npcs=30)
    }
    DefaultAssay(object) <- "RNA"
    object <- RunUMAP(object, reduction.name=glue("{nredu_list[[1]]}_UMAP"), reduction=glue("{nredu_list[[1]]}"), dims=dims_list[[1]], verbose=F)
    object <- RunUMAP(object, reduction.name=glue("{nredu_list[[2]]}_UMAP"), reduction=glue("{nredu_list[[2]]}"), dims=dims_list[[2]], verbose=F)
    object <- RunUMAP(object, reduction.name=glue("{nredu_list[[3]]}_UMAP"), reduction=glue("{nredu_list[[3]]}"), dims=dims_list[[3]], verbose=F)
  }

  DefaultAssay(object) <- ASSAY_list[[2]]
  object <- FindVariableFeatures(object, nfeatures=5000)
  object <- ScaleData(object, features=VariableFeatures(object))

  return(object)
}

MOFARed_run <- function(object,
                         ASSAY_list=list("RNA", "atac", "adt"),
                         nredu_list=list("rpca", "lsi", "apca"),
                         dims_list=list(1:50, 1:50, 1:30),
                         name = "skin"
){

  mute(library(MOFA2))
  ASSAY_1=ASSAY_list[[1]]
  ASSAY_2=ASSAY_list[[2]]
  ASSAY_3=ASSAY_list[[3]]

  nredu_1 = nredu_list[[1]]
  nredu_2 = nredu_list[[2]]
  nredu_3 = nredu_list[[3]]

  dims_1 = dims_list[[1]]
  dims_2 = dims_list[[2]]
  dims_3 = dims_list[[3]]

   mofa <- create_mofa_from_matrix(
    list(t(Embeddings(object[[nredu_1]])[, dims_1]),
         t(Embeddings(object[[nredu_2]])[, dims_2]),
         t(Embeddings(object[[nredu_3]])[, dims_3])
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
  colnames(mofaUMAP) <- paste0("UMAP_", 1:2)
  mofaUMAP <- mofaUMAP[colnames(object), ]

  object[["MOFARed_UMAP"]] <- CreateDimReducObject(embeddings=as.matrix(mofaUMAP), key="mofared", assay=ASSAY_1)

  #factors <- 1:get_dimensions(mofa)[["K"]]
  Z <- get_factors(mofa, factors = factors, groups = "all")[[1]]
  Z <- Z[colnames(object), ]
  object[["MOFARed"]] <- CreateDimReducObject(embeddings=as.matrix(Z), key="factor", assay=ASSAY_1)

  return(object)

}

MOFA_run <- function(object,
                    ASSAY_list ,
                    name = "skin"
){

  mute(library(MOFA2))
  ASSAY_1=ASSAY_list[[1]]
  ASSAY_2=ASSAY_list[[2]]
  ASSAY_3=ASSAY_list[[3]]

  mofa <- create_mofa(object, assays = c(ASSAY_1,ASSAY_2, ASSAY_3))
  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- 15
  mofa <- prepare_mofa(mofa,
    model_options = model_opts
  )
  mofa <- run_mofa(mofa)

  plot_factor_cor(mofa)
  factors <- 1:get_dimensions(mofa)[["K"]]
  if(name == "pbmc"){
    factors <- factors[!factors%in%c(4,7)]
  }

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


schema_run <- function(object,
                       ASSAY_list = list(),
                       dims=1:30,
                       name="skin"
){

  embedd <- read.csv(glue("{dataset}/save/schema.csv"), row.names = 1)
  object[["schema"]] <- CreateDimReducObject(embeddings=as.matrix(embedd[colnames(object), ]), key="schema", assay=ASSAY_list[[1]])
  object <- RunUMAP(object, reduction="schema", reduction.name=glue("schema_UMAP"), dims=dims, verbose=F)
  return(object)
}

MOJITOO_run <- function(object,
                         ASSAY_list=list("RNA", "atac", "adt"),
                         nredu_list=list("rpca", "lsi", "apca"),
                         dims_list=list(1:50, 1:50, 1:30),
                         reduction.name = "MOJITOO"

){

  mute(library(MOJITOO))
  object <- mojitoo(
       object=object,
       reduction.list=nredu_list,
       dims.list=dims_list,
       reduction.name='MOJITOO',
       #assay=ASSAY_list[[1]]
  )

  DefaultAssay(object) <- ASSAY_list[[1]]
  embedd <- Embeddings(object[[reduction.name]])
  object <- RunUMAP(object, reduction="MOJITOO", reduction.name="MOJITOO_UMAP", dims=1:ncol(embedd), verbose=F)
  return(object)
}

DIABLORed_run <- function(object,
                         ASSAY_list=list("RNA", "atac", "adt"),
                         nredu_list=list("rpca", "lsi", "apca"),
                         dims_list=list(1:50, 1:50, 1:30),
                         reduction.name = "DIABLORed"

){

  dims_1 <- dims_list[[1]]
  dims_2 <- dims_list[[2]]
  dims_3 <- dims_list[[3]]

  library(mixOmics)
  Y = 1:ncol(object) ## set all Y to be
  X1 = Embeddings(object = object[[nredu_list[[1]]]])[, dims_1]
  X2 = Embeddings(object = object[[nredu_list[[2]]]])[, dims_2]
  X3 = Embeddings(object = object[[nredu_list[[3]]]])[, dims_3]

  X <- list(X1, X2, X3)
  names(X) <- c(ASSAY_list[[1]], ASSAY_list[[2]], ASSAY_list[[3]])
  diablo <- block.plsda(X, Y, design='full', ncomp=30)
  dimension <- cbind(diablo$variates[[ASSAY_list[[1]]]],
                     diablo$variates[[ASSAY_list[[2]]]],
                     diablo$variates[[ASSAY_list[[3]]]])
  dims <- ncol(dimension)
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

      X <- list(RNA=mtx_list[[1]], ATAC=mtx_list[[2]])
      labels <- data.frame(labels=object$celltype) # the prior labels of cells
      scAI_outs <- create_scAIobject(raw.data = X, do.sparse = F)
      scAI_outs <- preprocessing(scAI_outs, assay = NULL)
      scAI_outs <- addpData(scAI_outs, pdata = labels, pdata.name = "Cell types")

      start_time <- Sys.time()
      scAI_outs <- run_scAI(scAI_outs, K = 20, nrun = 1, do.fast = T)
      end_time <- Sys.time()
      save_object(scAI_outs, file.path(dataset, "save/scAI_outs.Rds"))
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
    cluster_name <- glue("{nredu}_clusters")
    m <- Embeddings(object, reduction=nredu)[keep_names, dims]
    if(nredu == "MOJITOO"){# MOJITOO silhouette score normalize it to -1 to 1
      N11 <- function(x) ((2*((x-min(x))/(max(x)-min(x)))) - 1)
      m <- apply(m, 2, N11)
    }
    d <- parDist(m, method=method)
    silh_obj <- cluster::silhouette(icelltypes, d)
    object@tools[[glue("silh_{method}")]][[nredu]] <- silh_obj
  }
  if ("silh" %in% names(object@tools)){
    object@tools$silh <- NULL
  }
  return(object)
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
  return(object)
}

getRangeARI <- function(object,
                        ASSAY,
                        reduction,
                        vres=seq(0.1, 2, by=0.1)

){
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
       }
       object@meta.data[, key] <- NULL
    }
  }
  return(object)
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
