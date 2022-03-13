library(cli)
library(Matrix)
library(stringi)
library(stringr)
library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)

library(org.Hs.eg.db)
library(org.Mm.eg.db)

library(EnsDb.Mmusculus.v75)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(liger)
library(rliger)


unshared_bin_mtx <- function(peaks, organ="hg38", cpus=30){
  library(foreach)
  library(doParallel)
  library(tictoc)
  tic()
  registerDoParallel(cores=cpus)


  peak_ranges <- lapply(rownames(peaks),
                        function(x)c(x,rev(stri_reverse(str_split_fixed(stri_reverse(x), stri_reverse("-"), 3)))))
  df_peak <- as.data.frame.matrix(do.call(rbind, peak_ranges))
  names(df_peak) <- c("name", "chr", "start", "end")
  if(organ=="mm9" | organ=="mm10"){
    keep_chromosome <- c(1:19, "X", "Y")
    ensdb <- switch(organ, "mm9"=EnsDb.Mmusculus.v75, "mm10"=EnsDb.Mmusculus.v79)
  }else if(organ=="hg19"| organ=="hg38"){
    keep_chromosome <- c(1:22, "X", "Y")
    ensdb <- switch(organ, "hg19"= EnsDb.Hsapiens.v75, "hg38"=EnsDb.Hsapiens.v86)
  }else{
    stop("wrong organ")
  }
  df_peak <- df_peak %>% dplyr::filter(chr %in% paste0("chr", keep_chromosome))

  gene.coords <- genes(ensdb)
  gene.coords <-  keepSeqlevels(gene.coords, keep_chromosome, pruning.mode="coarse")

  #seqlevelsStyle(gene.coords) <- "UCSC"
  ## https://github.com/Bioconductor/GenomeInfoDb/issues/27
  ucsc.levels <- str_replace(string=paste("chr",seqlevels(gene.coords),sep=""), pattern="chrMT", replacement="chrM")
  seqlevels(gene.coords) <- ucsc.levels
  gr_tile <-  tileGenome(seqinfo(gene.coords), tilewidth=100000, cut.last.tile.in.chrom=T)

  if(organ=="mm9" | organ=="mm10"){
    black <- read.table("data/cellranger-reference/Blacklists/mm10-blacklist.v2.bed.gz", sep="\t")
  }else if(organ=="hg19"){
    black <- read.table("data/cellranger-reference/Blacklists/hg19-blacklist.v2.bed.gz", sep="\t")
  }else if(organ=="hg38"){
    black <- read.table("data/cellranger-reference/Blacklists/hg38-blacklist.v2.bed.gz", sep="\t")
  }

  names(black) <- c("chr", "start", "end", "reason" )
  #black$chr <-  stringr::str_sub(black$chr, start=4)
  gr_black <- makeGRangesFromDataFrame(black)
  gr_ov <- findOverlaps(gr_tile, gr_black)
  assertthat::assert_that(length(queryHits(gr_ov) )> 0)
  gr_tile <- gr_tile[-queryHits(gr_ov), ]

  gr_tile$name <- paste(as.character(seqnames(gr_tile)),
                        as.character(start(ranges(gr_tile))),
                        as.character(end(ranges(gr_tile))),
                        sep="-")
  gr_peak <- makeGRangesFromDataFrame(df_peak, keep.extra.columns=T)
  fo <- findOverlaps(gr_tile, gr_peak)

  tic()
  qh <-queryHits(fo)
  uqh <- unique(qh)
  sh <-subjectHits(fo)
  ret_list = foreach(i=1:length(uqh)) %dopar%{
    if(i %% 1000 == 0)
      print(i)
    idx <- uqh[i]## get the first one
    hit_idxs <- which(qh == idx)
    hits <- sh[hit_idxs]
    peaks_idx <- gr_peak[hits, ]$name
    if(length(peaks_idx)>1){
       return(colSums((peaks[peaks_idx, ]>0)))
    }else{
      return( as.numeric(peaks[peaks_idx, ]>0))
    }
  }
  toc()
  mtx <- as(do.call(rbind, ret_list), "dgCMatrix")
  rownames(mtx) <- gr_tile[uqh, ]$name
  rsums <- rowSums(mtx)
  row0 <- which(rsums == 0)
  if(length(row0) > 0){
    mtx <- mtx[-row0, ]
  }

  message(date(), " finished")
  return(mtx)
}


shared_counts <- function(fragment_file, dataset='pbmc', organ="hg38", sort_fragment=F){
  ##shell
  gene_file <- glue("tools/save/{organ}_genes.bed")
  promoter_file <- glue("tools/save/{organ}_promoter.bed")



  if(sort_fragment){
    message("sorting fragment")
    if(endsWith(fragment_file, "gz")){
      cmd <- glue("zcat {fragment_file} | sort -k1,1 -k2,2n -k3,3n  > {dataset}/save/atac_fragments.sort.bed")
    }else{
      cmd <- glue("sort -k1,1 -k2,2n -k3,3n {fragment_file} > {dataset}/save/atac_fragments.sort.bed")
    }
    #cmd <- glue("bedtools sort -i {fragment_file} > {dataset}/save/atac_fragments.sort.bed")
    system(cmd)
  }else{
    cmd <- glue("cp  {fragment_file}  {dataset}/save/atac_fragments.sort.bed")
    system(cmd)

  }

  cmd1 = glue("sort -k 1,1 -k2,2n -k3,3n {gene_file} > {dataset}/save/{organ}_genes.bed")
  cmd2 = glue("sort -k 1,1 -k2,2n -k3,3n {promoter_file} >{dataset}/save/{organ}_promoter.bed")
  system(cmd1)
  system(cmd2)


  cmd1 = glue('bedmap --ec --delim "\t" --echo --echo-map-id {dataset}/save/{organ}_promoter.bed {dataset}/save/atac_fragments.sort.bed > {dataset}/save/atac_promoters_bc.bed')
  cmd2 = glue('bedmap --ec --delim "\t" --echo --echo-map-id {dataset}/save/{organ}_genes.bed {dataset}/save/atac_fragments.sort.bed > {dataset}/save/atac_genes_bc.bed')


  message(date(), " mapping promoter")
  system(cmd1)
  message(date(), " mapping gene")
  system(cmd2)


  genes.bc <- read.table(file = glue("{dataset}/save/atac_genes_bc.bed"), sep = "\t", as.is = c(4,7), header = FALSE)
  promoters.bc <- read.table(file = glue("{dataset}/save/atac_promoters_bc.bed"), sep = "\t", as.is = c(4,7), header = FALSE)

  bc <- genes.bc[,7]
  bc_split <- strsplit(bc,";")
  bc_split_vec <- unlist(bc_split)
  bc_unique <- unique(bc_split_vec)
  bc_counts <- table(bc_split_vec)
  bc_filt <- names(bc_counts)[bc_counts > 1500]
  barcodes <- bc_filt

  gene.counts <- makeFeatureMatrix(genes.bc, barcodes)
  promoter.counts <- makeFeatureMatrix(promoters.bc, barcodes)

  gene.counts <- gene.counts[order(rownames(gene.counts)),]
  promoter.counts <- promoter.counts[order(rownames(promoter.counts)),]
  rownames(promoter.counts) %in% rownames(gene.counts)
  mrows <- match(rownames(promoter.counts), rownames(gene.counts))
  sum_counts <- gene.counts
  sum_counts[mrows, ] <- sum_counts[mrows, ]  + promoter.counts


  message(date(), " finished")
  return(sum_counts)

}

