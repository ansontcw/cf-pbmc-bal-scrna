# source code of scds package

bcds <- function(sce, ntop=500, srat=1, verb=FALSE, retRes=FALSE,
                 nmax="tune", varImp=FALSE, estNdbl = FALSE){
  
  #- check first argument (sce)
  if(!is(sce,"SingleCellExperiment"))
    stop('First argument (sce) needs to be of class SingleCellExperiment')
  if( !("counts" %in% names(assays(sce))) )
    stop('First argument (sce) needs to have a "counts" assay')
  if(!is(counts(sce),"sparseMatrix"))
    counts(sce) = Matrix::Matrix(counts(sce),sparse=TRUE)
  
  #- select variable genes
  if(verb) message("-> selecting genes\n")
  ind1            = Matrix::rowSums(counts(sce)>0)>0.01*ncol(sce)
  lc              = Matrix::t(log1p(counts(sce)[ind1,]))
  lc              = lc/Matrix::rowMeans(lc)
  vrs             = apply(lc,2,stats::var)
  hvg             = order(vrs,decreasing=TRUE)[seq_len(ntop)]
  lc              = lc[,hvg]
  lc              = lc/Matrix::rowMeans(lc)
  
  #- "simulated" doublets
  if(verb) message("-> simulating doublets\n")
  p1  = sample(seq_len(ncol(sce)),srat*ncol(sce),replace=TRUE)
  p2  = sample(seq_len(ncol(sce)),srat*ncol(sce),replace=TRUE)
  lc2 = Matrix::t(log1p(counts(sce)[ind1,p1][hvg,] + counts(sce)[ind1,p2][hvg,]))
  lc2 = lc2/Matrix::rowMeans(lc2)
  X   = rbind(lc,lc2)
  
  #- learn classifier
  if(verb) message("-> training classifier\n")
  colnames(X) = paste("GEN",seq_len(ncol(X)),sep="_") #- ranger picky
  y           = c(rep(-1,nrow(lc)),rep(1,nrow(lc2)))
  mm  = xgb.DMatrix(X,label=(y+1)/2)
  #- fixed rounds:
  if(nmax != "tune"){
    res    = NA
    varImp = FALSE
    retRes = FALSE
    pre    = xgboost(mm,nrounds=nmax,tree_method="hist",
                     nthread = 16, early_stopping_rounds = 2, subsample=0.5,
                     objective = "binary:logistic",verbose=0)
    sce$bcds_score = stats::predict(pre, newdat= mm[seq_len(ncol(sce)),])
    #- get doublet calls
    if(estNdbl){
      dbl.pre = stats::predict(pre, newdat= mm[seq(ncol(sce)+1,nrow(X)),])
      est_dbl = get_dblCalls_ALL(sce$bcds_score,dbl.pre,rel_loss=srat)
      if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
      metadata(sce)$bcds$ndbl = est_dbl
      metadata(sce)$bcds$sim_scores = dbl.pre
      sce$bcds_call = sce$bcds_score >= est_dbl["balanced","threshold"]
    }
    #- learning rounds with CV:
  } else {
    res = xgb.cv(data =mm, nthread = 16, nrounds = 500, objective = "binary:logistic",
                 nfold=5,metrics=list("error"),prediction=TRUE,
                 early_stopping_rounds=2, tree_method="hist",subsample=0.5,verbose=0)
    ni  = res$best_iteration
    ac  = res$evaluation_log$test_error_mean[ni] + 1*res$evaluation_log$test_error_std[ni]
    ni  = min(which( res$evaluation_log$test_error_mean <= ac  ))
    nmax = ni
    pre = xgboost(mm,nrounds=nmax,tree_method="hist",
                  nthread = 16, early_stopping_rounds = 2, subsample=0.5,
                  objective = "binary:logistic",verbose=0)
    sce$bcds_score = res$pred[seq_len(ncol(sce))]
    #- get doublet calls
    if(estNdbl){
      dbl.pre = stats::predict(pre, newdat= mm[seq(ncol(sce)+1,nrow(X)),])
      est_dbl = get_dblCalls_ALL(sce$bcds_score,dbl.pre,rel_loss=srat)
      if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
      metadata(sce)$bcds$ndbl = est_dbl
      metadata(sce)$bcds$sim_scores = dbl.pre
      sce$bcds_call = sce$bcds_score >= est_dbl["balanced","threshold"]
      
    }
  }
  #- variable importance
  if(varImp){
    if(verb) message("-> calculating variable importance\n")
    vimp = xgb.importance(model=pre)
    vimp$col_index = match(vimp$Feature,colnames(X))
  }
  #- result
  if(retRes){
    hvg_bool                    = (seq_len(nrow(sce))) %in% which(ind1)
    hvg_bool[hvg_bool]          = (seq_len(sum(hvg_bool))) %in% hvg
    hvg_ord                     = rep(NA,nrow(sce))
    hvg_ord[hvg_bool]           = which(ind1)[hvg]
    rowData(sce)$bcds_hvg_bool  = hvg_bool
    rowData(sce)$bcds_hvg_ordr  = hvg_ord
    metadata(sce)$bcds_res_cv   = res
    metadata(sce)$bcds_res_all  = pre
    metadata(sce)$bcds_nmax     = nmax
    if(varImp){
      vimp$gene_index             = hvg_ord[hvg_bool][vimp$col_index]
      metadata(sce)$bcds_vimp     = vimp[seq_len(100),-c(1,5)]
    }
  }
  
  if(verb) message("-> done.\n\n")
  return(sce)
  
}


get_dblCalls_ROC <- function(scrs_real, scrs_sim, rel_loss=1){
  #=============================================================
  
  p    = length(scrs_sim)/(length(scrs_sim)+length(scrs_real))
  rc   = roc(response=c(rep(0,length(scrs_real)),rep(1,length(scrs_sim))),predictor=c(scrs_real,scrs_sim))
  sens = rc$sensitivities
  spec = rc$specificities
  #- use youden to find ROC cutoff
  r       = p/rel_loss/(1-p)
  cut_ind = which.max(sens+r*spec-1)
  thresh  = rc$thresholds[cut_ind]
  
  ndbl    = sum(scrs_real >= thresh)
  fram    = mean(scrs_sim < thresh) #- fraction of sim doublets missed
  
  res = c(ndbl,thresh,fram)
  names(res) = c("number","threshold","fnr") #- false negative rate
  return(res)
}

#' Derive doublet calls from doublset scores
#'
#' Given score vectors for real data and artificial doubles, derive doublet calls based on determining doublet score cutoffs.
#'
#' @param scrs_real numeric vector, the scores for the real/original data
#' @param scrs_sim numeric vector, the scores for the artificial doublets
#' @param type character or numeric, describes how the score threshold for calling doublets is determined. Either \code{"balanced"} or a number between zero and one that indicates the fraction of artificial doublets missed when making calls. Default: \code{"balanced"}.
#' @return numeric, vector containing the (estimated) number of doublets, the score threshold and the fraction of artificial doublets missed (false negative rate, of sorts)
#' @importFrom stats optimize ecdf uniroot quantile
get_dblCalls_dist <- function(scrs_real,scrs_sim, type="balanced"){
  #==================================================================
  
  #- do "balanced errpr"
  
  if(type == "balanced"){
    #======================
    
    es  = ecdf(scrs_sim)
    er  = ecdf(scrs_real)
    rtf = function(thresh) 1-er(thresh) #- right tail; decreases mono with arg
    ltf = function(thresh) es(thresh)   #- left tail; increases mono with arg
    
    res        = uniroot(function(x)ltf(x) - rtf(x), c(min(scrs_real),max(scrs_real)))
    res_val    = res$root
  } else {
    if(!is.numeric(type)) stop("invalid type argument\n")
    res_val = quantile(scrs_sim,prob=type)
  }
  res_ndl    = sum(scrs_real >= res_val)
  res_fnr    = mean(scrs_sim < res_val)
  res        = c(res_ndl,res_val,res_fnr)
  names(res) = c("number","threshold","fnr") #- false negative rate
  return(res)
}

#' Wrapper for getting doublet calls
#'
#' @param scrs_real numeric vector, the scores for the real/original data
#' @param scrs_sim numeric vector, the scores for the artificial doublets
#' @param rel_loss numeric scalar, relative weight of a false positive classification compared with a false negative. Default:1 (same loss for fp and fn).
#' @return numeric, matrix containing the (estimated) number of doublets, the score threshold and the fraction of artificial doublets missed (false negative rate, of sorts) as columns and four types of estimating: "youden", "balanced" and a false negative rate of artificial doublets of 0.1 and 0.01, respecitvely.
get_dblCalls_ALL <- function(scrs_real,scrs_sim,rel_loss=1){
  #=========================================================
  est_dbl = rbind( get_dblCalls_ROC( scrs_real,scrs_sim,rel_loss),
                   get_dblCalls_dist(scrs_real,scrs_sim,"balanced"),
                   get_dblCalls_dist(scrs_real,scrs_sim,0.1),
                   get_dblCalls_dist(scrs_real,scrs_sim,0.01))
  rownames(est_dbl)  = c("youden","balanced","0.1","0.01")
  return(est_dbl)
}

#' Find doublets/multiplets in UMI scRNA-seq data;
#'
#' Annotates doublets/multiplets using co-expression based approach
#'
#' @param sce single cell experiment (\code{SingleCellExperiment}) object to analyze; needs \code{counts} in assays slot.
#' @param ntop integer, indimessageing number of top variance genes to consider. Default: 500
#' @param binThresh integer, minimum counts to consider a gene "present" in a cell. Default: 0
#' @param verb progress messages. Default: FALSE
#' @param retRes logical, whether to return gene pair scores & top-scoring gene pairs? Default: FALSE.
#' @param estNdbl logical, should the numer of doublets be estimated from the data. Enables doublet calls. Default:FALSE. Use with caution.
#' @return sce input sce object \code{SingleCellExperiment} with doublet scores added to colData as "cxds_score" column.
#' @importFrom Matrix Matrix rowSums rowMeans t
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assay<- assays assays<- assayNames rowData rowData<- colData colData<-
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom stats pbinom
#' @export
#' @examples
#' data("sce_chcl")
#' ## create small data set using only 100 cells
#' sce_chcl_small = sce_chcl[, 1:100]
#' sce_chcl_small = cxds(sce_chcl_small)
cxds <- function( sce, ntop=500, binThresh=0, verb=FALSE, retRes=FALSE, estNdbl = FALSE){
  #===========================================================================================
  
  #- check sce argument
  if(!is(sce,"SingleCellExperiment"))
    stop('First argument (sce) needs to be of class SingleCellExperiment')
  if( !("counts" %in% names(assays(sce))) )
    stop('First argument (sce) needs to have a "counts" assay')
  
  #- 1. Binarize counts
  if(verb) message("\n-> binarizing counts\n")
  B = counts(sce) > binThresh
  
  #- 2. Select high-variance genes
  if(verb) message("-> selecting genes\n")
  ps  = Matrix::rowMeans(B)
  vrs = ps*(1-ps)
  hvg = order(vrs,decreasing=TRUE)[seq_len(ntop)]
  Bp  = B[hvg,]
  
  #- score all pairs
  #- Simple model for (01) and (10) observations per gene pair
  if(verb) message("-> scoring gene pairs\n")
  ps    = Matrix::rowMeans(Bp)
  prb   = outer(ps,1-ps)
  prb   = prb + t(prb)
  obs   = Bp %*% (1-Matrix::t(Bp))
  obs   = obs + Matrix::t(obs)
  #- log p-vals for the observed numbers of 01 and 10
  S = stats::pbinom(as.matrix(obs)-1,prob=prb,size=ncol(Bp),
                    lower.tail=FALSE,log=TRUE)
  #- scores
  if(verb) message("-> calcluating cell scores\n")
  scores = Matrix::colSums(Bp * (S%*%Bp))
  sce$cxds_score = -scores
  #- estimate number of doublets
  if(estNdbl){
    if(verb) message("-> estimating number of doublets\n")
    nsamp       = ncol(sce)
    p1          = sample(seq_len(ncol(sce)),nsamp,replace=TRUE)
    p2          = sample(seq_len(ncol(sce)),nsamp,replace=TRUE)
    Bp_sim      = Bp[,p1] + Bp[,p2] #- still logical
    scores_sim  = as.numeric(Matrix::colSums(Bp_sim * (S%*%Bp_sim)))
    est_dbl     = get_dblCalls_ALL(-scores,-scores_sim, rel_loss=1)
    if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
    metadata(sce)$cxds$ndbl       = est_dbl
    metadata(sce)$cxds$sim_scores = -scores_sim
    sce$cxds_call = sce$cxds_score >= est_dbl["balanced","threshold"]
  }
  
  if(retRes){
    if(verb) message("-> prioritizing gene pairs\n")
    res = list(scores=-scores, S=-S, hvg=hvg, binThresh=binThresh)
    #- rank gene pairs by wighted average doublet contribution
    tmp = Matrix::t(Matrix::t(Bp) *(-scores))
    tmp = tmp %*% Matrix::t(Bp)
    tmp = tmp * S
    S = as.matrix(tmp)
    
    colnames(S) = seq_len(ncol(S))
    rownames(S) = colnames(S)
    topPrs = matrix(NA,nrow=100,ncol=2)
    for(i in seq_len(100)){
      pr = which(S == min(S) ,arr.ind=TRUE)[c(1,2)]
      topPrs[i,] = as.integer(colnames(S)[pr])
      S = S[-pr,]
      S = S[,-pr]
    }
    res$topPairs = topPrs
    
    #- put result in sce
    if(verb) message("-> done.\n\n")
    hvg_bool                     = (seq_len(nrow(sce))) %in% res$hvg
    hvg_ord                      = rep(NA,nrow(sce))
    hvg_ord[hvg_bool]            = res$hvg
    rowData(sce)$cxds_hvg_bool   = hvg_bool
    rowData(sce)$cxds_hvg_ordr   = hvg_ord
    metadata(sce)$cxds$S         = res$S
    metadata(sce)$cxds$topPairs  = res$topPairs
    metadata(sce)$cxds$binThresh = res$binThresh
    #sce$cxds_score               = res$scores
  } else{
    if(verb) message("-> done.\n\n")
    #sce$cxds_score = -scores;
  }
  return(sce)
}

#' Extract top-scoring gene pairs from an SingleCellExperiment where cxds has been run
#'
#' @param sce single cell experiment to analyze; needs "counts" in assays slot.
#' @param n integer. The number of gene pairs to extract. Default: 100
#' @importFrom Matrix t
#' @return matrix Matrix with two colulmns, each containing gene indexes for gene pairs (rows).
cxds_getTopPairs <- function(sce,n=100){
  #=======================================
  
  ind = rowData(sce)$cxds_hvg_ordr[!is.na(rowData(sce)$cxds_hvg_ordr)]
  Bp  = counts(sce)[ind,] > metadata(sce)$cxds$binThresh
  imp = Matrix::t(Matrix::t(Bp) * sce$cxds_score)
  imp = imp  %*% Matrix::t(Bp)
  imp = imp * metadata(sce)$cxds$S
  imp = as.matrix(imp)
  res = list(imp=imp)
  
  colnames(imp) = seq_len(ncol(imp))
  rownames(imp) = colnames(imp)
  topPrs = matrix(NA,nrow=n,ncol=2)
  for(i in seq_len(n)){
    pr = which(imp == max(imp) ,arr.ind=TRUE)[c(1,2)]
    topPrs[i,] = as.integer(colnames(imp)[pr])
    imp = imp[-pr,]
    imp = imp[,-pr]
  }
  res$topPairs = topPrs
}

cxds_bcds_hybrid <- function(sce, cxdsArgs=NULL, bcdsArgs=NULL, verb=FALSE, estNdbl=FALSE,force=FALSE){
  #==============================================================================
  
  #- enforce consistency wrt doublet calls
  cxdsArgs$estNdbl = estNdbl
  bcdsArgs$estNdbl = estNdbl
  
  #- re-run only if scores are missing or estNdbl has not been performed.
  #  FIXME: cxdsArgs and bcdsArgs might be inconsistent between what was passed to this function and what was used in the \code{sce} that is passed.
  
  run_cxds  = is.null(colData(sce)$cxds_score)
  run_bcds  = is.null(colData(sce)$bcds_score)
  run_cxds  = run_cxds | (estNdbl & (is.null(metadata(sce)$cxds$ndbl)))
  run_bcds  = run_bcds | (estNdbl & (is.null(metadata(sce)$bcds$ndbl)))
  
  if(force) {
    run_bcds = TRUE
    run_cxds = TRUE
  }
  
  if(run_cxds) {
    if(verb) message("-> running cxds\n")
    sce = do.call(cxds,c(list(sce=sce),cxdsArgs))
  }
  if(run_bcds) {
    if(verb) message("-> running bcds\n")
    sce = do.call(bcds,c(list(sce=sce),bcdsArgs))
  }
  #- sqish and average scores
  squish  = function(x) { x = x-min(x) ; x = x/max(x) ; return(x) }
  
  s.cxds   = sce$cxds_score %>% squish
  s.bcds   = sce$bcds_score %>% squish
  s.hybrid = s.cxds + s.bcds
  #- put in sce
  sce$hybrid_score = s.hybrid
  
  if(estNdbl){
    #- hybrid score of artificial doublets
    s.cxds.sim = metadata(sce)$cxds$sim_scores
    s.bcds.sim = metadata(sce)$bcds$sim_scores
    if(length(s.cxds.sim) > length(s.bcds.sim)){
      s.cxds.sim = s.cxds.sim[sample(seq_len(length(s.bcds.sim)))]
    }
    if(length(s.cxds.sim) < length(s.bcds.sim)){
      s.bcds.sim = s.bcds.sim[sample(seq_len(length(s.cxds.sim)))]
    }
    s.hybrid.sim = (s.cxds.sim %>% squish) + (s.bcds.sim %>% squish)
    
    #- annotate
    est_dbl = get_dblCalls_ALL(sce$hybrid_score,s.hybrid.sim)
    if(is.null(metadata(sce)$hybrid)) metadata(sce)$hybrid = list()
    metadata(sce)$hybrid$ndbl = est_dbl
    metadata(sce)$hybrid$sim_scores = s.hybrid.sim
    sce$hybrid_call = sce$hybrid_score >= est_dbl["balanced","threshold"]
  }
  
  return(sce)
}

