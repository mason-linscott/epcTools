#' BM likelihood function
#'
#' Performs a single maximum likelihood search under an Brownian Motion process
#'
#' @param cache A cache object created using epcTools.
#'
#' @param p Numeric vector of starting parameters to be loaded corresponding to theta and sig2. If NULL start.searcher function will be run instead.
#'
#' @param method A string corresponding to the optimx::optimx algorithm to be run, recommended 'Nelder-Mead' for first search followed by 'BFGS'.
#'
#' @param c_t A custom theta function under an EPC process, NULL for OU models.
#'
#' @param c_a A custom alpha function under an EPC process, NULL for OU models.
#'
#' @return A table of parameters (theta and sig2) with a likelihood, execution time, model convergence parameters from optimx::optimx (ignore), and AIC.
#'
#' @export
bm.lik <- function(cache,p,method,c_t=NULL,c_a=NULL) {
  if (missing(p) && is.null(c_t) && is.null(c_a)) {
    p=as.numeric(start.searcher(cache))
  }
  pcmModel_bm <- PCMBase::PCM(1, model="BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
  PCMtree<-PCMBase::PCMTree(cache$tree)
  metaI <- PCMBase::PCMInfo(cache$PCM.X, PCMtree, pcmModel_bm)

  BM.lik.pars2pcmbase <- function(pars,cache){
    p <- c(pars$theta, pars$sig2)
    return(p)
  }
  lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, PCMtree, model=pcmModel_bm,metaI=metaI)

  lik_func<-function(p,cache,c_t=NULL,c_a=NULL){
    pars=list()
    pars$theta <- p[1]
    pars$sig2 <- p[2]

    if(any(pars$sig2<0)) {
      lik<- 99999999999
      return(lik) }
    else{
      lik<-lik_funk(BM.lik.pars2pcmbase(pars,cache))[1]*-1
      return(lik)

    }
  }

  result<-optimx::optimx(c(p[1],p[2]),lik_func,cache=cache,method=method)
  colnames(result)[1:3]<-c("theta","sig2","likelihood")
  result$AIC<-(2*result$likelihood+2*(length(p)))
  return(result)
}

