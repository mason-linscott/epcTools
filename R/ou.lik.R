#' OU likelihood function
#'
#' Performs a single maximum likelihood search under an Ornstein-Uhlenbeck process
#'
#' @param cache A cache object created using epcTools.
#'
#' @param p Numeric vector of starting parameters to be loaded corresponding to theta, sig2, and alpha. If NULL start.searcher function will be run instead.
#'
#' @param method A string corresponding to the optimx algorithm to be run, recommended 'Nelder-Mead' for first search followed by 'BFGS'.
#'
#' @param c_t A custom theta function under an EPC process, NULL for OU models.
#'
#' @param c_a A custom alpha function under an EPC process, NULL for OU models.
#'
#' @return A table of parameters (theta, sig2, and alpha) with a likelihood, execution time, model convergence parameters from optimx (ignore), and AIC.
#'
#' @export
ou.lik <- function(cache,p,method,c_t=NULL,c_a=NULL) {
  if (missing(p) && is.null(c_t) && is.null(c_a)) {
    p=as.numeric(start.searcher(cache))
  }
  pcmModel_ou <- PCM(1, model="OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
  PCMtree<-PCMTree(cache$tree)
  metaI <- PCMInfo(cache$PCM.X, PCMtree, pcmModel_ou)

  ou.lik.pars2pcmbase <- function(pars,cache){
    .p <- c(pars$alpha, pars$theta, pars$sig2)
    p <- c(.p, 0)
    return(p)
  }
  lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, PCMtree, model=pcmModel_ou,metaI=metaI)

  lik_func<-function(p,cache,c_t=NULL,c_a=NULL){
    pars=list()
    pars$theta <- p[1]
    pars$sig2 <- p[2]
    pars$alpha <- p[3]

    if(any(pars$alpha<0) || pars$sig2<0) {
      lik<- 99999999999
      return(lik) }
    else{
      lik<-lik_funk(ou.lik.pars2pcmbase(pars,cache))[1]*-1
      return(lik)

    }
  }

  result<-optimx::optimx(c(p[1],p[2],p[3]),lik_func,cache=cache,method=method)
  colnames(result)[1:4]<-c("theta","sig2","alpha","likelihood")
  result$AIC<-(2*result$likelihood+2*(length(p)))
  return(result)
}

