#' Basic likelihood function
#'
#' Returns a likelihood value given a cache object and vector of parameters
#'
#' @param cache A cache object created using epcTools.
#'
#' @param pars Numeric vector of parameters to be loaded.
#'
#' @return A numeric likelihood value
#'
#' @export
base.lik <- function(cache,p,c_t=NULL,c_a=NULL) {

  lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, cache$PCMtree, model=cache$pcmModel_l,metaI=cache$metaI)
  pars=list()
  if (all(c("theta","alpha") %in%  cache$epc_params)) {

    if (cache$m_type=="linear") {
      pars$b0_t <- p[1]
      pars$b1_t <- p[2]
      pars$sig2 <- p[3]
      pars$b0_a <- p[4]
      pars$b1_a <- p[5]
      pars$alpha <- pars$b0_a + pars$b1_a*cache$environment
      pars$theta <- pars$b0_t + pars$b1_t*cache$environment
    }

    if (cache$m_type=="step") {
      pars$step_thresh_t <- p[1]
      pars$step_high_t <- p[2]
      pars$step_low_t <- p[3]
      pars$sig2 <- p[4]
      pars$step_thresh_a <- p[5]
      pars$step_high_a <- p[6]
      pars$step_low_a <- p[7]
      pars$theta <- ifelse(cache$environment>pars$step_thresh_t,pars$step_high_t,pars$step_low_t)
      pars$alpha <- ifelse(cache$environment>pars$step_thresh_a,pars$step_high_a,pars$step_low_a)

    }

    if (cache$m_type=="custom") {
      pars$sig2<-p[1]
      pars$theta<-cache$custom_model_t(cache$environment,p[2:(1+c_t)])
      pars$alpha<-cache$custom_model_a(cache$environment,p[(2+c_t):length(p)])
    }

  }
  if (!"theta" %in% cache$epc_params && "alpha" %in% cache$epc_params) {

    if (cache$m_type=="linear") {
      pars$theta <- p[1]
      pars$sig2 <- p[2]
      pars$b0_a <- p[3]
      pars$b1_a <- p[4]
      pars$alpha <- pars$b0_a + pars$b1_a*cache$environment
    }

    if (cache$m_type=="step") {
      pars$theta <- p[1]
      pars$sig2 <- p[2]
      pars$step_thresh_a <- p[3]
      pars$step_high_a <- p[4]
      pars$step_low_a <- p[5]
      pars$alpha <- ifelse(cache$environment>pars$step_thresh_a,pars$step_high_a,pars$step_low_a)

    }

    if (cache$m_type=="custom") {
      pars$theta<-p[1]
      pars$sig2<-p[2]
      pars$alpha<-cache$custom_model_a(cache$environment,p[3:(length(p))])
    }

  }
  if ("theta" %in% cache$epc_params && !"alpha" %in% cache$epc_params) {

    if (cache$m_type=="linear") {
      pars$b0_t <- p[1]
      pars$b1_t <- p[2]
      pars$sig2 <- p[3]
      pars$alpha <- p[4]
      pars$theta <- pars$b0_t + pars$b1_t*cache$environment
    }

    if (cache$m_type=="step") {
      pars$step_thresh_t <- p[1]
      pars$step_high_t <- p[2]
      pars$step_low_t <- p[3]
      pars$sig2 <- p[4]
      pars$alpha <- p[5]
      pars$theta <- ifelse(cache$environment>pars$step_thresh_t,pars$step_high_t,pars$step_low_t)
    }

    if (cache$m_type=="custom") {
      pars$sig2=p[1]
      pars$alpha=p[2]
      pars$theta<-cache$custom_model_t(cache$environment,p[3:(length(p))])
    }

  }

  if ("OU" %in% cache$epc_params) {

    pcmModel_ou <- PCM(1, model="OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
    PCMtree<-PCMTree(cache$tree)
    metaI <- PCMInfo(cache$PCM.X, PCMtree, pcmModel_ou)

    ou.lik.pars2pcmbase <- function(pars,cache){
      .p <- c(pars$alpha, pars$theta, pars$sig2)
      p <- c(.p, 0)
      return(p)
    }
    lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, PCMtree, model=pcmModel_ou,metaI=metaI)

    pars$theta <- p[1]
    pars$sig2 <- p[2]
    pars$alpha <- p[3]

    if(pars$alpha<0 || pars$sig2<0) {
      lik<- 99999999999
      return(lik) }
    else{
      lik<-lik_funk(ou.lik.pars2pcmbase(pars,cache))[1]*-1
      return(lik)

    }
  }

  if ("BM" %in% cache$epc_params) {
    pcmModel_bm <- PCM(1, model="BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
    PCMtree<-PCMTree(cache$tree)
    metaI <- PCMInfo(cache$PCM.X, PCMtree, pcmModel_bm)

    BM.lik.pars2pcmbase <- function(pars,cache){
      p <- c(pars$theta, pars$sig2)
      return(p)
    }
    lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, PCMtree, model=pcmModel_bm,metaI=metaI)

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



  if(any(pars$alpha<0) || pars$sig2<0) {
    lik<- 99999999999
    return(lik) }
  else{
    lik<-lik_funk(lik.pars2pcmbase(pars,cache))[1]*-1
    return(lik)
  }
}

