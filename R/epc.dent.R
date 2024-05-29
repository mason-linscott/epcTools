#' EPC dentist search function
#'
#' Performs dentist search for all model parameters for BM, OU, and EPC models
#'
#' @param cache A cache object created using epcTools.
#'
#' @param pars Numeric vector of starting parameters to be loaded corresponding to parameters specified in cache object (format: theta, sig2, alpha; if any parameter varies with the environment the model parameters follow the same order). If NULL start.searcher function will be run instead.
#'
#' @param best_log Best likelihood value for dentist to search around
#'
#' @param c_t A custom theta function under an EPC process, NULL if linear or step is being used.
#'
#' @param c_a A custom alpha function under an EPC process, NULL if linear or step is being used.
#'
#' @return  Dentist output with 95% confidence interval values for all parameters.
#'
#' Note: This function can run into recursion limits on individual systems which may ultimately prevent it's use by epc.max.lik. Use the skip.dentist=TRUE flag on epc.max.lik if you run into this issue.
#'
#' @export
epc.dent<- function(cache,pars,best_log,c_t=NULL,c_a=NULL,reps=NULL){
  #lik_funk is the actual likelihood function - params go in, likelihood comes out.
  #all the rest is for getting parameters ready for optimx::optimx, or for finding starting parameters
  lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, cache$PCMtree, model=cache$pcmModel_l,
                                           metaI=cache$metaI)

  #lik_func is for feeding in parameters to optimx::optimx depending on model type
    if (any(c("BM","OU") %in% cache$epc_params)) {
      switch(cache$epc_params,
             "BM"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pcmModel_bm <- PCMBase::PCM(1, model="BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
               PCMtree<-PCMBase::PCMTree(cache$tree)
               metaI <- PCMBase::PCMInfo(cache$PCM.X, PCMtree, pcmModel_bm)
               BM.lik.pars2pcmbase <- function(pars_bm,cache){
                 p <- c(pars_bm$theta, pars_bm$sig2)
                 return(p)
               }
               lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, PCMtree, model=pcmModel_bm,metaI=metaI)
               pars_epc$theta <- pars[1]
               pars_epc$sig2 <- pars[2]
               if(any(pars_epc$sig2<0)) {
                 lik<- 99999999999
                 return(lik) }
               else{
                 lik<-lik_funk(BM.lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 return(lik)
               }
             }
             },

             "OU"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pcmModel_ou <- PCMBase::PCM(1,
model="OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
               PCMtree<-PCMBase::PCMTree(cache$tree)
               metaI <- PCMBase::PCMInfo(cache$PCM.X, PCMtree, pcmModel_ou)
               OU.lik.pars2pcmbase <- function(pars_ou,cache){
                 .p <- c(pars_ou$alpha, pars_ou$theta, pars_ou$sig2)
                 p <- c(.p, 0)
                 return(p)
               }

               lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, PCMtree, model=pcmModel_ou,metaI=metaI)
               pars_epc$theta <- pars[1]
               pars_epc$sig2 <- pars[2]
               pars_epc$alpha <- pars[3]
               if(any(pars_epc$alpha<0) || pars_epc$sig2<0) {
                 lik<- 99999999999
                 return(lik) }
               else{
                 lik<-lik_funk(OU.lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 return(lik)

               }
             }
            })
    }

    if (all(c("theta","alpha") %in%  cache$epc_params)) {

      switch(cache$m_type,

             "linear"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$b0_t <- pars[1]
               pars_epc$b1_t <- pars[2]
               pars_epc$sig2 <- pars[3]
               pars_epc$b0_a <- pars[4]
               pars_epc$b1_a <- pars[5]
               pars_epc$alpha <- pars_epc$b0_a + pars_epc$b1_a*cache$environment
               pars_epc$theta <- pars_epc$b0_t + pars_epc$b1_t*cache$environment
               if(any(pars_epc$alpha < 0) || pars_epc$sig2 < 0) {
                 lik<- 99999999999
                 return(lik)
               } else {
                 lik<-lik_funk(lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 pars_epc$lik <- lik
                 return(lik)
               }
             }},


             "step"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$step_thresh_t <- pars[1]
               pars_epc$step_high_t <- pars[2]
               pars_epc$step_low_t <- pars[3]
               pars_epc$sig2 <- pars[4]
               pars_epc$step_thresh_a <- pars[5]
               pars_epc$step_high_a <- pars[6]
               pars_epc$step_low_a <- pars[7]
               pars_epc$theta <- ifelse(cache$environment>pars_epc$step_thresh_t,pars_epc$step_high_t,pars_epc$step_low_t)
               pars_epc$alpha <- ifelse(cache$environment>pars_epc$step_thresh_a,pars_epc$step_high_a,pars_epc$step_low_a)
               if (any(pars_epc$step_stresh_t > max(cache$environment)) || any(pars_epc$step_stresh_t < min(cache$environment)) ||
                   any(pars_epc$step_stresh_a > max(cache$environment)) || any(pars_epc$step_stresh_a < min(cache$environment))) {
                 lik<- 99999999999
                 return(lik)
               }
               }},

             "custom"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$sig2<-pars[1]
               pars_epc$theta<-cache$custom_model_t(cache$environment,pars[2:(1+c_t)])
               pars_epc$alpha<-cache$custom_model_a(cache$environment,pars[(2+c_t):length(p)])
               if(any(pars_epc$alpha < 0) || pars_epc$sig2 < 0) {
                 lik<- 99999999999
                 return(lik)
               } else {
                 lik<-lik_funk(lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 pars_epc$lik <- lik
                 return(lik)
               }
               }})
    }

    if (!"theta" %in% cache$epc_params && "alpha" %in% cache$epc_params) {

      switch(cache$m_type,

             "linear"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$theta <- pars[1]
               pars_epc$sig2 <- pars[2]
               pars_epc$b0_a <- pars[3]
               pars_epc$b1_a <- pars[4]
               pars_epc$alpha <- pars_epc$b0_a + pars_epc$b1_a*cache$environment
               if(any(pars_epc$alpha < 0) || pars_epc$sig2 < 0) {
                 lik<- 99999999999
                 return(lik)
               } else {
                 lik<-lik_funk(lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 pars_epc$lik <- lik
                 return(lik)
               }
             }},

             "step"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$theta <- pars[1]
               pars_epc$sig2 <- pars[2]
               pars_epc$step_thresh_a <- pars[3]
               pars_epc$step_high_a <- pars[4]
               pars_epc$step_low_a <- pars[5]
               pars_epc$alpha <- ifelse(cache$environment>pars_epc$step_thresh_a,pars_epc$step_high_a,pars_epc$step_low_a)
               if (any(pars_epc$step_stresh_t > max(cache$environment)) || any(pars_epc$step_stresh_t < min(cache$environment)) ||
                   any(pars_epc$step_stresh_a > max(cache$environment)) || any(pars_epc$step_stresh_a < min(cache$environment))) {
                 lik<- 99999999999
                 return(lik)
               }
               }},

             "custom"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$theta<-pars[1]
               pars_epc$sig2<-pars[2]
               pars_epc$alpha<-cache$custom_model_a(cache$environment,pars[3:(length(p))])
               if(any(pars_epc$alpha < 0) || pars_epc$sig2 < 0) {
                 lik<- 99999999999
                 return(lik)
               } else {
                 lik<-lik_funk(lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 pars_epc$lik <- lik
                 return(lik)
               }
               }
             })

    }

    if ("theta" %in% cache$epc_params && !"alpha" %in% cache$epc_params) {

      switch(cache$m_type,

             "linear"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$b0_t <- pars[1]
               pars_epc$b1_t <- pars[2]
               pars_epc$sig2 <- pars[3]
               pars_epc$alpha <- pars[4]
               pars_epc$theta <- as.numeric(pars_epc$b0_t) + as.numeric(pars_epc$b1_t)*as.numeric(cache$environment)
               if(any(pars_epc$alpha < 0) || pars_epc$sig2 < 0) {
                 lik<- 99999999999
                 return(lik)
               } else {
                 lik<-lik_funk(lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 pars_epc$lik <- lik
                 return(lik)
               }

               }
             },

             "step"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$step_thresh_t <- pars[1]
               pars_epc$step_high_t <- pars[2]
               pars_epc$step_low_t <- pars[3]
               pars_epc$sig2 <- pars[4]
               pars_epc$alpha <- pars[5]
               pars_epc$theta <- ifelse(cache$environment>pars_epc$step_thresh_t,pars_epc$step_high_t,pars_epc$step_low_t)
               if (any(pars_epc$step_stresh_t > max(cache$environment)) || any(pars_epc$step_stresh_t < min(cache$environment)) ||
                   any(pars_epc$step_stresh_a > max(cache$environment)) || any(pars_epc$step_stresh_a < min(cache$environment))) {
                 lik<- 99999999999
                 return(lik)
               }

             }},
             "custom"={
               lik_func<-function(pars,cache,c_t,c_a){
                 pars_epc=list()
               pars_epc$sig2=pars[1]
               pars_epc$alpha=pars[2]
               pars_epc$theta<-cache$custom_model_t(cache$environment,pars[3:(length(p))])
               if(any(pars_epc$alpha < 0) || pars_epc$sig2 < 0) {
                 lik<- 99999999999
                 return(lik)
               } else {
                 lik<-lik_funk(lik.pars2pcmbase(pars_epc,cache))[1]*-1
                 pars_epc$lik <- lik
                 return(lik)
               }
               }})
    }

  nsteps_epc<-length(pars)*500
  dent_result<-dentist::dent_walk(par=pars, fn=lik_func, best_neglnL=best_log, nsteps=nsteps_epc, print_freq=500, cache=cache)

  return(dent_result)

}

