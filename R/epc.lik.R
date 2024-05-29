#' EPC likelihood function
#'
#' Performs a single maximum likelihood search under an EPC process for one or more OU model parameters
#'
#' @param cache A cache object created using epcTools.
#'
#' @param p Numeric vector of starting parameters to be loaded corresponding to parameters specified in cache object (format: theta, sig2, alpha; if any parameter varies with the environment the model parameters follow the same order). If NULL start.searcher function will be run instead.
#'
#' @param method A string corresponding to the optimx algorithm to be run, recommended 'Nelder-Mead' for first search followed by 'BFGS'.
#'
#' @param c_t A custom theta function under an EPC process, NULL if linear or step is being used.
#'
#' @param c_a A custom alpha function under an EPC process, NULL if linear or step is being used.
#'
#' @return A table of parameters (theta,sig2, and alpha or any combination of EPC model parameters) with a likelihood, execution time, model convergence parameters from optimx::optimx (ignore), and AIC.
#'
#' @export
epc.lik<- function(cache,p,method,c_t=NULL,c_a=NULL){
  #lik_funk is the actual likelihood function - params go in, likelihood comes out.
  #all the rest is for getting parameters ready for optimx::optimx, or for finding starting parameters
  lik_funk <- PCMBase::PCMCreateLikelihood(cache$PCM.X, cache$PCMtree, model=cache$pcmModel_l,
                                           metaI=cache$metaI)
  if (missing(p) && is.null(c_t) && is.null(c_a)) {
    p=as.numeric(start.searcher(cache))
  }

  #lik_func is for feeding in parameters to optimx::optimx depending on model type
  lik_func<-function(p,cache,c_t,c_a){
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

    if (any(pars$step_stresh_t > max(cache$environment)) || any(pars$step_stresh_t < min(cache$environment)) ||
        any(pars$step_stresh_a > max(cache$environment)) || any(pars$step_stresh_a < min(cache$environment))) {
      lik<- 99999999999
      return(lik)
    }
    #print(pars$step_thresh)

    if(any(pars$alpha < 0) || pars$sig2 < 0) {
      lik<- 99999999999
      return(lik)
    } else {
      lik<-lik_funk(lik.pars2pcmbase(pars,cache))[1]*-1
      pars$lik <- lik
      #print(pars$step_thresh)
      #print(pars$lik)
      return(lik)
    }
  }

  if (all(c("theta","alpha") %in%  cache$epc_params)) {
    if (cache$m_type=="linear") {
      result<-optimx::optimx(c(p[1],p[2],p[3],p[4],p[5]),lik_func,cache=cache,method=method)
      colnames(result)[1:6]<-c("b0_t","b1_t","sig2","b0_a","b1_a","likelihood")
    }

    if (cache$m_type=="step") {
      result<-optimx::optimx(c(p[1],p[2],p[3],p[4],p[5],p[6],p[7]),lik_func,cache=cache,method=method)
      colnames(result)[1:8]<-c("step_thresh_t","step_high_t","step_low_t","sig2","step_thresh_a","step_high_a","step_low_a","likelihood")
    }

    if (cache$m_type=="custom") {
      result<-optimx::optimx(p,lik_func,c_t=c_t,c_a=c_a,cache=cache,method=method)
      #colnames(result)[1:6]<-c("step_thresh_t","step_high_t","step_low_t","sig2","alpha","likelihood")
    }

  }
  if (!"theta" %in% cache$epc_params && "alpha" %in% cache$epc_params) {
    if (cache$m_type=="linear") {
      result<-optimx::optimx(c(p[1],p[2],p[3],p[4]),lik_func,cache=cache,method=method)
      colnames(result)[1:5]<-c("theta","sig2","b0_a","b1_a","likelihood")
    }

    if (cache$m_type=="step") {
      result<-optimx::optimx(c(p[1],p[2],p[3],p[4],p[5]),lik_func,cache=cache,method=method)
      colnames(result)[1:6]<-c("theta","sig2","step_thresh_a","step_high_a","step_low_a","likelihood")
    }

    if (cache$m_type=="custom") {
      result<-optimx::optimx(p,lik_func,c_a=c_a,cache=cache,method=method)
      #colnames(result)[1:6]<-c("step_thresh_t","step_high_t","step_low_t","sig2","alpha","likelihood")
    }


  }
  if ("theta" %in% cache$epc_params && !"alpha" %in% cache$epc_params) {
    if (cache$m_type=="linear") {
      result<-optimx::optimx(c(p[1],p[2],p[3],p[4]),lik_func,cache=cache,method=method)
      colnames(result)[1:5]<-c("b0_t","b1_t","sig2","alpha","likelihood")
    }

    if (cache$m_type=="step") {
      result<-optimx::optimx(c(p[1],p[2],p[3],p[4],p[5]),lik_func,cache=cache,method=method)
      colnames(result)[1:6]<-c("step_thresh_t","step_high_t","step_low_t","sig2","alpha","likelihood")
    }

    if (cache$m_type=="custom") {
      result<-optimx::optimx(p,lik_func,c_t=c_t,cache=cache,method=method)
      #colnames(result)[1:6]<-c("step_thresh_t","step_high_t","step_low_t","sig2","alpha","likelihood")
    }

  }
  result$AIC<-(2*result$likelihood+2*(length(p)))

  return(result)
}

