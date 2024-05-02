#' EPC simulation function
#'
#' Simulates traits under an EPC process for one or more OU model parameters
#'
#' @param mtree A phylogenetic tree
#'
#' @param epc_params EPC params included in model. Values can be 'alpha' , 'theta', 'OU', or 'BM'. 'BM' and 'OU' are to simulate non EPC models Brownian motion and Ornstein-Uhlenbeck, respectively.
#'
#' @param m_type A string corresponding to the EPC environment relationship, can be 'linear', 'step', or 'custom'.
#'
#' @param X A vector corresponding to environmental values at each time slice
#'
#' @param n_slice Number of time slices across the tree, must equal the 'X' parameter.
#'
#' @param reps_per_tree Number of times to simulate data under the given parameters.
#'
#' @param theta A numeric value for theta for simulation.
#'
#' @param sig2 A numeric value for sig2 for simulation.
#'
#' @param alpha A numeric value for alpha for simulation.
#'
#' @param b0_a A numeric value for the intercept b0 under a linear relationship between the environment and alpha.
#'
#' @param b1_a A numeric value for the slope b1 under a linear relationship between the environment and alpha.
#'
#' @param b0_t A numeric value for the intercept b0 under a linear relationship between the environment and theta.
#'
#' @param b1_t A numeric value for the slope b1 under a linear relationship between the environment and theta.
#'
#' @param step_thresh_t A numeric value for the threshold of a step relationship between the environment and theta.
#'
#' @param step_high_t A numeric value of theta when the environment exceeds step_thresh_t'.
#'
#' @param step_low_t A numeric value of theta when the environment falls below step_thresh_t'.
#'
#' @param step_thresh_a A numeric value for the threshold of a step relationship between the environment and alpha.
#'
#' @param step_high_a A numeric value of alpha when the environment exceeds step_thresh_a'.
#'
#' @param step_low_a A numeric value of alpha when the environment falls below step_thresh_a'.
#'
#' @param custom.model.t A user specified function that outputs a value of theta given some environmental value.
#'
#' @param custom_params_t Parameters to simulate for a custom EPC theta relationship given 'custom.model.t'.
#'
#' @param custom.model.a A user specified function that outputs a value of alpha given some environmental value.
#'
#' @param custom_params_a Parameters to simulate for a custom EPC alpha relationship given 'custom.model.at'.
#'
#' @param custom_param_names Names of custom parameters for custom EPC - alpha/theta models.
#'
#' @return A cache object that includes the information needed for a maximum likelihood search.
#'
#' @export

sim.EPC<-function(mtree,epc_params,m_type,X,n_slice,reps_per_tree=1,
                  theta=NULL,sig2,alpha=NULL,
                  b0_a=NULL,b1_a=NULL,b0_t=NULL,b1_t=NULL,
                  step_thresh_t=NULL,step_high_t=NULL,step_low_t=NULL,
                  step_thresh_a=NULL,step_high_a=NULL,step_low_a=NULL,
                  custom.model.t=NULL,custom_params_t=NULL,custom.model.a=NULL,
                  custom_params_a=NULL,custom_param_names=NULL)
{
  if("theta" %in%  epc_params && is.null(step_thresh_t) && is.null(custom.model.t) && is.null(b0_t)) {
    stop("Missing EPC model/model parameters for theta")
  }

  if("alpha" %in%  epc_params && is.null(step_thresh_a) && is.null(custom.model.a) && is.null(b0_a)) {
    stop("Missing EPC model/model parameters for alpha")
  }

  cache_list<-list() #list of list stored parameters and data; this is the output
  #resolution<-3*10p[[1*(n_slice)
  #sliceEnv<-getSliceEnv(X,n_slice,resolution)[,1]
  #sliceEnv[11]<-sliceEnv[10]
  sliceEnv<-X # specified time slices
  max_time<-max(nodeHeights(mtree)) # get max time for slicing

  for(j in 1:reps_per_tree){
    rmap <- makeEpochSimmap(mtree, seq(0,max_time, length.out=(n_slice+1))[2:n_slice]) #make the epoch simmap

    ##pars## a list of parameters that will be used to make the EPC models
    pars <- list()
    pars$sb <- rmap$pars$sb
    pars$loc <- rmap$pars$loc
    pars$t2 <- rmap$pars$t2
    pars$k <- rmap$pars$k
    pars$ntheta <- n_slice

    #Parameters for EPC model
    pars$theta <- theta
    pars$alpha <- alpha
    pars$sig2 <- sig2
    pars$b0_a <- b0_a
    pars$b1_a <- b1_a
    pars$b0_t <- b0_t
    pars$b1_t <- b1_t
    pars$step_thresh_a <- step_thresh_a
    pars$step_high_a <- step_high_a
    pars$step_low_a <- step_low_a
    pars$step_thresh_t <- step_thresh_t
    pars$step_high_t <- step_high_t
    pars$step_low_t <- step_low_t

    #EPC_params
    if ("theta" %in%  epc_params) {
      if (m_type=="linear") {
        pars$theta <- pars$b0_t + pars$b1_t * sliceEnv
      }

      if (m_type=="step") {
        pars$theta<-ifelse(sliceEnv>pars$step_thresh_t,pars$step_high_t,pars$step_low_t)
      }

      if (m_type=="custom") {
        pars$theta<-custom.model.t(sliceEnv[1],custom_params_t)
        for (i in 1:length(custom_params_t)) {
          pars[[paste0("c_t", i)]]<-custom_params_t[i]
        }
      }
    }

    if ("alpha" %in%  epc_params) {
      if (m_type=="linear") {
        pars$alpha <- pars$b0_a + pars$b1_a * sliceEnv
      }

      if (m_type=="step") {
        pars$alpha<-ifelse(sliceEnv>pars$step_thresh_a,pars$step_high_a,pars$step_low_a)
      }

      if (m_type=="custom") {
        pars$alpha<-custom.model.a(sliceEnv,custom_params_a)
        for (i in 1:length(custom_params_a)) {
          pars[[paste0("c_a", i)]]<-custom_params_a[i]
        }
      }

    }


    dat <- setNames(rnorm(Ntip(mtree), 0, 1), mtree$tip.label) #dummy data
    cache <- bayou:::.prepare.ou.univariate(mtree, dat, SE=0, pred=NULL) #creates all the paremters in one nice cache to be called later

    ## Functions that converts simmap to singleton tree.
    siTree <- phytools::map.to.singleton(rmap$smap)
    siTree$edge.regime <- names(siTree$edge.length)
    siTree$edge.jump <- rep(0L, length(siTree$edge.length))
    cache$PCMtree <- PCMTree(siTree)

    ## Steps: Make model
    pcmModel <- MixedGaussian(1, modelTypes=rep("OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x", (n_slice)),
                              mapping=1:(n_slice))

    pars2pcmbase <- function(pars){
      if (all(c("theta","alpha") %in%  epc_params)) {
        .p <- lapply(1:(n_slice), function(x) c(pars$alpha[x], pars$theta[x], pars$sig2))
        p <- c(pars$theta[1], do.call(c, .p), 0)
        return(p)
      }

      if (!"theta" %in% epc_params && "alpha" %in% epc_params) {
        .p <- lapply(1:(n_slice), function(x) c(pars$alpha[x], pars$theta, pars$sig2))
        p <- c(pars$theta[1], do.call(c, .p), 0)
        return(p)
      }

      if ("theta" %in% epc_params && !"alpha" %in% epc_params) {
        .p <- lapply(1:(n_slice), function(x) c(pars$alpha, pars$theta[x], pars$sig2))
        p <- c(pars$theta[1], do.call(c, .p), 0)
        return(p)
      }
    }


    PCMBase::PCMParamLoadOrStore(pcmModel, pars2pcmbase(pars), offset=0L,
                                 k=1, R = (n_slice), load = TRUE)

    pcmModel_l <- MixedGaussian(1, modelTypes=rep("OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x", (n_slice)),
                                mapping=1:(n_slice))


    PCMBase::PCMParamLoadOrStore(pcmModel_l, pars2pcmbase(pars), offset=0L,
                                 k=1, R = (n_slice), load = TRUE)


    cache$pcmModel <- pcmModel #the actual model object
    cache$pcmModel_l <- pcmModel #the actual model object
    cache$PCM.X <- t(cache$dat) # the actual data object
    metaI <- PCMInfo(t(cache$dat), cache$PCMtree, cache$pcmModel) #meta information required by PCMBase
    simdat <- PCMSim(cache$PCMtree, model=cache$pcmModel, X0=cache$pcmModel$X0[1]) #here is the simulated data
    cache$PCM.X <- t(simdat[1,1:Ntip(mtree)]) #retain simulated data in list
    cache$metaI <- PCMInfo(t(simdat[1,1:Ntip(mtree)]), cache$PCMtree, cache$pcmModel) #retain simulated data in list
    cache$rmap<-rmap #retain simulated simmap in list
    cache$n_slice<-n_slice #retain number of slices for lik calc
    cache$environment<-sliceEnv
    cache$tree<-mtree
    cache$m_type<-m_type
    cache$pars<-pars
    cache$epc_params<-epc_params
    cache$custom_model_t<-custom.model.t
    cache$custom_params_t<-custom_params_t
    cache$custom_model_a<-custom.model.a
    cache$custom_params_a<-custom_params_a
    cache_list[[j]]<-cache
  }
  cache_list
}

