#' EPC cache building function
#'
#' Simulates traits under an EPC process for one or more OU model parameters
#'
#' @param mtree A phylogenetic tree
#'
#' @param epc_params EPC params included in model. Values can be 'alpha' , 'theta', 'OU', or 'BM'. 'BM' and 'OU' are to simulate Brownian motion or a single peak Ornstein-Uhlenbeck model, respectively.
#'
#' @param m_type A string corresponding to the EPC environment relationship, can be 'linear', 'step', or 'custom'.
#'
#' @param environment_data A two column table where column 1 is time and column 2 is the environment. Time is ordered such that 0 corresponds to the present and larger values indicate further back in time.
#'
#' @param trait_data A vector of numeric trait data where each value has a corresponding species name
#'
#' @param n_slice Number of time slices across the tree, must equal the 'X' parameter.
#'
#' @param custom.model.t A user specified function that outputs a value of theta given some environmental value.
#'
#' @param custom_params_t Parameters to simulate for a custom EPC theta relationship given 'custom.model.t'.
#'
#' @param custom.model.a A user specified function that outputs a value of alpha given some environmental value.
#'
#' @param custom_params_a Parameters to simulate for a custom EPC alpha relationship given 'custom.model.a'.
#'
#' @param custom_param_names Names of custom parameters for custom EPC - alpha/theta models.
#'
#' @return A cache object that includes the information needed for a maximum likelihood search.
#'
#' @export

cache.epc<-function(mtree,epc_params,m_type,environment_data,trait_data,n_slice,env_treat="bin",
                  custom.model.t=NULL,custom_params_t=NULL,custom.model.a=NULL,
                  custom_params_a=NULL,custom_param_names=NULL)
{
  cache_list<-list() #list of list stored parameters and data; this is the output

  #reconcile environmental curve and tree by removing tips outside of the environmental curve
  env_range<-range(environment_data[,1])
  tree_range<-range(nodeHeights(mtree))
  heights<-nodeHeights(mtree)
  tip_labels = mtree$tip.label
  outside_env_tips = tip_labels[which(heights[1:length(tip_labels)] < env_range[1] | heights[1:length(tip_labels)] > env_range[2])]
  if (length(outside_env_tips)>0) {
    mtree = ape::drop.tip(mtree, outside_env_tips)
    cat(paste0("Dropped ", length(outside_env_tips)," tips outside of the environments range. If this is unexpected check your data."))
    tip_labels=mtree$tip.label
  }

  trait_data <- trait_data[rownames(trait_data) %in% tip_labels, , drop=FALSE ]

  switch(env_treat,
         "bin"={
           sliceEnv<-getBinEnv(environment_data,n_slice)
         },
         "smooth"={
           sliceEnv<-getSliceEnv(environment_data,n_slice)
         })

  max_time<-max(nodeHeights(mtree)) # get max time for slicing

  rmap <- makeEpochSimmap(mtree, base::seq(0,max_time, length.out=(n_slice+1))[2:n_slice]) #make the epoch simmap

    dat <- as.matrix(trait_data)
    cache <- .prepare.ou.univariate(mtree, dat, SE=0, pred=NULL) #creates all the paremters in one nice cache to be called later

    ## Functions that converts simmap to singleton tree.
    siTree <- phytools::map.to.singleton(rmap$smap)
    siTree$edge.regime <- names(siTree$edge.length)
    siTree$edge.jump <- rep(0L, length(siTree$edge.length))
    cache$PCMtree <- PCMBase::PCMTree(siTree)

    ## Steps: Make model
    pcmModel <- PCMBase::MixedGaussian(1, modelTypes=rep("OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x", (n_slice)),
                                       mapping=1:(n_slice))

    # PCMBase::PCMParamLoadOrStore(pcmModel, pars2pcmbase(pars), offset=0L,
    #                              k=1, R = (n_slice), load = TRUE)

    pcmModel_l <- PCMBase::MixedGaussian(1, modelTypes=rep("OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x", (n_slice)),
                                         mapping=1:(n_slice))


    # PCMBase::PCMParamLoadOrStore(pcmModel_l, pars2pcmbase(pars), offset=0L,
    #                              k=1, R = (n_slice), load = TRUE)


    cache$pcmModel <- pcmModel #the actual model object
    cache$pcmModel_l <- pcmModel #the actual model object
    cache$PCM.X <- t(cache$dat) # the actual data object
    metaI <- PCMBase::PCMInfo(t(cache$dat), cache$PCMtree, cache$pcmModel) #meta information required by PCMBase
    cache$rmap<-rmap #retain simulated simmap in list
    cache$n_slice<-n_slice #retain number of slices for lik calc
    cache$environment<-sliceEnv
    cache$tree<-mtree
    cache$metaI<-metaI
    cache$m_type<-m_type
    cache$epc_params<-epc_params
    cache$custom_model_t<-custom.model.t
    cache$custom_params_t<-custom_params_t
    cache$custom_model_a<-custom.model.a
    cache$custom_params_a<-custom_params_a
    cache_list<-cache

  cache_list
}

