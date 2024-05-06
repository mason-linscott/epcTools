#' EPC likelihood function
#'
#' Performs repeated maximum likelihood searches under an EPC process for one or more OU model parameters. Then, performs a dentist run to determine 95% CI for model parameters.
#'
#' @param cache A cache object created using epcTools.
#'
#' @param p Numeric vector of starting parameters to be loaded corresponding to parameters specified in cache object (format: theta, sig2, alpha; if any parameter varies with the environment the model parameters follow the same order). If NULL start.searcher function will be run instead.
#'
#' @param c_t A custom theta function under an EPC process, NULL if linear or step is being used.
#'
#' @param c_a A custom alpha function under an EPC process, NULL if linear or step is being used.
#'
#' @return A list including a table of parameters (theta,sig2, and alpha or any combination of EPC model parameters) with likelihood, execution time, model convergence parameters from optimx::optimx (ignore), and AIC. Dentist output with 95% confidence interval values for all parameters.
#'
#' @export
epc.max.lik<-function(cache,p,c_t=NULL,c_a=NULL,skip_dentist=FALSE){

  options(warn=-1)

  #set the max lik search function based on cache
  if (any(c("theta","alpha") %in%  cache$epc_params)){
    max.lik<-epc.lik
  }

  if (any(c("BM") %in%  cache$epc_params)){
    max.lik<-epcTools::bm.lik
  }

  if (any(c("OU") %in%  cache$epc_params)){
    max.lik<-epcTools::ou.lik
  }

  #run start searcher and any user parameters to get initial likelihoods
  if (missing(p) && is.null(c_t) && is.null(c_a)) {
    p<-as.numeric(start.searcher(cache))
    init_lik<-max.lik(cache,p,method="Nelder-Mead",c_a=NULL,c_t=NULL)
  }

  else {

    init_lik<-max.lik(cache,p,method="Nelder-Mead",c_a=NULL,c_t=NULL)

  }

  lik2_ref<-max.lik(cache,p=as.numeric(init_lik[1:length(p)]),method="BFGS",c_a=NULL,c_t=NULL)

  lik3_start<-max.lik(cache,method="Nelder-Mead",c_a=NULL,c_t=NULL)

  lik4_ref<-max.lik(cache,p=as.numeric(lik3_start[1:length(p)]),method="BFGS",c_a=NULL,c_t=NULL)

  lik5_start<-max.lik(cache,method="Nelder-Mead",c_a=NULL,c_t=NULL)

  lik6_ref<-max.lik(cache,p=as.numeric(lik5_start[1:length(p)]),method="BFGS",c_a=NULL,c_t=NULL)

  lik7_start<-max.lik(cache,method="Nelder-Mead",c_a=NULL,c_t=NULL)

  lik8_ref<-max.lik(cache,p=as.numeric(lik7_start[1:length(p)]),method="BFGS",c_a=NULL,c_t=NULL)

  final_table<-rbind(init_lik,lik2_ref,lik3_start,lik4_ref,lik5_start,lik6_ref,lik7_start,lik8_ref)

  final_table <- final_table[order(final_table$likelihood),]

  if (skip_dentist==FALSE) {
    dent_est <- epc.dent(cache,pars=as.numeric(head(final_table,n=1)[1:length(p)]), as.numeric(final_table[1,length(p)+1]))

    colnames(dent_est$all_ranges)<-colnames(final_table[1:length(p)])

    colnames(dent_est$results)<-c("lik",colnames(final_table)[1:length(p)])

    result_list<-list()

    result_list$fit <- data.frame(EPC_model=cache$epc_params, Type=cache$m_type, Lik=final_table$likelihood[1], Npar=length(p), AIC=final_table$AIC[1])

    result_list$dentist_output<- dent_est

    result_list$search_results<-final_table[-((length(p)+2):(length(p)+7))]

    class(result_list)<-c("epc","list")

    return(result_list)

  }
  if (skip_dentist==TRUE){
    return(head(final_table,n=1))
  }
}

