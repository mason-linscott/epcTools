#' Consistency simulator of EPC results
#'
#' Given parameter results from an epc.max.lik, simulate many EPC cache objects with same parameters to determine EPC model identifiability
#'
#' @param params A single numeric row of parameter results from epc.max.lik,
#'
#' @param cache An EPC cache object used for epc.max.lik
#'
#' @param reps Number of simulations to perform
#'
#' @param cores Number of cores to use for simulations, do not exceed your computers limits!
#'
#' @return A list with maximum likelihood search results for each reptition under all model types and all results weighted by AIC weights.
#'
#' @export
consist.sim <- function(params,cache,reps,cores,c_t=NULL,c_a=NULL) {
  tryCatch({
    theta<-params$theta
    sig2<-params$sig2
    alpha<-params$alpha

    b0_t<-params$b0_t
    b1_t<-params$b1_t
    b0_a<-params$b0_a
    b1_a<-params$b1_a

    step_thresh_t<-params$step_thresh_t
    step_high_t<-params$step_high_t
    step_low_t<-params$step_low_t
    step_thresh_a<-params$step_thresh_a
    step_high_a<-params$step_high_a
    step_low_a<-params$step_low_a
    if (cache$m_type=="custom") {
      custom.model.a<-cache$custom.model.a
      custom.model.t<-cache$custom.model.t
      custom_params_t<-params[2:(2+c_t)]
      custom_params_a<-params[(2+c_t):(2+c_t+c_a)]
    }

    custom.model.a<-NULL
    custom.model.t<-NULL
    custom_params_t<-NULL
    custom_params_a<-NULL


  },
  error=function(e) {
    return(NULL)
  }
  )
  power_sim<-sim.EPC(mtree=cache$tree,cache$epc_params,cache$m_type,cache$environment,cache$n_slice,reps_per_tree=reps,
                     theta=theta,sig2=sig2,alpha=alpha,
                     b0_a=b0_a,b1_a=b1_a,b0_t=b0_t,b1_t=b1_t,
                     step_thresh_t=step_thresh_t,step_high_t=step_high_t,step_low_t=step_low_t,
                     step_thresh_a=step_thresh_a,step_high_a=step_high_a,step_low_a=step_low_a,
                     custom.model.t=custom.model.t,custom.model.a=custom.model.a,
                     custom_params_a=custom_params_a,custom_params_t=custom_params_t,custom_param_names=NULL)

  #Create consist starting parameters; ADD STEP LATER 04/23/24
  if(!is.null(alpha)) {
    alpha_s<-alpha
    b0a_s<-alpha
    b1a_s<-0}
  else {
    alpha_s<-mean(b0_a + b1_a*cache$environment)
    b0a_s<-b0_a
    b1a_s<-b1_a
  }

  if(!is.null(theta)) {
    theta_s<-theta
    b0t_s<-theta
    b1t_s<-0}
  else {
    theta_s<-mean(b0_a + b1_a*cache$environment)
    b0t_s<-b0_t
    b1t_s<-b1_t
  }

  sig_s <- sig2

  cl <- makeCluster(cores) #not to overload your computer
  registerDoParallel(cl)

  BM_consist <- foreach(f=1:reps, .combine=rbind,.packages=c('phytools','bayou','geiger','BAMMtools','dplyr','PCMBase','castor')) %dopar% {
    require(epcTools)
    power_sim[[f]]$epc_params<-"BM"
    tempMatrix = as.data.frame(epc.max.lik(power_sim[[f]],c(theta_s,sig_s),skip_dentist=TRUE)) #calling a function
    #do other things if you want
    tempMatrix #Equivalent to bm_1f_50 = rbind(bm_1f_50, tempMatrix)
  }
  print("BM OK")
  OU_consist <- foreach(f=1:reps, .combine=rbind,.packages=c('phytools','bayou','geiger','BAMMtools','dplyr','PCMBase','castor')) %dopar% {
    require(epcTools)
    power_sim[[f]]$epc_params<-"OU"
    tempMatrix = as.data.frame(epc.max.lik(power_sim[[f]],c(theta_s,sig_s,alpha_s),skip_dentist=TRUE)) #calling a function
    #do other things if you want
    tempMatrix #Equivalent to bm_1f_50 = rbind(bm_1f_50, tempMatrix)
  }
  print("OU OK")

  EPC_a_consist <- foreach(f=1:reps, .combine=rbind,.packages=c('phytools','bayou','geiger','BAMMtools','dplyr','PCMBase','castor')) %dopar% {
    require(epcTools)
    power_sim[[f]]$epc_params<-"alpha"
    tempMatrix = as.data.frame(epc.max.lik(power_sim[[f]],c(theta_s,sig_s,b0a_s,b1a_s),skip_dentist=TRUE)) #calling a function
    #do other things if you want
    tempMatrix #Equivalent to bm_1f_50 = rbind(bm_1f_50, tempMatrix)
  }
  print("EPC A OK")

  EPC_t_consist <- foreach(f=1:reps, .combine=rbind,.packages=c('phytools','bayou','geiger','BAMMtools','dplyr','PCMBase','castor')) %dopar% {
    require(epcTools)
    power_sim[[f]]$epc_params<-"theta"
    tempMatrix = as.data.frame(epc.max.lik(power_sim[[f]],c(b0t_s,b1t_s,sig_s,alpha_s),skip_dentist=TRUE)) #calling a function
    #do other things if you want
    tempMatrix #Equivalent to bm_1f_50 = rbind(bm_1f_50, tempMatrix)
  }
  print("EPC T OK")

  EPC_at_consist <- foreach(f=1:reps, .combine=rbind,.packages=c('phytools','bayou','geiger','BAMMtools','dplyr','PCMBase','castor')) %dopar% {
    require(epcTools)
    power_sim[[f]]$epc_params<-c("alpha","theta")
    tempMatrix = as.data.frame(epc.max.lik(power_sim[[f]],c(b0t_s,b1t_s,sig_s,b0a_s,b1a_s),skip_dentist=TRUE)) #calling a function
    #do other things if you want
    tempMatrix #Equivalent to bm_1f_50 = rbind(bm_1f_50, tempMatrix)
  }

  print("EPC AT OK")
  doParallel::stopCluster(cl)

  power_results<-list()
  power_results$BM<-BM_consist
  power_results$OU<-OU_consist
  power_results$EPC_a<-EPC_a_consist
  power_results$EPC_t<-EPC_t_consist
  power_results$EPC_at<-EPC_at_consist

  AIC_data<-data.frame(BM=BM_consist$AIC,OU=OU_consist$AIC,EPC_a=EPC_a_consist$AIC,EPC_t=EPC_t_consist$AIC,EPC_at=EPC_at_consist$AIC)

  power_results$AIC<-AIC_data

  power_results$AIC_weights<-as.data.frame(calculate_aic_weights(AIC_data))

  # Function to calculate AIC weights

  return(power_results)
}

