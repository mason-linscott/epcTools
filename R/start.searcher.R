#' Starting parameter searcher
#'
#' Finds starting parameters for BM, OU, and EPC models
#'
#' @param cache A cache object created using epcTools.
#'
#' @return  A vector of starting parameters corresponding to information in the cache object.
#'
#' @export
start.searcher<-function(cache) {

  #subfunctions for creating b0 and b1 for alpha and theta
  b0b1a_gen<-function(cache) {
    strong_b1_a_positive<-median(max_alpha/abs(cache$environment))

    strong_b0<-max_alpha
    weak_b0<-0

    b1_a_base_seq<-stats::rnorm(100,0,strong_b1_a_positive)*stats::runif(1,0.9,1)
    b0_a_base_seq<-stats::runif(100,0,strong_b0)*stats::runif(1,0.9,1)

    b0b1_pre<-expand.grid(b0_a_base_seq,b1_a_base_seq)
    b0b1_fin<-data.frame(b0a=numeric(),b1a=numeric())
    pre_alpha<-list()

    for(i in 1:nrow(b0b1_pre)) {
      range_accept<-c(0,max_alpha)
      pre_alpha[[i]]<-(b0b1_pre[i,1]+b0b1_pre[i,2]*cache$env)
      if (all(spatstat.utils::inside.range(pre_alpha[[i]],range_accept))==TRUE) {
        b0b1_fin[i,]<-b0b1_pre[i,]
      }
    }

    b0b1_fin<-stats::na.omit(b0b1_fin)
    rownames(b0b1_fin) <- NULL

    b0b1_unique <- b0b1_fin %>%
      group_by(b0a) %>%
      sample_n(1) %>%
      ungroup()

    b0b1_unique <-as.data.frame(b0b1_unique)
    return(b0b1_unique)
  }

  b0b1t_gen<-function(cache) {
    strong_b1_t_positive<-median(max(cache$PCM.X)/abs(cache$environment))
    strong_b0<-mean(cache$PCM.X)

    b1_t_base_seq<-stats::rnorm(100,0,strong_b1_t_positive)*stats::runif(1,0.9,1)
    b0_t_base_seq<-stats::runif(100,0,strong_b0)*stats::runif(1,0.9,1)

    b0b1_pre<-expand.grid(b0_t_base_seq,b1_t_base_seq)
    b0b1_fin<-data.frame(b0t=numeric(),b1t=numeric())
    pre_theta<-list()

    for(i in 1:nrow(b0b1_pre)) {
      range_accept<-c(abs(strong_b0)*-1,abs(strong_b0))
      pre_theta[[i]]<-(b0b1_pre[i,1]+b0b1_pre[i,2]*cache$env)
      if (all(spatstat.utils::inside.range(pre_theta[[i]],range_accept))==TRUE) {
        b0b1_fin[i,]<-b0b1_pre[i,]
      }
    }

    b0b1_fin<-stats::na.omit(b0b1_fin)
    rownames(b0b1_fin) <- NULL

    b0b1_unique <- b0b1_fin %>%
      group_by(b0t) %>%
      sample_n(1) %>%
      ungroup()

    b0b1_unique <-as.data.frame(b0b1_unique)
    return(b0b1_unique)
  }

  tip_root_length<-adephylo::distRoot(cache$phy)
  max_height<-max(nodeHeights(cache$phy))
  max_alpha<-log(2)/quantile(BAMMtools:::NU.branching.times(cache$phy),probs = base::seq(0, 1, 0.05))[4]

  if ("OU" %in% cache$epc_params)  {
    alpha_base_seq<-base::seq(0,max_alpha,max_alpha/10)[1:10]*stats::runif(1,0.9,1)
    theta_mean=mean(cache$PCM.X)
    theta_base_seq=base::seq(theta_mean,theta_mean*10,theta_mean)*stats::runif(1,0.9,1)
    sigma_mean=var(as.numeric(cache$PCM.X))/max(nodeHeights(cache$PCMtree))
    sigma_base_seq=c(sigma_mean*stats::runif(1,0.9,1),sigma_mean*5*stats::runif(1,0.9,1),sigma_mean*10*stats::runif(1,0.9,1),sigma_mean*50*stats::runif(1,0.9,1),sigma_mean*100*stats::runif(1,0.9,1),sigma_mean*200*stats::runif(1,0.9,1),sigma_mean*500*stats::runif(1,0.9,1))
    combos_n<-expand.grid(theta_base_seq,sigma_base_seq,alpha_base_seq)
    if (nrow(combos_n)>500) {
      combos_t<-sample_n(combos_n,500)
    }
    else {
      combos_t<-combos_n
    }
  }

  if ("BM" %in% cache$epc_params) {
    theta_mean=mean(cache$PCM.X)*stats::runif(1,0.9,1)
    theta_base_seq=base::seq(theta_mean,theta_mean*10,theta_mean)*stats::runif(1,0.9,1)
    sigma_mean=var(as.numeric(cache$PCM.X))/max(nodeHeights(cache$PCMtree))
    sigma_base_seq=c(sigma_mean*stats::runif(1,0.9,1),sigma_mean*5*stats::runif(1,0.9,1),sigma_mean*10*stats::runif(1,0.9,1),sigma_mean*50*stats::runif(1,0.9,1),sigma_mean*100*stats::runif(1,0.9,1),sigma_mean*200*stats::runif(1,0.9,1),sigma_mean*500*stats::runif(1,0.9,1))
    combos_n<-expand.grid(theta_base_seq,sigma_base_seq)
    if (nrow(combos_n)>500) {
      combos_t<-sample_n(combos_n,500)
    }
    else {
      combos_t<-combos_n
    }
  }

  if (all(c("theta","alpha") %in%  cache$epc_params)) {
    if (cache$m_type=="linear") {
      b0b1a_unique<-b0b1a_gen(cache)
      b0b1t_unique<-b0b1t_gen(cache)

      sigma_mean=var(as.numeric(cache$PCM.X))/max(nodeHeights(cache$PCMtree))
      sigma_base_seq=c(sigma_mean*stats::runif(1,0.9,1),sigma_mean*5*stats::runif(1,0.9,1),sigma_mean*10*stats::runif(1,0.9,1),sigma_mean*50*stats::runif(1,0.9,1),sigma_mean*100*stats::runif(1,0.9,1),sigma_mean*200*stats::runif(1,0.9,1),sigma_mean*500*stats::runif(1,0.9,1))

      combos_n_pre<-expand.grid(b0b1t_unique[,1],sigma_base_seq,b0b1a_unique[,1])

      colnames(combos_n_pre)<-c("b0t","sigma","b0a")

      combos_n_pre_add <- dplyr::left_join(combos_n_pre, b0b1a_unique, by = "b0a")
      combos_n_pre_add <- dplyr::left_join(combos_n_pre_add , b0b1t_unique, by = "b0t")
      combos_n<-as.data.frame(combos_n_pre_add)

      combos_n<-combos_n[c(1,5,2,3,4)] #reorder

      rownames(combos_n) = base::seq(length=nrow(combos_n))
      if (nrow(combos_n)>500) {
        combos_t<-sample_n(combos_n,500)
      }
      else {
        combos_t<-combos_n
      }
    }

    if (cache$m_type=="step") {
      thresh_t_base_seq=base::seq(min(cache$environment),max(cache$environment),(max(cache$environment)-min(cache$environment))/10)
      step_t_base_seq=stats::runif(10,min(cache$PCM.X),max(cache$PCM.X))
      step2_t_base_seq=stats::runif(10,min(cache$PCM.X),max(cache$PCM.X))
      present_tmp<-data.frame(time=tip_root_length[tip_root_length > max(tip_root_length) - 0.0000000001])
      tip_data<-data.frame(tip_data=t(cache$PCM.X))
      present_data<-base::merge(present_tmp, tip_data, by = 'row.names', all =FALSE)
      sigma_mean=var(as.numeric(present_data$tip_data))/max(nodeHeights(cache$PCMtree))
      sigma_base_seq=c(sigma_mean,sigma_mean*5,sigma_mean*10,sigma_mean*50,sigma_mean*100,sigma_mean*200)
      thresh_a_base_seq=base::seq(min(cache$environment),max(cache$environment),(max(cache$environment)-min(cache$environment))/10)
      step_a_base_seq=base::seq(0,0.5,0.05)
      step2_a_base_seq=base::seq(0,0.5,0.05)
      combos_n<-expand.grid(thresh_t_base_seq,step_t_base_seq,step2_t_base_seq,sigma_base_seq,thresh_a_base_seq,step_a_base_seq,step2_a_base_seq)
      if (nrow(combos_n)>500) {
        combos_t<-sample_n(combos_n,500)
      }
      else {
        combos_t<-combos_n
      }
    }
  }

  if (!"theta" %in% cache$epc_params && "alpha" %in% cache$epc_params) {
    if (cache$m_type=="linear") {

      b0b1a_unique<-b0b1a_gen(cache)

      theta_mean=mean(cache$PCM.X)*stats::runif(1,0.7,1)
      theta_base_seq=base::seq(theta_mean/2,theta_mean*2,theta_mean)*stats::runif(1,0.7,1)
      sigma_mean=var(as.numeric(cache$PCM.X))/max(nodeHeights(cache$PCMtree))
      sigma_base_seq=c(sigma_mean,sigma_mean*5*stats::runif(1,0.9,1),sigma_mean*10*stats::runif(1,0.9,1),sigma_mean*50*stats::runif(1,0.9,1),sigma_mean*100*stats::runif(1,0.9,1),sigma_mean*200*stats::runif(1,0.9,1),sigma_mean*500*stats::runif(1,0.9,1))

      combos_n_pre<-expand.grid(theta_base_seq,sigma_base_seq,b0b1a_unique[,1])

      colnames(combos_n_pre)<-c("theta","sigma","b0a")

      combos_n_pre_add <- dplyr::left_join(combos_n_pre, b0b1a_unique, by = "b0a")

      combos_n <- combos_n_pre_add %>%
        group_by(b0a) %>%
        sample_n(1) %>%
        ungroup()

      combos_n<-as.data.frame(combos_n)
      rownames(combos_n) = base::seq(length=nrow(combos_n))
      if (nrow(combos_n)>500) {
        combos_t<-sample_n(combos_n,500)
      }
      else {
        combos_t<-combos_n
      }
    }

    if (cache$m_type=="step") {
      theta_mean=mean(cache$PCM.X)
      theta_base_seq=base::seq(theta_mean/2,theta_mean*2,theta_mean)
      sigma_mean=var(as.numeric(cache$PCM.X))
      sigma_base_seq=stats::runif(10,sigma_mean/max(nodeHeights(cache$PCMtree)),sigma_mean*5)
      thresh_a_base_seq=base::seq(min(cache$environment),max(cache$environment),(max(cache$environment)-min(cache$environment))/10)
      step_a_base_seq=base::seq(0,0.5,0.05)
      step2_a_base_seq=base::seq(0,0.5,0.05)
      combos_n<-expand.grid(theta_base_seq,sigma_base_seq,thresh_a_base_seq,step_a_base_seq,step2_a_base_seq)
      if (nrow(combos_n)>500) {
        combos_t<-sample_n(combos_n,500)
      }
      else {
        combos_t<-combos_n
      }    }
  }

  if ("theta" %in% cache$epc_params && !"alpha" %in% cache$epc_params) {
    if (cache$m_type=="linear") {


      b0b1t_unique<-b0b1t_gen(cache)
      sigma_mean=var(as.numeric(cache$PCM.X))/max(nodeHeights(cache$PCMtree))
      sigma_base_seq=c(sigma_mean,sigma_mean*5*stats::runif(1,0.9,1),sigma_mean*10*stats::runif(1,0.9,1),sigma_mean*50*stats::runif(1,0.9,1),sigma_mean*100*stats::runif(1,0.9,1),sigma_mean*200*stats::runif(1,0.9,1))
      alpha_base_seq<-base::seq(0,max_alpha,max_alpha/10)[1:10]

      combos_n_pre<-expand.grid(b0b1t_unique[,1],sigma_base_seq,alpha_base_seq)
      colnames(combos_n_pre)<-c("b0t","sigma","alpha")

      combos_n_pre_add <- dplyr::left_join(combos_n_pre, b0b1t_unique, by = "b0t")
      combos_n<-combos_n_pre_add[c(1,4,2,3)] #reorder

      if (nrow(combos_n)>500) {
        combos_t<-sample_n(combos_n,500)
      }
      else {
        combos_t<-combos_n
      }
    }
    if (cache$m_type=="step") {
      thresh_t_base_seq=base::seq(min(cache$environment),max(cache$environment),(max(cache$environment)-min(cache$environment))/10)
      step_t_base_seq=stats::runif(10,min(cache$PCM.X),max(cache$PCM.X))
      step2_t_base_seq=stats::runif(10,min(cache$PCM.X),max(cache$PCM.X))
      sigma_mean=var(as.numeric(cache$PCM.X))
      sigma_base_seq=stats::runif(10,sigma_mean/max(nodeHeights(cache$PCMtree)),sigma_mean*5)
      alpha_base_seq<-stats::runif(5,0,0.5)
      combos_n<-expand.grid(thresh_t_base_seq,step_t_base_seq,step2_t_base_seq,sigma_base_seq,alpha_base_seq)
      if (nrow(combos_n)>500) {
        combos_t<-sample_n(combos_n,500)
      }
      else {
        combos_t<-combos_n
      }    }
  }

  suppressWarnings({
    lik_values<-list()
    for (i in 1:nrow(combos_t)) {
      skip_to_next <- FALSE
      tryCatch(lik_values[i]<-epcTools::base.lik(cache,as.numeric(combos_t[i,])), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { lik_values[i] <-99999999999 }
    }
    combos_t$lik<-lik_values
  })

  #print(start_params)
  combos_t<-as.data.frame(lapply(combos_t, unlist))
  combos_re<-unique(expand.grid(head(combos_t[order(combos_t$lik),],n=3L))[,1:ncol(combos_t)-1])

  suppressWarnings({
    lik_values<-list()
    for (i in 1:nrow(combos_re)) {
      lik_values[i]<-epcTools::base.lik(cache,as.numeric(combos_re[i,]))
    }
    combos_re$lik<-lik_values
  })
  start_params<-combos_re[which.min(combos_re$lik),]
  return(start_params[-length(start_params)])
}

