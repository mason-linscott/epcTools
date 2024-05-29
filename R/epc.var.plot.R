#' EPC environment-trait relationship visualization function
#'
#' Visualizes epc.max.lik outputs and reconstructed relationships such as the stationary variance through time or shifts in optimal values
#'
#' @param max_output Output from epc.max.lik
#'
#' @param cache An EPC cache object.
#'
#' @return An epc.var.plot that depicts the optima and stationary variance shifting through time under the maximum likelihood estimate
#'
#' @export



epc.var.plot<-function(max_output,cache){

  if (!is.null(max_output$dentist_output)) {

    sig2_core<-max_output$dentist_output$all_ranges$sig2[1]

    if (all(c("theta","alpha") %in%  cache$epc_params)) {

      if (cache$m_type=="linear") {
        theta_core<-max_output$dentist_output$all_ranges$b0_t[1]+max_output$dentist_output$all_ranges$b1_t[1]*cache$environment
        alpha_core<-max_output$dentist_output$all_ranges$b0_a[1]+max_output$dentist_output$all_ranges$b1_a[1]*cache$environment

      }

      if (cache$m_type=="step") {
        theta_core <- ifelse(cache$environment>max_output$dentist_output$all_ranges$step_thresh_t[1],max_output$dentist_output$all_ranges$step_high_t[1],max_output$dentist_output$all_ranges$step_low_t[1])
        alpha_core <- ifelse(cache$environment>max_output$dentist_output$all_ranges$step_thresh_a[1],max_output$dentist_output$all_ranges$step_high_a[1],max_output$dentist_output$all_ranges$step_low_a[1])
      }
    }

    if (!"theta" %in% cache$epc_params && "alpha" %in% cache$epc_params) {

      if (cache$m_type=="linear") {
        alpha_core<-max_output$dentist_output$all_ranges$b0_a[1]+max_output$dentist_output$all_ranges$b1_a[1]*cache$environment
        theta_core<-rep(max_output$dentist_output$all_ranges$theta[1],length(alpha_core))
      }

      if (cache$m_type=="step") {
        alpha_core <- ifelse(cache$environment>max_output$dentist_output$all_ranges$step_thresh_a[1],max_output$dentist_output$all_ranges$step_high_a[1],max_output$dentist_output$all_ranges$step_low_a[1])
        theta_core<-rep(max_output$dentist_output$all_ranges$theta[1],length(alpha_core))
      }
    }


    if ("theta" %in% cache$epc_params && !"alpha" %in% cache$epc_params) {

      if (cache$m_type=="linear") {
        theta_core<-max_output$dentist_output$all_ranges$b0_t[1]+max_output$dentist_output$all_ranges$b1_t[1]*cache$environment
        alpha_core<-rep(max_output$dentist_output$all_ranges$alpha[1],length(theta_core))
      }

      if (cache$m_type=="step") {
        theta_core <- ifelse(cache$environment>max_output$dentist_output$all_ranges$step_thresh_t[1],max_output$dentist_output$all_ranges$step_high_t[1],max_output$dentist_output$all_ranges$step_low_t[1])
        alpha_core<-rep(max_output$dentist_output$all_ranges$alpha[1],length(theta_core))
      }

    }
  }

  else
  {
    sig2_core<-max_output$sig2[1]

    if (all(c("theta","alpha") %in%  cache$epc_params)) {

      if (cache$m_type=="linear") {
        theta_core<-max_output$b0_t[1]+max_output$b1_t[1]*cache$environment
        alpha_core<-max_output$b0_a[1]+max_output$b1_a[1]*cache$environment

      }

      if (cache$m_type=="step") {
        theta_core <- ifelse(cache$environment>max_output$step_thresh_t[1],max_output$step_high_t[1],max_output$step_low_t[1])
        alpha_core <- ifelse(cache$environment>max_output$step_thresh_a[1],max_output$step_high_a[1],max_output$step_low_a[1])
      }
    }

    if (!"theta" %in% cache$epc_params && "alpha" %in% cache$epc_params) {

      if (cache$m_type=="linear") {
        alpha_core<-max_output$b0_a[1]+max_output$b1_a[1]*cache$environment
        theta_core<-rep(max_output$theta[1],length(alpha_core))
      }

      if (cache$m_type=="step") {
        alpha_core <- ifelse(cache$environment>max_output$step_thresh_a[1],max_output$step_high_a[1],max_output$step_low_a[1])
        theta_core<-rep(max_output$theta[1],length(alpha_core))
      }
    }


    if ("theta" %in% cache$epc_params && !"alpha" %in% cache$epc_params) {

      if (cache$m_type=="linear") {
        theta_core<-max_output$b0_t[1]+max_output$b1_t[1]*cache$environment
        alpha_core<-rep(max_output$alpha[1],length(theta_core))
      }

      if (cache$m_type=="step") {
        theta_core <- ifelse(cache$environment>max_output$step_thresh_t[1],max_output$step_high_t[1],max_output$step_low_t[1])
        alpha_core<-rep(max_output$alpha[1],length(theta_core))
      }

    }
  }

  stat_var<-sig2_core/(2*alpha_core)

  mycol <- rgb(0, 0, 255, max = 255, alpha = 20, names = "blue50")
  tree_max<-max(nodeHeights(cache$phy)[,2])
  time_periods<-seq(0,tree_max,tree_max/cache$n_slice)
  xlim_plot<-c(0,max(nodeHeights(cache$phy)))

  ##BOX
  tops<-theta_core+(stat_var/2)
  bots<-theta_core-(stat_var/2)
  rights<-time_periods[-(length(time_periods))]
  lefts<-time_periods[-1]
  phenogram(cache$PCMtree, t(cache$PCM.X)[,1], ftype="off",ylim=range(tops,bots,cache$PCM.X),color=mycol)
  rect(lefts,bots,rights,tops)
}
