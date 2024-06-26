#' Summarize EPC run
#'
#' A summary function for the output of epc.max.lik
#'
#' @param epc_out An epc.max.lik output to be plotted
#'
#' @return  A summary of parameters and model performance for an epc.max.lik object
#'
#' @export
summary.epc<-function(epc_out){
  cat(paste0("An EPC - ", epc_out$fit$EPC_model, " model was ran where those parameters had an ",epc_out$fit$Type, " relationship with the environment", "\n\nMaximum likelihood estimate:\n", epc_out$fit$Lik, "\n\nNumber of Parameters:\n", epc_out$fit$Npar,"\n\nAIC:\n", epc_out$fit$AIC))
  cat("\n\nParameters Estimates: \n")
  print(epc_out$dentist_output$all_ranges)
}
