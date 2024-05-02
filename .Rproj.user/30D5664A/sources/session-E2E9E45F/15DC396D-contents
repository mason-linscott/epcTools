#' Function for loading parameters into pcmbase format
#'
#' @param pars Numeric vector of parameters to be loaded.
#'
#' @param cache A cache object created using epcTools.

#' @return List of pars for loading into pcmbase.
#' @export

lik.pars2pcmbase <- function(pars,cache){
  if (all(c("theta","alpha") %in%  cache$epc_params)) {
    .p <- lapply(1:(cache$n_slice), function(x) c(pars$alpha[x], pars$theta[x], pars$sig2))
    p <- c(pars$theta[1], do.call(c, .p), 0)
    return(p)
  }

  if (!"theta" %in% cache$epc_params && "alpha" %in% cache$epc_params) {
    .p <- lapply(1:(cache$n_slice), function(x) c(pars$alpha[x], pars$theta, pars$sig2))
    p <- c(pars$theta[1], do.call(c, .p), 0)
    return(p)
  }

  if ("theta" %in% cache$epc_params && !"alpha" %in% cache$epc_params) {
    .p <- lapply(1:(cache$n_slice), function(x) c(pars$alpha, pars$theta[x], pars$sig2))
    p <- c(pars$theta[1], do.call(c, .p), 0)
    return(p)
  }
}

