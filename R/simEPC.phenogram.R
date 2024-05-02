#' Plot simulated EPC phenogram
#'
#' A plot function that will plot a phenogram for simulated EPC cache objects
#'
#' @param caches A cache object created using epcTools.
#'
#' @param n Number of times to simulate and plot the data.
#'
#' @return  A phenogram plot of trait values through time.
#'
#' @export
simEPC.phenogram<-function(cache,n) {
  mycol <- rgb(0, 0, 255, max = 255, alpha = 6, names = "blue50")
  simdat <- PCMSim(cache$PCMtree, model=cache$pcmModel, X0=cache$pcmModel$X0[1])
  phenogram(cache$PCMtree, simdat[1,], ftype="off",color=mycol, xlab=("Time from root (My)"), ylab=("Trait value"))
  for (i in 1:n) {
    simdat <- PCMSim(cache$PCMtree, model=cache$pcmModel, X0=cache$pcmModel$X0[1])
    phenogram(cache$PCMtree, simdat[1,], ftype="off",color=mycol,add=TRUE)
    #lines(seq(0,100, length.out=11), log(2)/pars$alpha, lty=2, col="red")
  }
  #lines(seq(0,cache$max_time, length.out=cache$n_slice), log(2)/cache$pars$alpha, lty=2, col="red")
}


