#' Plot dentist parameter outputs
#'
#' A lightly modified plot.dentist function for the epcTools framework
#'
#' @param epc_out An epc.max.lik output to be plotted
#'
#' @param local.only If only the 95% CI should be plotted.
#'
#' @return  A plot corresponding to plot.dentist with 95% confidence values for each parameter plotted.
#'
#' @export
plot.epc<-function(epc_out,local.only=FALSE, ...){
  nparams <- ncol(epc_out$dentist_output$results)-1
  nplots <- nparams + (nparams^2 - nparams)/2
  results <- epc_out$dentist_output$results
  threshold <- epc_out$dentist_output$best_neglnL + epc_out$dentist_output$delta
  results$color <- ifelse(results[,1]<=threshold, "black", "gray")
  results_outside <- subset(results, results$color=="gray")
  results_inside <- subset(results, results$color=="black")
  graphics::par(mfrow=c(ceiling(nplots/nparams), nparams))
  for (i in base::sequence(nparams)) {
    if(local.only){
      xlim=range(results_inside[,i+1])
      ylim=c(range(results_inside[,1]))
    }else{
      xlim=range(c(results_inside[,i+1]), results_outside[,i+1])
      ylim=c(range(epc_out$dentist_out$best_neglnL, epc_out$dentist_output$best_neglnL + (epc_out$dentist_output$delta*5)))
    }
    plot(results[,i+1], results[,1], pch=20, col=results$color, main=colnames(results)[i+1], xlab=colnames(results)[i+1], ylab="Negative Log Likelihood", bty="n", xlim=xlim, ylim=ylim, ...)
    graphics::abline(h=threshold, col="blue")
    graphics::points(results[which.min(results[,1]), i+1], results[which.min(results[,1]),1], pch=21, col="red")
  }
  for (i in base::sequence(nparams)) {
    for (j in base::sequence(nparams)) {
      if(j>i) {
        if(local.only){
          xlim=range(results_inside[,i+1])
          ylim=range(results_inside[,j+1])
        }else{
          xlim=range(c(results_inside[,i+1]), results_outside[,i+1])
          ylim=range(c(results_inside[,j+1]), results_outside[,j+1])
        }
        plot(results_outside[,i+1], results_outside[,j+1], pch=20, col=results_outside$color, xlab=colnames(results)[i+1], ylab=colnames(results)[j+1], bty="n", main=paste0(colnames(results)[j+1], " vs. ", colnames(results)[i+1]), xlim=xlim, ylim=ylim, ...)
        graphics::points(results_inside[,i+1], results_inside[,j+1], pch=20, col=results_inside$color)
        graphics::points(results[which.min(results[,1]), i+1], results[which.min(results[,1]),j+1], pch=21, col="red")
      }
    }
  }
  graphics::par(mfrow=c(1,1))
}

