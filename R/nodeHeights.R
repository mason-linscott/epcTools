#' Support function to get node heights
#'
#' @param x tree
#'
#' @export

.nodeHeights <- function(tree, scale="YBP"){
  branchTimes <- signif(BAMMtools:::NU.branching.times(tree),5)
  TH <- max(branchTimes)
  n.tip <- length(tree$tip.label)
  branchTimes <- c(setNames(rep(0, n.tip), 1:n.tip), branchTimes)
  nH <- cbind(branchTimes[tree$edge[,1]], branchTimes[tree$edge[,2]])
  if(scale=="YBP"){
    return(nH)
  }
  if(scale=="ABS"){
    return(TH-nH)
  }
}
