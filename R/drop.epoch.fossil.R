#' Drop all tips from a time epoch
#'
#' @param phy A phylogenetic tree.
#'
#' @param m Lower bound of time to drop tips from
#'
#' @param M Upper bound of time to drop tips from
#'
#' @return  A tree with dropped extant tips
#'
#' @export

drop.epoch.fossil <- function(phy, m, M, tol = 1e-8)
{
  n <- Ntip(phy)
  x<-castor::get_all_distances_to_root(phy)[1:n]
  t<- max(x) -x
  drop.tip(phy, phy$tip.label[which(t > m + tol & t < M)])
}



