#' Drop a set of random fossil tips
#'
#' @param phy A phylogenetic tree.
#'
#' @param m Number of tips to drop
#'
#' @return  A tree with dropped extant tips
#'
#' @export

drop.random.fossil <- function(phy, m, tol = 1e-8)
{
  n <- Ntip(phy)
  x<-castor::get_all_distances_to_root(phy)[1:n]
  drop.tip(phy, phy$tip.label[sample(which(x < max(x) - tol),m)[1:m]])
}


