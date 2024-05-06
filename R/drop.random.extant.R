#' Drop random extant tip
#'
#' @param phy A phylogenetic tree.
#'
#' @param m Number of tips to drop
#'
#' @return  A tree with dropped extant tips
#'
#' @export

drop.random.extant <- function(phy, m, tol = 1e-8)
{
  n <- Ntip(phy)
  x<-castor::get_all_distances_to_root(phy)[1:n]
  drop.tip(phy, phy$tip.label[sample(which(near(x,max(x))=="TRUE"),m)[1:m]])
}

