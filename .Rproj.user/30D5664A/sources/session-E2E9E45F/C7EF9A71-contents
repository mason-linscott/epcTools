#' Make epoch simmaps
#'
#' Cuts a tree up into epochs with the same mapped simmap state across the tree
#'
#' From Josef Uyeda
#'
#' @param tree A phylogenetic tree.
#'
#' @param slices A vector where tree will be cut into slices
#'
#' @return An epoch simmap tree
#'
#' @export
makeEpochSimmap <- function(tree, slices){
  originalorder <- attributes(tree)$order
  potree <- reorder.phylo(tree, "postorder")
  nH1 <- .nodeHeights(potree, scale="ABS")
  nH2 <- .nodeHeights(potree, scale="YBP")
  sbs <- lapply(slices, function(x) getBranchesSlice(x, potree, nH1,nH2))
  .sb <- lapply(sbs, function(x) x$branches)
  N <- sapply(.sb, length)
  pars <- list()
  pars$sb <- do.call(c, .sb)
  pars$loc <- do.call(c, lapply(sbs, function(x) x$loc))
  pars$t2 <- do.call(c, lapply(1:length(slices), function(x) rep(x+1, N[x])))
  pars$k <- length(pars$sb)
  emap <- pars3simmap(pars, potree)$tree
  class(emap) <- c("simmap", "phylo")
  emap <- phytools::reorderSimmap(emap, originalorder)
  return(list(smap=emap, pars=pars))
}

