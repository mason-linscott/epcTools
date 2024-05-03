#' Basic likelihood function
#'
#' Returns a likelihood value given a cache object and vector of parameters
#'
#' @param slice A vector of slice times for the tree.
#'
#' @param trees A phylogenetic tree.
#'
#' @return A sliced tree for loading into epoch simmap
#'
#' @export
getBranchesSlice <- function(slice, tree, nH1,nH2){
  present <- unname(which(apply(nH1, 1, function(x) slice < x[2] & slice > x[1])))
  locs <- (slice) - nH1[present,1]
  return(list(branches=present, nodes=tree$edge[present,2], locl=unname(locs)))
}

