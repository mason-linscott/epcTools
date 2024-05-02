#' Support function to load pars for epoch time slicing
#'
#' @param pars Numeric vector of parameters to be used to determine simmap positions
#'
#' @param tree A phylogenetic tree.
#'
#' @return Vector of values to be passed for simmap creation
#'
#' @export


pars3simmap<-function (pars, tree)
{
  tree <- reorder(tree, "postorder")
  sb <- pars$sb
  loc <- pars$loc
  t2 <- pars$t2
  if (!all(pars$sb %in% 1:nrow(tree$edge)))
    stop("Invalid parameter list. Specified branches not found in the tree")
  if (!all(pars$loc < tree$edge.length[pars$sb]))
  {
    loc<-pars$loc[which(pars$loc < tree$edge.length[pars$sb])]
    sb<-pars$sb[which(pars$loc < tree$edge.length[pars$sb])]
    t2 <- pars$t2[which(pars$loc < tree$edge.length[pars$sb])]
  }
  Th <- NULL
  nbranch <- length(tree$edge.length)
  maps <- lapply(tree$edge.length, function(x) {
    y <- x
    names(y) <- 1
    y
  })
  dup <- which(duplicated(sb))
  if (pars$k > 0) {
    if (length(dup) > 0) {
      maps[sb[-dup]] <- lapply(1:length(sb[-dup]), .addshift2map,
                               maps = maps, sb = sb[-dup], loc = loc[-dup],
                               t2 = t2[-dup])
    }
    else {
      maps[sb] <- lapply(1:length(sb), .addshift2map,
                         maps = maps, sb = sb, loc = loc, t2 = t2)
    }
    for (i in dup) {
      maps[[sb[i]]] <- .addshift2map(i, maps = maps, sb = sb,
                                     loc = loc, t2 = t2)
    }
    nopt <- rep(1, nbranch)
    for (i in nbranch:1) {
      if (i %in% sb) {
        opt <- as.integer(names(maps[[i]])[length(maps[[i]])])
        nopt[tree$edge[i, 2]] <- opt
        names(maps[[i]])[1] <- nopt[tree$edge[i, 1]]
      }
      else {
        names(maps[[i]])[1] <- nopt[tree$edge[i, 1]]
        nopt[tree$edge[i, 2]] <- nopt[tree$edge[i, 1]]
      }
    }
    shiftdown <- nopt[tree$edge[, 1]]
    new.maps <- lapply(1:nbranch, function(x) {
      names(maps[[x]])[1] <- shiftdown[x]
      maps[[x]]
    })
    new.maps <- maps
    for (j in 1:nbranch) {
      names(new.maps[[j]])[1] <- shiftdown[j]
    }
    anc.theta <- unlist(lapply(new.maps[sb], function(x) as.integer(names(x)[length(x) -
                                                                               1])), F, F)
    o <- rev(order(sb, loc * -1))
    shifted.maps <- new.maps[sb[o]]
    t1 <- rep(NA, length(t2))
    for (i in 1:length(t2)) {
      nm <- as.integer(names(maps[[sb[o][i]]]))
      t1[nm[2:length(nm)] - 1] <- nm[1:(length(nm) - 1)]
      Th[t2[o[i]]] <- Th[t1[o[i]]]
    }
  }
  else {
    new.maps <- maps
  }
  new.tree <- tree
  new.tree$maps <- new.maps
  new.pars <- pars
  col <- c(1, rainbow(pars$k))
  names(col) <- 1:(pars$k + 1)
  return(list(tree = new.tree, pars = new.pars, col = col))
}

.addshift2map <- function(x,maps=maps,sb=sb,loc=loc,t2=t2){
  m <- maps[[sb[x]]]
  cs.m <- cumsum(m)
  o <- min(which(cs.m>loc[x]))
  if(o==1){
    m[o] <- loc[x]
  } else {
    m[o] <- loc[x]-cs.m[o-1]
  }
  new.m <- cs.m[o]-loc[x]
  names(new.m) <- t2[x]
  M <- c(m[1:o],new.m,m[1:length(m) > o])
  return(M)
}

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

