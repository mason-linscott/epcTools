#' Support function to load pars for epoch time slicing taken from bayou
#'
#' @param x tree
#'
#' @export

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
