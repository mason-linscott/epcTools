#' EPC environment-trait relationship visualization function
#'
#' Visualizes epc.max.lik outputs and reconstructed relationships such as the stationary variance through time or shifts in optimal values
#'
#' @param TS A time series
#'
#' @return A phylogeny constructed from the time series.
#'
#' @export
ts.tree.maker<-function(TS) {

  TS_increasing<-sort(TS[,1],decreasing=TRUE)

  node_number<-length(TS_increasing)*2-2

  edge_lengths<-data.frame( matrix(ncol = 1, nrow =2*(length(TS_increasing)-1)))
  colnames(edge_lengths) <- "edge_lengths"
  j=1
  for(i in 1:(length(TS_increasing)-1)) {
    edge_lengths[j,1]<-TS_increasing[i]-TS_increasing[i+1]
    j=j+1
    edge_lengths[j,1]<-0.000001
    j=j+1
  }

  j=1
  edges<-matrix(ncol=2,nrow=length(TS_increasing)*2-2)
  for (i in seq(1,nrow(edges),by=2)) {
    edges[i,1]<-j+length(TS_increasing)
    j=j+1
    edges[i,2]<-j+length(TS_increasing)
  }

  j=1
  for (i in seq(2,nrow(edges),by=2)) {
    edges[i,1]<-j+length(TS_increasing)
    edges[i,2]<-j
    j=j+1
  }

  edges[nrow(edges)-1,2]<-length(TS_increasing)

  tip_labels<-paste0("t", 1:length(TS_increasing))
  N_node<-length(TS_increasing)-1
  TS_tree<-list(edge.length=as.numeric(edge_lengths$edge_lengths),Nnode=N_node,edge=edges,tip.labels=tip_labels)
  class(TS_tree)<-"phylo"
  TS_tree<-reorder(TS_tree,"cladewise")
  return(TS_tree)
}
