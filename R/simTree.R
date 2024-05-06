#' Simulate a phylogeny
#'
#' @param n_extant Number of extant tips.
#'
#' @param tree_size Size of tree to simulate
#'
#' @param max_time Max length of tree desired in time
#'
#' @param n_tree Number of trees to simulate
#'
#' @return A phylogeny
#'
#' @export
simTree<-function(n_extant,tree_size,max_time,n_tree)
{
  tree_list<-list() #store trees in list
  n_fossil=tree_size-n_extant
  b_modify<-0.1*(n_extant/tree_size)
  p_birth=0.2+b_modify
  p_death=0.2

  for(g in 1:n_tree) {

    if (n_extant>1)
    {
      pre_tree<-ape::rphylo(n=n_extant,birth=p_birth,death=p_death,fossils = TRUE) #simulate a tree with fossils
      while(Ntip(pre_tree)<tree_size) #make sure simtree is big enough
      {
        pre_tree<-ape::rphylo(n=n_extant,birth=p_birth,death=p_death,fossils = TRUE)
      }
      tree<-drop.random.fossil(pre_tree,(Ntip(pre_tree)-(tree_size)))
    }
    if (n_extant==tree_size) #drop fossils if simulating ultrametric trees
    {
      tree<-drop.random.fossil(pre_tree,(Ntip(pre_tree)-(tree_size)))
    }


    if (n_extant==0) #drop extant tips if simulating only fossils
    {
      n_extant=50
      pre_tree<-try(ape::rphylo(n=n_extant,birth=p_birth,death=p_death,fossils = TRUE),silent=TRUE)
      while((Ntip(pre_tree)-n_extant)<tree_size) #make sure simtree is big enough
      {
        pre_tree<-try(ape::rphylo(n=n_extant,birth=p_birth,death=p_death,fossils = TRUE),silent=TRUE)
      }
      near_tree<-drop.random.extant(pre_tree,n_extant)
      tree<-drop.random.fossil(near_tree, (Ntip(near_tree)-tree_size))

    }
    tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*max_time #scale to time

    tree_list[[g]]<-tree

  }
  if (length(tree_list)==1)
  {
    tree<-tree_list[[1]]
    return(tree)
  }
  else{
    return(tree_list)
  }
}

