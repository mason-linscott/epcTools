#' Diatom timeseries data
#'
#' This is a default dataset included in the package.
#'
#' @docType data
#' @usage readRDS(diatom_data)
#' @format A list containing a timeseries, trait, environment, and a phylogenetic data from Yellowstone lake
#' \describe{
#'   \item{diatomTS}{A vector of times where diatom traits were measured}
#'   \item{diatom_traits}{A data frame of three traits that correspond to the timeseries measurements: spines, ribs, and diameter}
#'   \item{diatom_environment}{A dataframe where column 1 is time and column 2 is the environmental value (in this case biogenic silica flux)}
#'   \item{diatom_tree}{A phylogeny created using ts.tree.maker() on the timeseries dataset but included here to save a little time}

#' }
"diatom_data"
