#' Creates phylogenetic tree and generates a single trait for each species.
#'
#' @param b Birth rate.
#' @param d Death rate.
#' @param n_taxa A number of species in the tree.
#' @param model A model, under which trait should be generated ("BM" or "OU").
#' @param sigma A standard-deviation of the random component for each branch.
#' @param alpha A strength of the selective constraint for each branch in OU model.
#'
#' @return List consisting of phylogenetic tree, matrix of phylogenetic distances, trait and matrix of trait distances.
#' @export
#'
#' @examples
#' sim_tree_n_traits(b = 0.1, d = 0, n_taxa = 500, model = "BM", sigma = 0.1, alpha = 0.1)

sim_tree_n_traits <- function(b = 0.1, d = 0, n_taxa = 500, model = "BM", sigma = 0.1, alpha = 0.1){
  tree <- geiger::sim.bdtree(b = b, d = d, stop = "taxa", n = n_taxa)
  dist.phy <- cophenetic(tree)
  colnames(dist.phy) <- rownames(dist.phy) <- 1:n_taxa

  tr <- ape::rTraitCont(tree, model = model, sigma = sigma, alpha = alpha)
  trait <- (tr - min(tr)) / (max(tr) - min(tr)) * 100
  names(trait) <- 1:length(trait)
  tr.dist <- as.matrix(dist(trait))

  res <- list(phylo = tree, phylo_dist = dist.phy, trait = trait, trait_dist = tr.dist)

  return(res)
}
