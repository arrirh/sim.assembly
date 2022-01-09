#' Generates and evolve community under given parameters.
#'
#' @param arena_size Size of square landscape.
#' @param arena_roughness Level of spatial autocorrelation of landscape.
#' @param phylo_b Birth rate of phylogenetic tree.
#' @param phylo_d Death rate of phylogenetic tree.
#' @param phylo_n_taxa A number of species in the tree.
#' @param trait_model A model, under which trait should be generated ("BM" or "OU").
#' @param trait_sigma A standard-deviation of the random component for each branch in trait simulation.
#' @param trait_alpha A strength of the selective constraint for each branch in OU model.
#' @param com_range A number, range in which modeling species interact.
#' @param betaEnv A number, coefficient for environmental filtering.
#' @param betaComp A number, coefficient for limiting similarity.
#' @param com_nu A number, rate of speciation.
#' @param com_sigma A number, niche width.
#' @param sim_n_iter A number of iterations of community simulation.
#' @param be_verbose Boolean, returns messages during the execution.
#' @param sim_verb_iter A number of iterations in community simulation after which the message is returned.
#'
#' @return A list, community data.
#' @export
#'
#' @examples
#' create_com(arena_size = 258, arena_roughness = 0.75,
#' phylo_b = 0.1, phylo_d = 0, phylo_n_taxa = 500,
#' trait_model = "BM", trait_sigma = 0.1, trait_alpha = 0.1,
#' com_range = 3, betaEnv = 0, betaComp = 0,
#' com_nu = 0.001, com_sigma = 10,
#' sim_n_iter = 2000, be_verbose = T, sim_verb_iter = 500)

create_com <- function(arena_size = 258, arena_roughness = 0.75,
                       phylo_b = 0.1, phylo_d = 0, phylo_n_taxa = 500,
                       trait_model = "BM", trait_sigma = 0.1, trait_alpha = 0.1,
                       com_range = 3, betaEnv = 0, betaComp = 0,
                       com_nu = 0.001, com_sigma = 10,
                       sim_n_iter = 2000, be_verbose = T, sim_verb_iter = 500){

  mat <- sim_arena(size = arena_size, roughness = arena_roughness)
  tree_n_traits <- sim_tree_n_traits(b = phylo_b, d = phylo_d, n_taxa = phylo_n_taxa,
                                     model = trait_model, sigma = trait_sigma, alpha = trait_alpha)

  com <- createCom(mat, tree_n_traits$trait, range = com_range,
                   betaEnv = betaEnv, betaComp = betaComp,
                   nu = com_nu, sigma = com_sigma)

  if (be_verbose){
    print('Community is generated. Initiating community simulation.')
  }
  com <- evolve_com(com, n_iter = sim_n_iter, be_verbose = be_verbose, verb_iter = sim_verb_iter)
  com <- append(com, list(tree_n_traits$phylo))
  names(com)[16] <- 'phylo'

  if (be_verbose){
    print('Community is simulated.')
  }

  return(com)

}
