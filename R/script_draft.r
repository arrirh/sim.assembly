#' Supporting function for scaling_pipeline. Creates a Dataframe containing all parameters in model.
#'
#' @param n A number of rows.
#'
#' @return Dataframe.
#' @export
#'
#' @examples
#' script_draft(8)


script_draft <- function(n){
  script <- data.frame(arena_size = numeric(n), arena_roughness = numeric(n),
                       phylo_b = numeric(n), phylo_d = numeric(n), phylo_n_taxa = numeric(n),
                       trait_model = character(n), trait_sigma = numeric(n), trait_alpha = numeric(n),
                       com_range = numeric(n), betaEnv = numeric(n), betaComp = numeric(n),
                       com_nu = numeric(n), com_sigma = numeric(n), sim_n_iter = numeric(n),
                       grain_size = numeric(n), abundance_weighted = logical(n),
                       scaling_iter = numeric(n))
  return(script)
}
