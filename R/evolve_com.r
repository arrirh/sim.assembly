#' Simulated previously generated community.
#'
#' @param com A list, returned by create_com.
#' @param n_iter A number of iterations of simulation.
#' @param be_verbose Boolean, returns messages during the execution.
#' @param verb_iter A number of iterations in community simulation after which the message is returned.
#'
#' @return A list, community data.
#' @export
#'
#' @examples
#' evolve_com(com, n_iter = 2000, be_verbose = T, verb_iter = 500)

evolve_com <- function(com, n_iter = 2000, be_verbose = T, verb_iter = 500){
  if (be_verbose){
    print(system.time(coms <- cc_simCom(com, n_iter, verb_iter)))
  } else {
    coms <- cc_simCom(com, n_iter, n_iter)
  }

  return(coms)
}
