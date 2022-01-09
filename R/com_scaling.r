#' Performs scaling of functional or phylogenetic structure of simulated community.
#'
#' @param com A list, returned by evolve_com.
#' @param grain_size A number, size of community data matrix.
#' @param use_tr If True, uses trait distances, else uses phylogenetic distances.
#' @param abundance_weighted Boolean, indicates, if abundances should be weighted.
#' @param scaling_iter A number of iterations.
#' @param be_verbose Boolean, returns messages during the execution.
#'
#' @return A list, result of scaling.
#' @export
#'
#' @examples
#' com_scaling(com, grain_size = 8, use_tr = T, abundance_weighted = T, scaling_iter = 100, be_verbose = T)


com_scaling <- function(com, grain_size = 8, use_tr = T, abundance_weighted = T, scaling_iter = 100, be_verbose = T){

  phylo <- com$phylo
  com <- com[-16]

  cda <- make_cda(com$id, grain_size)
  if (use_tr){
    dist <- as.matrix(dist(com$niche))
  } else {
    dist <- cophenetic(phylo)
  }


  if (be_verbose){
    print('Perfoming scaling.')
    print(system.time(resG <- scaling_grain(com$id, dist, abundance_weighted, scaling_iter)))
    print(system.time(resP <- scaling_pool(cda, dist, abundance_weighted, scaling_iter)))
    print(system.time(resC <- scaling_coherent(cda, dist, abundance_weighted, scaling_iter)))
  } else {
    resG <- scaling_grain(com$id, dist, abundance_weighted, scaling_iter)
    resP <- scaling_pool(cda, dist, abundance_weighted, scaling_iter)
    resC <- scaling_coherent(cda, dist, abundance_weighted, scaling_iter)
  }

  res <- list(resG = numeric(), resP = numeric(), resC = numeric())

  res$resG <- data.frame(r = resG[, 1], ses = resG[, 2])
  res$resP <- data.frame(r = 1:max(resP[, 1]), ses = tapply(resP[, 5], resP[, 1], mean, na.rm = T))
  res$resC <- data.frame(r = 0:max(resC[, 1]), ses = tapply(resC[, 5], resC[, 1], mean, na.rm = T))




  return(res)

}
