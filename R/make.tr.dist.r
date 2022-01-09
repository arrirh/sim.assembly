#' Makes matrix of functional distances.
#'
#' @param com A list, returned by evolve_com.
#'
#' @return Matrix of functional distances.
#' @export
#'
#' @examples
#' make.tr.dist(com)

make.tr.dist <- function(com)
{
  tr <- data.frame(com$niche)
  rownames(tr) <- as.character(1:(length(com$n)))
  dist.tr <- as.matrix(daisy(tr, metric = "euclidean"))

  return(dist.tr)
}
