#' Creates square landscape on which community will be simulated.
#'
#' @param size Size of square landscape.
#' @param roughness Level of spatial autocorrelation.
#'
#' @return Landscape
#' @export
#'
#' @examples
#' sim_arena(size = 256, roughness = 0.75)

sim_arena <- function(size = 256, roughness = 0.75){
  size <- size + 2
  md <- NLMR::nlm_mpd(ncol = size, nrow = size, roughness = roughness, torus = T)
  val <- raster::values(md)
  mat <- matrix(val, nrow = sqrt(length(val)))
  if ((size %% 2) == 0){
    mat <- mat[-nrow(mat), -ncol(mat)] * 100
  } else {
    mat <- mat * 100
  }

  return(mat)
}
