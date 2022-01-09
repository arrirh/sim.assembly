#' Randomizes given community data and calculates MPD.
#'
#' @param cdm Community data matrix.
#' @param dist A matrix of phylogenetic or functional distances.
#' @param ab Boolean, indicates, if abundances should be weighted.
#'
#' @return Value of MPD.
#' @export
#'
#' @examples
#' rand.mpd(cdm, dist, ab = F)

rand.mpd <- function(cdm, dist, ab = F)
{
  perm <- sample(nrow(dist))
  colnames(cdm) <- colnames(cdm)[perm]
  mpd(cdm, dist, ab)
}
