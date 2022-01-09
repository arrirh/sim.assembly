#' Makes community data matrices of given size from simulated community.
#'
#' @param com A list, returned by evolve_com.
#' @param size A number, size of each community data matrix.
#'
#' @return Community data matrices.
#' @export
#'
#' @examples
#' make.cdm(com, size = 8)

make.cdm <- function(com, size)
{
  n <- floor(com$L/size)
  cdm <- as.data.frame(matrix(0, ncol = length(com$niche), nrow = n^2))
  counter <- 1

  for (ii in 1:n) {
    for (jj in 1:n) {
      temp <- com$id[((ii-1)*size+1):(ii*size), ((jj-1)*size+1):(jj*size)]
      tab <- table(temp)
      id <- as.numeric(names(tab))
      cdm[counter, id] <- tab
      counter <- counter + 1
    }
  }

  return(cdm)
}

#' Makes community data array of given size from simulated community.
#'
#' @param com A list, returned by evolve_com.
#' @param size A number, size of each community data matrix.
#'
#' @return Community data array.
#' @export
#'
#' @examples
#' make.cda(com, size = 8)

make.cda <- function(com, size)
{
  n <- floor(com$L/size)
  cda <- array(data = 0, dim = c(n, n, length(com$niche)))

  for (ii in 1:n) {
    for (jj in 1:n) {
      temp <- com$id[((ii-1)*size+1):(ii*size), ((jj-1)*size+1):(jj*size)]
      tab <- table(temp)
      id <- as.numeric(names(tab))
      cda[ii, jj, id] <- tab
    }
  }

  return(cda)
}
