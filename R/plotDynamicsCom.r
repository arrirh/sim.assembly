#' Visualizes generated community.
#'
#' @param com A list, returned by evolve_com.
#' @param com2 A list, returned by evolve_com.
#'
#' @export
#'
#' @examples
#' plotDynamics(com, com2 = NULL)

plotDynamics <- function(com, com2 = NULL)
{
    opar <- par(mfrow = c(2,2), col = "darkblue")

    plot(com$S, type = "l")
    if (!is.null(com2)) lines(com2$S, col = "tomato")
    plot(com$H, type = "l")
    if (!is.null(com2)) lines(com2$H, col = "tomato")
    plot(com$N, type = "l")
    if (!is.null(com2)) lines(com2$N, col = "tomato")
    plot(com$H/log(com$S), type = "l", ylab = "E")
    if (!is.null(com2)) lines(com2$H/log(com2$S), col = "tomato")
    par(opar)
}

#' Visualizes generated community.
#'
#' @param com A list, returned by evolve_com.
#'
#' @export
#'
#' @examples
#' plotCom(com)

plotCom <- function(com)
{
  nnn <- com$id
  nnn[] <- com$niche[nnn[]]

    opar <- par(mfrow = c(2,2), mar = c(4, 4, 0.25, 0.25))
    image(x = 1:com$L, y = 1:com$L, z = com$env, col  = topo.colors(10), asp = 1, axes = F, ann = F)
    image(x = 1:com$L, y = 1:com$L, z = nnn, col  = topo.colors(10), asp = 1, axes = F, ann = F)

    plot(sort(com$n[com$n>0], decr = T),type = "o", log = "y", pch = 21, bg = "skyblue", ylab = "abundance", xlab = "rank")
    hist(log(com$n[com$n > 0]), col = "wheat", xlab = "log(abundance)", main = "")
    par(opar)
}

