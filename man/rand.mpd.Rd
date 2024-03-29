% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rand.mpd.r
\name{rand.mpd}
\alias{rand.mpd}
\title{Randomizes given community data and calculates MPD.}
\usage{
rand.mpd(cdm, dist, ab = F)
}
\arguments{
\item{cdm}{Community data matrix.}

\item{dist}{A matrix of phylogenetic or functional distances.}

\item{ab}{Boolean, indicates, if abundances should be weighted.}
}
\value{
Value of MPD.
}
\description{
Randomizes given community data and calculates MPD.
}
\examples{
rand.mpd(cdm, dist, ab = F)
}
