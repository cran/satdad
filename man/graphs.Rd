% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/satdad_Rfunctions.R
\name{graphs}
\alias{graphs}
\title{Graphs of the tail dependence structure for Mevlog models.}
\usage{
graphs(
  ds,
  names = NULL,
  n.MC = 1000,
  which = "taildependograph",
  random = FALSE,
  thick.td = 5,
  thick.ec = 5
)
}
\arguments{
\item{ds}{An object of class \code{ds}.}

\item{names}{A character vector of length \code{d} which replaces \code{as.character(1:d)} (the default ones).}

\item{n.MC}{Monte Carlo sample size. Default value is 1000. See details in  \code{\link[satdad]{tsic}}.}

\item{which}{A character string:  \code{taildependograph} (the default), \code{iecgraph}, or \code{both},}

\item{random}{A boolean. `FALSE` (the default): the vertex positions are fixed along a circle.  `TRUE`: some randomness is applied for positioning the vertices.}

\item{thick.td}{A numeric value for the maximal thickness of edges in \code{taildependograph}. Default value is 5.}

\item{thick.ec}{A numeric value for the maximal thickness of edges in \code{iecgraph}. Default value is 5.}
}
\value{
The function returns either the tail dependograph or the inverse extremal coefficients graph, or both, for the tail dependence structure `ds`.

The tail dependograph displays pairwise tail superset importance coefficients, which measure the extent to which pairs of components (and their supersets) contribute to the overall variance of the stable tail dependence function.
We refer to  Mercadier, C. and Roustant, O. (2019) for more details. These coefficients are computed using the `tsic` function with the `"ind = 2"` option.

The inverse extremal coefficients graph shows the inverse renormalized pairwise coefficients computed as \eqn{\theta_{ij}=1-\ell(1_i,1_j,\bold{0})/2}.
}
\description{
Tail dependograph and Inverse extremal coefficients graph  for Mevlog models.   A \code{Mevlog} model is a multivariate extreme value (symmetric or asymmetric) logistic model.
}
\details{
The tail dependence structure is set by a \code{ds} object. It thus corresponds to the stable tail dependence function \eqn{\ell}. The way to deduce the stable tail dependence function \eqn{\ell} from \code{ds} is explained in the Details section of \code{\link[satdad]{gen.ds}}.
}
\examples{

## Fix a 8-dimensional asymmetric tail dependence structure
ds8 <- gen.ds(d = 8)

## Plot the graphs that illustrate  characteristics of the tail dependence structure
graphs(ds = ds8, which = "both")

}
\references{
Mercadier, C. and Roustant, O. (2019)
The tail dependograph.
\emph{Extremes}, \bold{22}, 343--372.

Tiago de Oliveira, J. (1962/63)
Structure theory of bivariate extremes, extensions.
\emph{Estudos de Matematica, Estatistica, e Economicos}, 7:165--195.

Smith, R. L. (1990)
Max-stable processes and spatial extremes.
\emph{Dept. of Math., Univ. of Surrey}, Guildford GU2 5XH, England.
}
\seealso{
\code{\link[satdad]{tsic}}, \code{\link[satdad]{ec}}, \code{\link[satdad]{ellMevlog}}
}
\author{
Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
}
