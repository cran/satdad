% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/satdad_Rfunctions.R
\name{ellEmp}
\alias{ellEmp}
\title{Empirical stable tail dependence function.}
\usage{
ellEmp(sample, x, k)
}
\arguments{
\item{sample}{A \code{(n times d)} matrix.}

\item{x}{A \code{(N.x times d)} matrix.}

\item{k}{A vector of \code{N.k} integers smaller or equal to \code{n}.}
}
\value{
A \code{(N.k times N.x)} matrix is returned.
}
\description{
The stable tail dependence function of \code{sample} is estimated at each row of \code{x} and for all values of the threshold parameter \code{k}.
}
\examples{

## Fix a 5-dimensional asymmetric tail dependence structure
ds5 <- gen.ds(d = 5)

## Construct a 1000-sample of Mevlog random vector associated with ds5
sample5 <- rMevlog(n = 1000, ds = ds5)

## Select 3 vectors in R^5
x5 <- matrix(runif(5*3), ncol = 5)

## Select 4 values for the threshold parameter
k5 <- (2:5)*10

## Estimation of the stable tail dependence function
# We thus get a 4 x 3 matrix
ellEmp(sample = sample5, x = x5, k = k5)

## Theoretical values of the stable tail dependence function inherited from ds5
ellMevlog(x = x5, ds = ds5)

}
\references{
Huang, X. (1992).
Statistics of bivariate extremes.
PhD Thesis, Erasmus University Rotterdam, Tinbergen Institute Research series No. 22.

de Haan, L. and Resnick, S. I. (1993).
Estimating the limit distribution of multivariate extremes. Communications in Statistics.
Stochastic Models 9, 275--309.

Fougeres, A.-L., de Haan, L. and  Mercadier, C.  (2015).
Bias correction in multivariate extremes.
Annals of Statistics 43 (2), 903--934.
}
\seealso{
\code{\link[satdad]{ellMevlog}}, \code{\link[satdad]{gen.ds}}
}
\author{
Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
}
