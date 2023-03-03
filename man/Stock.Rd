\name{Stock}
\alias{Stock}
\docType{data}
\title{Dataset. Yearly maxima of Log Returns of ten stock indices 1990-2015.}
\description{
This dataset consists of a matrix with years as rows and columns as stock indices. They appear in the following order: "SP500", "DJ", "NASD", "SMI", "EURS", "CAC", "DAX",  "HSI",  "SSEC",  "NIKK". A cell gives the yearly maximum of Log Returns of the associated stock indices. The latter values have been extracted from the R package \code{qrmdata} of Hofert, M., Hornik, K. and McNeil, A.J. (2019).
}
\references{

Hofert, M., Hornik, K. and McNeil, A.J. (2019).
qrmdata: Data Sets for Quantitative Risk Management Practice.
R package version 2019-12-03-1
URL https://CRAN.R-project.org/package=qrmdata.


Mercadier, C. and Roustant, O. (2019)
The tail dependograph.
\emph{Extremes}, \bold{22}, 343--372.


}
\seealso{
\code{\link[satdad]{graphsEmp}}
}
\examples{
data(Stock)

## We reproduce below Figure 7(a) of Mercadier and Roustant (2019).

graphsEmp(Stock, k = 26, which = "taildependograph", names = colnames(Stock))

## We reproduce below Figure 8(a) of Mercadier and Roustant (2019).

graphsEmp(Stock, k = 26, which = "taildependograph", names = colnames(Stock), select = 9)

## We reproduce below Figure 8(b) of Mercadier and Roustant (2019).

graphsEmp(Stock, k = 26, which = "taildependograph", names = colnames(Stock), select = 20)
}
\keyword{datasets}
\author{
Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
}
