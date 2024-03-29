% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/satdad_Rfunctions.R
\name{tsic}
\alias{tsic}
\title{Tail superset importance coefficients for Mevlog models.}
\usage{
tsic(ds, ind = 2, n.MC = 1000, sobol = FALSE, norm = FALSE)
}
\arguments{
\item{ds}{An object of class \code{ds}.}

\item{ind}{A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = 2}, all pairwise coefficients are computed.}

\item{n.MC}{Monte Carlo sample size. Default value is 1000. See Details.}

\item{sobol}{A boolean. `FALSE` (the default). If `TRUE`:  the index is normalized by the theoretical global variance.}

\item{norm}{A boolean. `FALSE` (the default): original tsic is computed. `TRUE`:  tsic is normalized by its upper bound.}
}
\value{
The function returns a list of two elements
\itemize{
\item{\code{subsets}} A list of subsets from  \eqn{\{1,...,d\}}.

When \code{ind} is given as an integer, \code{subsets} is the list of subsets from  \eqn{\{1,...,d\}} with cardinality \code{ind}.

When \code{ind} is a list, it corresponds to \code{subsets}.

When \code{ind = "with.singletons"}  subsets is the list of all non empty subsets in \eqn{\{1,...,d\}}.

When \code{ind = "all"}   subsets is the list of all subsets in \eqn{\{1,...,d\}} with cardinality larger or equal to 2.

\item{\code{tsic}} A vector of tail superset importance coefficients associated with the list \code{subsets}. When \code{norm = TRUE}, then tsic  are normalized in the sense that the original values are divided by corresponding upper bounds.

}
}
\description{
Tail superset importance coefficients  for Mevlog models.  A \code{Mevlog} model is a multivariate extreme value (symmetric or asymmetric) logistic model.
}
\details{
The tail dependence structure is specified using a \code{ds} object, which corresponds to the stable tail dependence function  \eqn{\ell}.
The process for deducing  the stable tail dependence function \eqn{\ell} from \code{ds} is explained in the Details section of \code{\link[satdad]{gen.ds}}.

A tail superset importance coefficient (tsic) is a measure of the importance of a subset of components (and their supersets) in contributing to the global variance decomposition of  \eqn{\ell}.
The tsic  is computed using Monte Carlo methods based on the integral formula (3) in Mercadier and Roustant (2019).
Recall that Formula (9) in Liu and Owen (2006) provides an integral representation of the superset importance coefficient.

The tail dependograph is plotted using pairwise tsic values, which are computed using the function \code{tsic} and the \code{ind = 2} option.

The upper bound for a tsic associated with subset  \eqn{I} is given by Theorem 2 in Mercadier and Ressel (2021).
If  \eqn{|I|} is the cardinality of subset  \eqn{I}, then the upper bound is  \eqn{2 (|I| !)^2}/\eqn{((2|I|+2)!)}.

The tail dependence structure is set by a \code{ds} object. It thus corresponds to the stable tail dependence function \eqn{\ell}.
}
\examples{

## Fix a 5-dimensional asymmetric tail dependence structure
ds5 <- gen.ds(d = 5)

## Compute pairwise tsic
tsic(ds = ds5, ind = 2)

## Plot the tail dependograph
graphs(ds = ds5)

## Compute tsic on two specific subsets
tsic(ds = ds5, ind = list(1:4, 3:5))

## Compute normalized version of tsic
tsic(ds5,  ind = list(1:4, 3:5), norm = TRUE)

## Compute Sobol and normalized version of tsic
tsic(ds5,  ind = list(1:4, 3:5), norm = TRUE, sobol = TRUE)

}
\references{
Liu, R. and Owen, A. B. (2006)
Estimating mean dimensionality of analysis of variance decompositions.
\emph{J. Amer. Statist. Assoc.}, \bold{101(474)}:712--721.

Mercadier, C. and Ressel, P. (2021)
Hoeffding–Sobol decomposition of homogeneous co-survival functions: from Choquet representation to extreme value theory application.
Dependence Modeling, \bold{9(1)}, 179--198.

Mercadier, C. and Roustant, O. (2019)
The tail dependograph.
\emph{Extremes}, \bold{22}, 343--372.

Smith, R. L. (1990)
Max-stable processes and spatial extremes.
\emph{Dept. of Math., Univ. of Surrey}, Guildford GU2 5XH, England.

Tiago de Oliveira, J. (1962/63)
Structure theory of bivariate extremes, extensions.
\emph{Estudos de Matematica, Estatistica, e Economicos}, 7:165--195.
}
\seealso{
\code{\link[satdad]{graphs}}, \code{\link[satdad]{ellMevlog}}
}
\author{
Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
}
