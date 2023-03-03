\name{France}
\alias{France}
\docType{data}
\title{Dataset. Yearly Maxima of Temperature and coordinates of 21 French cities 1946-2000.}
\description{
  The \code{France} dataset is a list of two elements
  \itemize{
  \item{$ymt} a data frame of 55 rows and 21 columns, constructed after extraction from www.ecad.eu. The value at row \eqn{i} and colum \eqn{j} is the yearly maximum of temperature for the year \eqn{1946+i-1} in the \eqn{j-th} French city.
  \item{$coord} a list of two elements: Latitude \code{$lat} and Longitude \code{$lon} of 21 French cities.
  }

The name of the rows of \code{$ymt} are the year of the study: 1946--2000. The column names of \code{$ymt} are those of 21 French cities listed below.

[1] "MARSEILLE OBS. PALAIS-LONCHAMP"
 [2] "BOURGES AERODROME"
 [3] "BLAGNAC AEROP. TOULOUSE-BLAGNAC"
 [4] "MERIGNAC AEROPORT DE BORDEAUX"
 [5] "DEOLS CHATEAUROUX AERODROME DE DEOLS"
 [6] "PERPIGNAN"
 [7] "BRON LYON AEROPORT"
 [8] "PARIS-14E PARC MONTSOURIS"
 [9] "RENNES"
[10] "STRASBOURG-ENTZHEIM"
[11] "NANCY"
[12] "ORLEANS"
[13] "BESANCON"
[14] "LA-ROCHELLE"
[15] "BEAUVAIS-TILLE"
[16] "LE MANS"
[17] "METZ-FRESCATY"
[18] "MONTELIMAR"
[19] "NIMES"
[20] "VICHY-CHARMEIL"
[21] "COGNAC"



}

\references{

Klein Tank, A.M.G. and Coauthors, (2002).
Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment.
\emph{Int. J. of Climatol.}, \bold{22}, 1441--1453.
Data and metadata available at www.ecad.eu


Mercadier, C. and Roustant, O. (2019)
The tail dependograph.
\emph{Extremes}, \bold{22}, 343--372.


}
\seealso{
\code{\link[satdad]{tsicEmp}}, \code{\link[satdad]{ecEmp}},  \code{\link[satdad]{graphsMapEmp}}, \code{\link[satdad]{supportAnalysisEmp}}

}

\examples{
data(France)
maps::map('france',col='gray')
points(France$coord$lon,France$coord$lat, pch = 20, col = 1)
text(France$coord$lon,France$coord$lat+0.3,labels=1:21,cex=.8)
}
\keyword{datasets}
\author{
Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
}
