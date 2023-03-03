## ----install, eval=FALSE------------------------------------------------------
#  install.packages("satdad")

## -----------------------------------------------------------------------------
library(satdad)

## -----------------------------------------------------------------------------
## Construction of a ds object without using gen.ds
ds5 <- vector("list")
ds5$d <- 5
ds5$type <- "alog" 
ds5$sub <- list(c(1,3),2:4,c(2,5))
ds5$asy <- list(c(1,.3),c(.5,1-.3,1), c(1-.5,1))
ds5$dep <- c(.2,.5,.3)

## -----------------------------------------------------------------------------
## Three constructions of ds object by using gen.ds
# only d is  given, sub, asy and dep are randomly sampled
ds10 <- gen.ds(d = 10)
# d and sub are given, asy and dep are randomly sampled
ds10 <- gen.ds(d = 10,  sub = list(1:2,1:7,3:5,7:10)) 
# d is  given, mnns indicates the cardinality of non singleton subsets in B
# sub, asy and dep are randomly sampled
ds10 <- gen.ds(d = 10, mnns = 4)

## -----------------------------------------------------------------------------
ds3 <- gen.ds(d = 3, type = "log")
ds3$dep

## -----------------------------------------------------------------------------
ds3 <- gen.ds(d = 3, type = "log", dep = .3)

## -----------------------------------------------------------------------------
n <- 1000
sample.frechet <- rMevlog(n, ds5) # standard Frechet margins
loc <- runif(5)
scale <- runif(5, 1, 2)
shape <- runif(5, -1, 1)
mar.gev <- cbind(loc, scale, shape)
sample.gev <- rMevlog(n, ds5, mar = mar.gev) # GEV margins all distinct
sample.samegev <- rMevlog(n, ds5, mar = c(-1,0.1,1)) # Gumbel margins

## ---- eval = FALSE------------------------------------------------------------
#  x5 <- runif(5)
#  ellMevlog(x5, ds5)
#  pMevlog(x5, ds5) # cdf under standard Frechet margins
#  pMevlog(x5, ds5, mar = c(1,1,0)) # cdf under standard Gumbel margins
#  dMevlog(x5, ds5) # pdf under standard Frechet margins

## -----------------------------------------------------------------------------
n <- 1000
sample.ext <- rArchimaxMevlog(n, ds5, dist = "ext")
lambda <- runif(1, 1, 2)
sample.exp <- rArchimaxMevlog(n, ds5, dist = "exp", dist.param = lambda)
shape <- runif(1, 1, 2)
scale <- runif(1, 1, 2)
sample.gamma <- rArchimaxMevlog(n, ds5, dist = "gamma", dist.param = c(shape, scale))

## -----------------------------------------------------------------------------
x <- runif(5)
ellMevlog(x, ds5)
ellArchimaxMevlog(x, ds5)
copArchimaxMevlog(x, ds5, dist = "ext")
copArchimaxMevlog(x, ds5,  dist = "exp", dist.param = lambda)
copArchimaxMevlog(x, ds5, dist = "gamma", dist.param = c(shape, scale))

## -----------------------------------------------------------------------------
res.tsic5 <- tsic(ds5)
as.character(res.tsic5$subsets)
res.tsic5$tsic

## -----------------------------------------------------------------------------
oldpar <- par(mfrow=c(1,2))
graphs(ds10) # (left) the nodes are plotted on an invisible circle
graphs(ds10, random = TRUE) # (right) the position of the nodes are random
par(oldpar)

## -----------------------------------------------------------------------------
oldpar <- par(mfrow=c(1,2))
graphs(ds3) # (left) the symmetric structure
graphs(ds5) # (right) the asymmetric structure contructed "manualy"
par(oldpar)

## -----------------------------------------------------------------------------
plotClev(ds5)

## -----------------------------------------------------------------------------
res.tic5 <- tic(ds5, ind = "with.singletons", sobol = TRUE)
sobol5 <- res.tic5$tic # which sum should be 1

## -----------------------------------------------------------------------------
res.ec10 <- ec(ds10)
as.character(res.ec10$subsets)
res.ec10$ec

## -----------------------------------------------------------------------------
oldpar <- par(mfrow=c(1,2))
graphs(ds5, which = "iecgraph")
graphs(ds10, which = "iecgraph")
par(oldpar)

## -----------------------------------------------------------------------------
res.ecEmp <- ecEmp(sample.ext, ind = "with.singletons", k = 100)
res.tsicEmp <- tsicEmp(sample.exp, ind = "all", k = 100)
res.ticEmp <- ticEmp(sample.gamma, ind = 4, k = 100)

## -----------------------------------------------------------------------------
graphsEmp(sample.ext, k = 100)
plotClevEmp(sample.exp, ind = "all", k = 100)

## -----------------------------------------------------------------------------
ds5$sub

## -----------------------------------------------------------------------------
supportAnalysisEmp(sample.gamma, k = 100)

## -----------------------------------------------------------------------------
library(graphicalExtremes)
g <- igraph::graph_from_edgelist(danube$flow_edges)
loc <- as.matrix(danube$info[,c('PlotCoordX', 'PlotCoordY')])
plot(g, layout = loc, vertex.color ="white", vertex.label.color = "darkgrey")

## -----------------------------------------------------------------------------
dan <- danube$data_clustered
graphsEmp(dan,  k=50, layout = loc)

## -----------------------------------------------------------------------------
lon <- as.numeric(unlist(danube$info[,"Long"]))
lat <- as.numeric(unlist(danube$info[,"Lat"]))*2
coord.dan <- list(lat = lat, lon = lon)
graphsMapEmp(dan, region = NULL, coord = coord.dan, k = 50, eps = 0.1)

## -----------------------------------------------------------------------------
plotClevEmp(dan,  k = 50,  ind = 2, labels = FALSE)

## -----------------------------------------------------------------------------
graphsEmp(dan,  k=50, layout = loc, select = 50, simplify = TRUE)
graphsMapEmp(dan, region = NULL, coord = coord.dan, k = 50, select = 50, eps = 0.1)

## -----------------------------------------------------------------------------
 ## Figure 9 (a) of Mercadier  and Roustant (2019).
graphsMapEmp(sample = France$ymt, k = 55,
        coord = France$coord,  region = 'France', thick.td = 3, select = 9)

## Figure 9 (b) of Mercadier  and Roustant (2019).
graphsMapEmp(sample = France$ymt, k = 55,
        coord = France$coord,  region = 'France', thick.td = 3, select = 30)

## Figure 9 (c) of Mercadier  and Roustant (2019).
graphsMapEmp(sample = France$ymt, k = 55, 
      coord = France$coord,  region = 'France', thick.td = 3)

## Figure 7(a) of Mercadier and Roustant (2019).
graphsEmp(Stock, k = 26, names = colnames(Stock), random = TRUE)

## Figure 8(a) of Mercadier and Roustant (2019).
graphsEmp(Stock, k = 26, names = colnames(Stock), random = TRUE, select = 9)

## Figure 8(b) of Mercadier and Roustant (2019).
graphsEmp(Stock, k = 26, names = colnames(Stock), random = TRUE, select = 20)

