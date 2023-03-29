#' r-p-d-ell- functions for Mevlog models.
#'
#' Random vectors generation (\code{rMevlog}),  cumulative distribution function (\code{pMevlog}), probability density function (\code{dMevlog}), stable tail dependence function (\code{ellMevlog}) for Mevlog models. A \code{Mevlog} model is a multivariate extreme value (symmetric or asymmetric) logistic model.
#'
#' @aliases rMevlog
#' @aliases pMevlog
#' @aliases dMevlog
#' @aliases ellMevlog
#'
#' @usage
#'
#' rMevlog(n, ds, mar = c(1,1,1))
#' pMevlog(x, ds, mar = c(1,1,1))
#' dMevlog(x, ds, mar = c(1,1,1))
#' ellMevlog(x, ds)
#'
#' @param n The number of observations.
#'
#' @param ds An object of class \code{ds}.
#'
#' @param mar A vector of length 3 or a \code{(d times 3)} matrix.  See details.
#'
#' @param x A vector of size \code{d} or a matrix with \code{d} columns.
#'
#' @details
#'
#' The tail dependence structure is set by a \code{ds} object.  See Section \bold{Value} in \code{\link[satdad]{gen.ds}}.
#'
#' The marginal information \code{mar} is given by a 3-dimensional vector (the order should be location, scale and shape) or a matrix with 3 columns depending on whether the components share the same characteristics or not.   When the marginal parameters differ, \code{mar} is a matrix containing \eqn{d} locations in the first column, \eqn{d} scales in the second column and \eqn{d} shapes in the third column.
#'
#' The (a)symmetric logistic models respectively are simulated in `rMevlog` using Algorithms 2.1 and 2.2 in Stephenson(2003).
#'
#' @returns
#'  \code{rMevlog}  returns a \code{(n times d)}  matrix  containing \code{n} realizations of a \code{d}-variate Mevlog random vector with margins \code{mar} and tail dependence structure \code{ds}.
#'
#'  \code{pMevlog} returns a  scalar (when \code{x} is a numeric vector) or a vector (when \code{x} is a numeric matrix, in which case the evaluation is done across the rows). The margins are provided by \code{mar} and the tail dependence structure through a \code{ds} object.
#'
#'  \code{dMevlog} returns a  scalar (when \code{x} is a numeric vector) or a vector (when \code{x} is a numeric matrix, in which case the evaluation is done across the rows). The margins are provided by \code{mar} and the tail  dependence structure through a \code{ds} object.
#'
#'  \code{ellMevlog} returns a  scalar (when \code{x} is a numeric vector) or a vector (when \code{x} is a numeric matrix, in which case the evaluation is done across the rows). The tail dependence structure is provided by a \code{ds} object.
#'
#' @references
#'
#' Gumbel, E. J. (1960)
#' Distributions des valeurs extremes en plusieurs dimensions.
#' \emph{Publ. Inst. Statist. Univ. Paris}, \bold{9}, 171--173.
#'
#' Stephenson, A. (2002)
#' evd: Extreme Value Distributions.
#' \emph{R News}, 2(2):31--32.
#'
#' Stephenson, A. (2003)
#' Simulating Multivariate Extreme Value Distributions of Logistic Type.
#' \emph{Extremes}, \bold{6}, 49--59.
#'
#' Tawn, J. A. (1990)
#' Modelling multivariate extreme value distributions.
#' \emph{Biometrika}, \bold{77}, 245--253.
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{gen.ds}}, \code{\link[satdad]{tsic}},  \code{\link[satdad]{ec}}, \code{\link[satdad]{graphs}}
#'
#' @examples
#' ## Fix a 3-dimensional symmetric tail dependence structure
#' ds3 <- gen.ds(d = 3, type = "log")
#'
#' ## The dependence parameter is given by
#' ds3$dep
#'
#' ## Generate a 1000-sample of Mevlog random vectors associated with ds3
#' ## The margins are kept as standard  Frechet
#' sample3 <- rMevlog(n = 1000, ds = ds3)
#'
#' ## Fix a 10-dimensional asymmetric tail dependence structure
#' # The option \cdoe{mns = 4} produces a support involving subsets of cardinality 4 plus singletons.
#' ds10 <- gen.ds(d = 10,  mnns = 4)
#' ## Margins differ from one to another
#' mar10 <- matrix(runif(10*3), ncol = 3)
#'
#' ## Generate a 50-sample of Mevlog random vectors associated with ds10 and mar10
#' sample10 <- rMevlog(n = 50, ds = ds10, mar = mar10)
#'
#' ## Continuing with ds3 ; we compute other attributes
#' ## The cumulative distribution function
#' pMevlog(x = rep(1,3), ds = ds3)
#' # should be similar to :
#' # evd::pmvevd(q = rep(1,3), dep = ds3$dep, model = "log", d = 3, mar = c(1,1,1))
#' ## The probability density function:
#' dMevlog(x = rep(1,3), ds = ds3, mar = c(1.2,1,0.5))
#' # should be similar to :
#' # evd::dmvevd(x = rep(1,3), dep = ds3$dep, model = "log", d = 3, mar = c(1.2,1,0.5))
#' ## The stable tail dependence function:
#' ellMevlog(x = rep(1,3), ds = ds3)
#'
#' @export rMevlog
#' @export dMevlog
#' @export pMevlog
#' @export ellMevlog
#'
rMevlog <- function(n, ds, mar = c(1,1,1))
{
  d <-  ds$d
  sub <- ds$sub
  dep <- ds$dep
  asy <- ds$asy
  type <- ds$type
  if (type=="log"){sample <- .algo21_numMat_cpp(n,d,dep)}else{sample <- .algo22_numMat_cpp(n,d,sub,dep,asy)}
  if(class(mar)[1]=="numeric"){
    if(length(mar)!=3){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")
    }else{
      loc <- mar[1];sc <- mar[2];sh <- mar[3];
      if(sh == 0){sample <- loc + sc * log(sample)
      }else{sample <- loc + sc * (sample^sh-1)/sh}
      }
  }else{
    if(!is.matrix(mar)){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")
    }else{
      if((ncol(mar)!=3)|(nrow(mar)!=d)){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")
      }else{
        loc <- mar[,1]; sc <- mar[,2]; sh <- mar[,3];
        for(j in 1:d){
          if(sh[j] == 0){sample[,j] <- loc[j] + sc[j] * log(sample[,j])
          }else{sample[,j] <- loc[j] + sc[j] * (sample[,j]^sh[j]-1)/sh[j]}
        }
      }
    }
  }
  return(sample)
}



ellMevlog <- function(x,ds){
  d <- ds$d
  sub <- ds$sub
  dep <- ds$dep
  asy <- ds$asy
  if(is.matrix(x)){
    if(ncol(x) != d){stop("The entry x should be a ds$d-dimensional vector or a matrix with ds$d columns.")}
    res <- .ellmevlogm_cpp(x,d,sub,dep,asy)
  }else{
    if(!is.numeric(x)){stop("The entry x should be a ds$d-dimensional vector or a matrix with ds$d columns.")
    }else{
      if(length(x) != d){stop("The entry x should be a ds$d-dimensional vector or a matrix with ds$d columns.")
      }else{
        res <- .ellmevlogv_cpp(x,d,sub,dep,asy)
      }
    }
  }
  return(res)
}


pMevlog <- function(x, ds, mar = c(1,1,1)){
  d <- ds$d
  if(is.matrix(x)){
    if(ncol(x) != d){stop("The entry v should be a ds$d-dimensional vector or a matrix with ds$d columns.")}
  }else{
    if(!is.numeric(x)){stop("The entry v should be a ds$d-dimensional vector or a matrix with ds$d columns.")
    }else{
      if(length(x) != d){stop("The entry v should be a ds$d-dimensional vector or a matrix with ds$d columns.")
      }
    }
  }
  x <- matrix(x, ncol = d)
  if(class(mar)[1] == "numeric"){
    if(length(mar) != 3){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")
    }else{loc <- mar[1]; sc <- mar[2]; sh <- mar[3];
    if(sh == 0){
      z <- exp(-(x-loc)/sc)
    }else{
      z <- (1+sh/sc*(x-loc))^(-1/sh)
    }
    }
  }else{
    if(!is.matrix(mar)){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")
    }else{
      if((ncol(mar) != 3)|(nrow(mar) != d)){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")
      }else{
        loc <- mar[,1]; sc <- mar[,2]; sh <- mar[,3]; z <- x
        for(j in 1:d){
          if(sh[j] == 0){
            z[,j] <- exp(-(x[,j]-loc[j])/sc[j])
          }else{
            z[,j] <- (1+sh[j]/sc[j]*(x[,j]-loc[j]))^(-1/sh[j])
          }
        }
      }
    }
  }
  return(exp(-ellMevlog(z, ds)))
}


dMevlog <- function(x, ds, mar = c(1,1,1)){
  d <- ds$d
  sub <- ds$sub
  dep <- ds$dep
  asy <- ds$asy
  if((class(mar)[1] == "numeric")*(length(mar) != 3)){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")}
  if((!class(mar)[1] == "numeric")*(!is.matrix(mar))){stop("The entry mar should be a 3-dimensional vector or a matrix of size ds$d*3.")}
  if(is.matrix(x)){
    if(ncol(x) != d){stop("The entry x should be a ds$d-dimensional vector or a matrix with ds$d columns.")}
    if(class(mar)[1] == "numeric"){
      if((mar[1] == 1)*(mar[2] == 1)*(mar[3] == 1)){res <- pMevlog(x,ds,mar)*.ratio_diff_pmevlogm_cpp(x,d,sub,dep,asy)
      }else{res <- pMevlog(x,ds,mar)*.ratio_diff_pmevlogm_marv_cpp(x,d,sub,dep,asy,mar)}
    }else{
      res <- pMevlog(x,ds,mar)*.ratio_diff_pmevlogm_marm_cpp(x,d,sub,dep,asy,mar)
    }
  }else{
    if(!is.numeric(x)){stop("The entry x should be a ds$d-dimensional vector or a matrix with ds$d columns.")
    }else{
      if(length(x) != d){stop("The entry x should be a ds$d-dimensional vector or a matrix with ds$d columns.")
      }else{
        if(class(mar)[1] == "numeric"){
          if((mar[1] == 1)*(mar[2] == 1)*(mar[3] == 1)){res <- pMevlog(x,ds,mar)*.ratio_diff_pmevlogv_cpp(x,d,sub,dep,asy)
          }else{res <- pMevlog(x,ds,mar)*.ratio_diff_pmevlogv_marv_cpp(x,d,sub,dep,asy,mar)}
        }else{
          res <- pMevlog(x,ds,mar)*.ratio_diff_pmevlogv_marm_cpp(x,d,sub,dep,asy,mar)
        }
      }
    }
  }
  return(res)
}


#'  Generate a Mevlog tail dependence structure.
#'
#'
#' The function \code{gen.ds} creates (possibly randomly) a tail dependence structure for a multivariate extreme value logistic (Mevlog) model.
#'
#' @aliases gen.ds
#'
#'
#' @usage
#'
#' gen.ds(d, type = "alog", sub = NULL, dep = NULL, asy = NULL, mnns = d)
#'
#' @param d The dimension.
#' @param type The type of the model; represented by a character string.  This is similar to the  option \code{model} of \code{\link[evd]{rmvevd}}. It must be either \code{"log"} or \code{"alog"}  (the default), for the symmetric logistic and the asymmetric logistic model respectively.
#' @param sub An optional list of subsets of \eqn{\{1,...,d\}} involved in the tail dependence structure. If \code{type = "log"}, then \code{sub} should be given by \eqn{(1,\ldots,d)}, which is the way the code \code{NULL} will be interpreted. If \code{type = "alog"} and \code{sub = NULL} then  a random  list of vectors,  subsets of  \eqn{\{1,...,d\}}, is created. The cardinality of non singleton subsets in \code{sub} is given by \code{mnns}.  If the user provides \code{sub}, it has to be a list of vectors,  subsets of \eqn{\{1,...,d\}}, where each component from \eqn{\{1,...,d\}} appears at least once; Otherwise, one should add the missing singleton(s).
#' @param dep An optional vector of dependence parameter(s). If \code{type = "log"}, \code{dep} should be a single value. Otherwise, if \code{type = "alog"} and if the list \code{sub} is provided, then the length of the vector \code{dep} should be equal to that of the list \code{sub} (or a single value that will be replicated the length of \code{sub} times). Among these values, the dependence parameters associated singletons have to be equal to one. Otherwise, the values of \code{dep} associated to singleton will be ignored (and set to one).  When  \code{dep = NULL} its values are randomly generated.
#' @param asy An optional list of asymmetric weights. If \code{type = "log"}, then \code{asy} should be the vector \eqn{(1,\ldots,1)}, which is the way the code \code{NULL} will be interpreted. If \code{type = "alog"} and if \code{sub} is provided, the length of the list \code{asy} should be in accordance with the length of \code{sub}.   If \code{asy = NULL} then the values are randomly generated. Note that \code{asy} satisfies the sum-to-one constraints.
#' @param mnns The default value is arbitrarily equal to \eqn{d}. When \code{sub = NULL}, the list \code{sub} is randomly generated, and its size is closely related to \code{mnns}. The latter represents the number of non singletons subsets included in \code{sub}.
#'
#' @details
#'
#' A multivariate extreme value logistic (Mevlog) model is symmetric or asymmetric.
#'    \itemize{
#'  \item{\code{type = "log". }} {It generates a multivariate symmetric logistic model. Such model is a well-known generalization of the bivariate extreme value logistic model introduced by Gumbel (1960). The parameter `dep`  (with \eqn{0 < `dep` \leq 1}) is the only parameter needed to write the following equation
#' \deqn{\ell(u) = ( \sum_{i=1}^d u_i^{1/\code{dep}} )^{\code{dep}}.}

#'  If the parameter \code{dep} is missing, the function \code{gen.ds}  will randomly generate its value from a standard uniform distribution.
#'  The list \code{asy} is reduced to a vector of ones whereas the list \code{sub} only contains the maximal vector \eqn{(1, \ldots, d)}.
#'
#' This is a special case of the multivariate asymmetric logistic model (\code{alog} case).}
#'
#' \item{\code{type = "alog". }}{It generates a  multivariate asymmetric logistic model, which has been first introduced by Tawn (1990). We have
#'   \deqn{\ell(u)=\sum_{b\in B}  (\sum_{i \in b} (\beta_{i,b}u_i)^{1/\alpha_b})^{\alpha_b}}
#'   where \eqn{B} is the power set of \eqn{\{1,...,d\}} (or a strict subset of the power set), the dependence parameters \eqn{\alpha_b} lie in \eqn{(0,1]} and the collection of asymmetric weights \eqn{\beta_{i,b}} are coefficients from [0,1] satisfying \eqn{\forall i \in \{1,\ldots,d\}, \sum_{b\in B: i \in b} \beta_{i,b}=1}.
#'   Missing asymmetric weights \eqn{\beta_{i,b}} are assumed to be zero. }
#'
#'   }
#'
#' The function \code{gen.ds}  generates here an object of class \code{ds} which corresponds in this package to the stable tail dependence function \eqn{\ell}. The class \code{ds} consists of:
#'
#'   \itemize{
#'   \item{the dimension \code{d}.}
#'   \item{the type \code{"log"} or \code{alog}.}
#'     \item{the list \code{sub} that corresponds to \eqn{B}.
#'     When \code{sub} is provided, the same list of subsets is returned, eventually sorted. When \code{sub = NULL} then \code{sub} is a list of subsets of the power set of \eqn{\{1,...,d\}}. When the option \code{mnns} is used, the latter integer indicates the cardinality of non singleton subsets in \eqn{B}.
#'    \item{the dependence parameter \code{dep} or the vector of dependence parameters \code{dep}. When missing, these coefficients are obtained from independent standard uniform sampling.
#'    }
#'    \item{the list \code{asy} of asymmetric weights \eqn{\beta_{i,b}} for \eqn{b \in B}  and \eqn{i \in b}. When missing, these coefficients are obtained from independent standard uniform sampling followed by renormalization in order to satisfy the sum-to-one constraints.}
#'    }
#'}
#'
#'
#' @return
#'
#' \code{gen.ds} returns an object representing a tail dependence structure for Mevlog models.
#' Such object is a list containing  the following components:
#'  \itemize{
#'  \item{\code{d}} The dimension.
#'
#'  \item{\code{type}} The type of the model either \code{"log"} or \code{"alog"}.
#'
#'  \item{\code{sub}} The list of subsets of \eqn{\{1,...,d\}} involved in the tail dependence support.
#'
#'  \item{\code{dep}}  The vector of dependence parameter(s).
#'
#'  \item{\code{asy}} The list of asymmetric weights.
#'}
#'
#'
#'@note
#' The first interest of the \code{gen.ds} function is to generate randomly a tail dependence structure. Since \code{sub} and \code{asy} become quickly very large lists as \eqn{d} increases, it is very convenient to obtain automatically well-defined tail dependence structures for  multivariate extreme value logistic models.
#'
#' The second interest of the \code{gen.ds} function is to produce partial models where all subsets do not necessarily contribute to the tail dependence support.
#'
#' The function \code{gen.ds} does not manage margins characteristics which will be handle by the option \code{mar} in the  \code{r-d-p-Mevlog} functions.
#'
#' @references
#' Gumbel, E. J. (1960)
#' Distributions des valeurs extremes en plusieurs dimensions.
#' \emph{Publ. Inst. Statist. Univ. Paris}, \bold{9}, 171--173.
#'
#' Stephenson, A. (2002)
#' evd: Extreme Value Distributions.
#' \emph{R News}, 2(2):31--32.
#'
#' Tawn, J. A. (1990)
#' Modelling multivariate extreme value distributions.
#' \emph{Biometrika}, \bold{77}, 245--253.
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#'   \code{\link[satdad]{ellMevlog}}, \code{\link[satdad]{graphs}}
#'
#' @examples
#'
#' ## Fix a 5-dimensional symmetric tail dependence structure
#' ## The dependence paramater is fixed to .7
#' (ds5 <- gen.ds(d = 5, dep = .7, type = "log"))
#'
#' ## Fix a 3-dimensional asymmetric tail dependence structure
#' ## The list sub and asy are provided ; The vector dep is randomly generated
#' (ds3 <- gen.ds(d = 3, sub = list(c(1,2), c(1,2,3)), asy = list(c(0.4,0.6), c(0.6,0.4,1))))
#' graphs(ds = ds3)
#'
#' ## Fix a 8-dimensional asymmetric tail dependence structure
#' ## The lists sub and asy, as the vector dep, are randomly generated
#' (ds8 <- gen.ds(d = 8))
#' graphs(ds = ds8)
#'
#' @export gen.ds
#'
#'
gen.ds <- function (d,  type = "alog", sub = NULL, dep = NULL, asy = NULL, mnns = d){
	if (type == "log"){
		if(!is.null(sub)){
			warn_sub <- "Since type = 'log', the list sub should be NULL or a list reduced to the vector 1:d."
			if(length(sub) != 1){stop(warn_sub)}
			if(length(sub[[1]]) != d){stop(warn_sub)}
			if(prod(unlist(sub) != 1:d)){stop(warn_sub)}
		}
		if(!is.null(asy)){
			warn_asy <- "Since type = 'log', the list asy should be NULL or a list reduced to the vector of ones."
			if(length(asy) != 1){stop(warn_asy)}
			if(length(asy[[1]]) != d){stop(warn_asy)}
			if(unlist(asy) != rep(1,d)){stop(warn_asy)}
		}
		warn_dep <- "Since type = 'log', the argument dep should be NULL or a scalar in (0,1)."
		if(!length(dep) %in% c(0,1)){stop(warn_dep)}
		if(!is.null(dep)){
			if((dep<0)|(dep>1)){stop(warn_dep)}
		}
		if(is.null(dep)){dep <- runif(1, min = 0.01)}
			sub <- list(1:d)
			asy <- list(rep(1,d))
		}else{
		if(is.null(sub)){
			if(!length(asy) %in% c(0,2^d-1)){stop("Since sub = NULL, the list asy should be NULL or a list of 2^d-1 elements")}
			if(!length(dep) %in% c(0,1,2^d-1)){stop("Since sub = NULL, the vector dep should be NULL or a scalar or a vector of 2^d-1 elements")}
		}else{
			if(!length(asy) %in% c(0,length(sub))){stop("The list asy should be NULL or should have the length of sub")}
			if(!length(dep) %in% c(0,1,length(sub))){stop("The list dep should be NULL or a scalar or a vector with length that of sub")}
			missind <- .find_missing_indices_cpp(d,sub)
			if(prod(missind)==0){missind=NA}
			if(!is.na(missind)){stop(paste("The list sub should be NULL or a list where all values from 1:d appears at least once \n  Here",missind,"is/are missing"))}
			if(sum(unlist(lapply(sub,function(v){is.unsorted(v)})))>0){stop(paste("Each vector of the list sub should be sorted in ascending order"))}
			if(!is.null(asy)){if(prod(unlist(lapply(sub,length)) != unlist(lapply(asy,length)))){stop("The lists sub and asy should have the same structure")}}
		}
		if(!is.null(asy)){
			if(sum(sapply(asy,sum)) != d){stop("The list asy should satisfy the sum-to-one constraints")}
			if(!prod(tapply(unlist(asy),factor(unlist(sub)),sum) == rep(1,d))){stop("The list asy should satisfy the sum-to-one constraints")}
		}
		if(!is.null(dep)){
			prodep <- prod(dep)
			if((prodep<0)|(prodep>1)){stop("The vector dep should be NULL or a vector with components in (0,1)")}
		}
		if(is.null(sub)){
			sub <- .subsets_cpp(d)
			sub <- sample(x=sample(sample(sub,length(sub)), length(sub)), size = mnns)
			missind <- .find_missing_indices_cpp(d, sub)
			if(prod(missind) == 0){sub <- sub}else{sub <- c(missind,sub)}
			sub <- .sort_sub_cpp(sub)
			if(is.null(asy)){
				asy <- .generate_asy_sub_cpp(d, sub)
			}
		if(is.null(dep)){
			ldep <- length(sub)
			dep <- runif(ldep, min = 0.01)
			ind_sing <- (unlist(lapply(sub,length)) == 1)
			dep[ind_sing] <- 1
			}
		}else{
			if(is.null(asy)){
				asy <- .generate_asy_sub_cpp(d, sub)
			}
		if(is.null(dep)){
			ldep <- length(sub)
			dep <- runif(ldep, min = 0.01)
			ind_sing <- (unlist(lapply(sub,length)) == 1)
			dep[ind_sing] <- 1
			}
		if(length(dep)==1){
			dep <- rep(dep, length(sub))
			}
		}
	}
	ds <- list(d = d, type = type, sub = sub, dep = dep, asy = asy)
	return(ds)
}


#' Extremal coefficients for  Mevlog models.
#'
#' Theoretical extremal coefficients for  Mevlog models. A \code{Mevlog} model is a multivariate extreme value (symmetric or asymmetric) logistic model.
#'
#' @param ds An object of class \code{ds}.
#'
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = 2}, all pairwise coefficients are computed.
#'
#' @param norm A boolean. `FALSE` (the default): ec is computed. `TRUE`:  inverse normalized ec is computed.
#'
#' @details
#'
#' The tail dependence structure is set by a \code{ds} object. It thus corresponds to the stable tail dependence function \eqn{\ell}. The way to deduce the stable tail dependence function \eqn{\ell} from \code{ds} is explained in the Details section of \code{\link[satdad]{gen.ds}}.
#'
#' @return
#' The function returns a list of two elements:
#' \itemize{
#' \item{\code{subsets}} A list of subsets from  \eqn{\{1,...,d\}}.
#'
#' When \code{ind} is given as an integer, \code{subsets} is the list of subsets from  \eqn{\{1,...,d\}} with cardinality \code{ind}. When \code{ind} is the list, it corresponds to \code{subsets}.
#'
#' When \code{ind = "with.singletons"}  subsets is the list of all non empty subsets in \eqn{\{1,...,d\}}.
#'
#' When \code{ind = "all"}   subsets is the list of all subsets in \eqn{\{1,...,d\}} with cardinality larger or equal to 2.
#'
#'
#' \item{\code{ec}} A vector of theoretical extremal coefficients associated with the list \code{subsets}.
#'
#' An extremal coefficient associated with the subset \eqn{I} is \eqn{\ell(1_I,0_{I^c})}. Its value lies in \eqn{(1, |I|)}.
#'
#' When \code{norm = TRUE}, then inverse normalized ec are computed by  \eqn{\dfrac{|I|-ec}{|I|-1}}.
#' }
#' @seealso
#'
#'  \code{\link[satdad]{ellMevlog}},  \code{\link[satdad]{gen.ds}}, \code{\link[satdad]{graphs}},  \code{\link[satdad]{tsic}}
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @export ec
#'
#' @references
#'
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#' Tiago de Oliveira, J. (1962/63)
#' Structure theory of bivariate extremes, extensions.
#' \emph{Estudos de Matematica, Estatistica, e Economicos}, 7:165--195.
#'
#' Smith, R. L. (1990)
#' Max-stable processes and spatial extremes.
#' \emph{Dept. of Math., Univ. of Surrey}, Guildford GU2 5XH, England.
#'
#' @examples
#'
#' ## Fix a 4-dimensional asymmetric tail dependence structure
#' ds4 <-  gen.ds(d = 4)
#' ## Compute all theoretical extremal coefficients
#' ec(ds = ds4, ind = "with.singletons")
#' ## Compute theoretical extremal coefficients associated with the support of ds4
#' ec(ds = ds4, ind = ds4$sub)
#'
#' ## Fix a 6-dimensional asymmetric tail dependence structure
#' ds6 <- gen.ds(d = 6, sub = list(1:2,2:5,5:6))
#' ## Compute all theoretical extremal coefficients on subsets with cardinality 5
#' ec(ds = ds6, ind = 5)
#' ## Compute inverse renormalized ec
#' ec(ds = ds6, ind = list(1:2,1:4,1:6), norm = TRUE)
#'
ec <- function(ds, ind = 2, norm = FALSE)
{
  d <- ds$d
  sub <- ds$sub
  dep <- ds$dep
  asy <- ds$asy
  type <- ds$type

  if(!is.list(ind)){
    if(ind == "all"){
      ind <-  .subsets_cpp(d)
      ind <- ind[-(1:d)]
    }else if(ind == "with.singletons"){
      ind <-  .subsets_cpp(d)
    }else{
      ind <- as.list(as.data.frame(combn(1:d, ind)))
      names(ind) <- NULL
    }
  }
  ind.sub <- ind
  ec <- .ecdsmevlog_cpp(d, sub, dep, asy, ind.sub)
  if(norm){u <- unlist(lapply(ind.sub, length)); ec <- (u-ec)/(u-1)}
  return(list(subsets = ind.sub, ec = ec))
}




#' Tail superset importance coefficients for Mevlog models.
#'
#' Tail superset importance coefficients  for Mevlog models.  A \code{Mevlog} model is a multivariate extreme value (symmetric or asymmetric) logistic model.
#'
#' @param ds An object of class \code{ds}.
#'
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = 2}, all pairwise coefficients are computed.
#'
#' @param n.MC Monte Carlo sample size. Default value is 1000. See Details.
#' @param sobol A boolean. `FALSE` (the default). If `TRUE`:  the index is normalized by the theoretical global variance.
#' @param norm A boolean. `FALSE` (the default): original tsic is computed. `TRUE`:  tsic is normalized by its upper bound.
#'
#' @details
#'
#' The tail dependence structure is specified using a \code{ds} object, which corresponds to the stable tail dependence function  \eqn{\ell}.
#' The process for deducing  the stable tail dependence function \eqn{\ell} from \code{ds} is explained in the Details section of \code{\link[satdad]{gen.ds}}.
#'
#' A tail superset importance coefficient (tsic) is a measure of the importance of a subset of components (and their supersets) in contributing to the global variance decomposition of  \eqn{\ell}.
#' The tsic  is computed using Monte Carlo methods based on the integral formula (3) in Mercadier and Roustant (2019).
#' Recall that Formula (9) in Liu and Owen (2006) provides an integral representation of the superset importance coefficient.
#'
#' The tail dependograph is plotted using pairwise tsic values, which are computed using the function \code{tsic} and the \code{ind = 2} option.
#'
#' The upper bound for a tsic associated with subset  \eqn{I} is given by Theorem 2 in Mercadier and Ressel (2021).
#' If  \eqn{|I|} is the cardinality of subset  \eqn{I}, then the upper bound is  \eqn{2 (|I| !)^2}/\eqn{((2|I|+2)!)}.
#'
#' The tail dependence structure is set by a \code{ds} object. It thus corresponds to the stable tail dependence function \eqn{\ell}.
#'
#'
#' @returns  The function returns a list of two elements
#' \itemize{
#' \item{\code{subsets}} A list of subsets from  \eqn{\{1,...,d\}}.
#'
#' When \code{ind} is given as an integer, \code{subsets} is the list of subsets from  \eqn{\{1,...,d\}} with cardinality \code{ind}.
#'
#' When \code{ind} is a list, it corresponds to \code{subsets}.
#'
#' When \code{ind = "with.singletons"}  subsets is the list of all non empty subsets in \eqn{\{1,...,d\}}.
#'
#' When \code{ind = "all"}   subsets is the list of all subsets in \eqn{\{1,...,d\}} with cardinality larger or equal to 2.
#'
#' \item{\code{tsic}} A vector of tail superset importance coefficients associated with the list \code{subsets}. When \code{norm = TRUE}, then tsic  are normalized in the sense that the original values are divided by corresponding upper bounds.
#'
#' }
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{graphs}}, \code{\link[satdad]{ellMevlog}}
#'
#' @export tsic
#'
#' @references
#'
#' Liu, R. and Owen, A. B. (2006)
#' Estimating mean dimensionality of analysis of variance decompositions.
#' \emph{J. Amer. Statist. Assoc.}, \bold{101(474)}:712--721.
#'
#' Mercadier, C. and Ressel, P. (2021)
#' Hoeffding–Sobol decomposition of homogeneous co-survival functions: from Choquet representation to extreme value theory application.
#' Dependence Modeling, \bold{9(1)}, 179--198.
#'
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#' Smith, R. L. (1990)
#' Max-stable processes and spatial extremes.
#' \emph{Dept. of Math., Univ. of Surrey}, Guildford GU2 5XH, England.
#'
#' Tiago de Oliveira, J. (1962/63)
#' Structure theory of bivariate extremes, extensions.
#' \emph{Estudos de Matematica, Estatistica, e Economicos}, 7:165--195.
#'
#' @examples
#'
#' ## Fix a 5-dimensional asymmetric tail dependence structure
#' ds5 <- gen.ds(d = 5)
#'
#' ## Compute pairwise tsic
#' tsic(ds = ds5, ind = 2)
#'
#' ## Plot the tail dependograph
#' graphs(ds = ds5)
#'
#' ## Compute tsic on two specific subsets
#' tsic(ds = ds5, ind = list(1:4, 3:5))
#'
#' ## Compute normalized version of tsic
#' tsic(ds5,  ind = list(1:4, 3:5), norm = TRUE)
#'
#' ## Compute Sobol and normalized version of tsic
#' tsic(ds5,  ind = list(1:4, 3:5), norm = TRUE, sobol = TRUE)
#'
tsic <- function(ds, ind = 2, n.MC = 1000, sobol = FALSE, norm = FALSE)
{
	d <- ds$d
	sub <- ds$sub
	dep <- ds$dep
	asy <- ds$asy
	type <- ds$type

	if(!is.list(ind)){
	  if(ind == "all"){
	    ind <-  .subsets_cpp(d)
	    ind <- ind[-(1:d)]
	  }else if(ind == "with.singletons"){
	    ind <-  .subsets_cpp(d)
	  }else{
	    ind <- as.list(as.data.frame(combn(1:d, ind)))
	    names(ind) <- NULL
	  }
	}
	ind.sub <- ind
	tsic <- .tsicdsmevlog_list_cpp(n.MC, d, sub, dep, asy, ind.sub)
	if(sobol){va <- .tsicdsmevlog_empty_cpp(n.MC, d, sub, dep, asy); tsic <- tsic/va}
	if(norm){u <- unlist(lapply(ind.sub,length));
	dimu <- 2*(factorial(u))^2/factorial(2*u+2)
	inv.b <- 1/dimu
	tsic <- tsic*inv.b}
	return(list(subsets = ind.sub, tsic = tsic))
}


#' Graphs of the tail dependence structure for Mevlog models.
#'
#' Tail dependograph and Inverse extremal coefficients graph  for Mevlog models.   A \code{Mevlog} model is a multivariate extreme value (symmetric or asymmetric) logistic model.
#'
#' @name graphs
#'
#' @param ds An object of class \code{ds}.
#'
#' @param names A character vector of length \code{d} which replaces \code{as.character(1:d)} (the default ones).
#'
#' @param n.MC Monte Carlo sample size. Default value is 1000. See details in  \code{\link[satdad]{tsic}}.
#'
#' @param which A character string:  \code{taildependograph} (the default), \code{iecgraph}, or \code{both},
#'
#' @param random A boolean. `FALSE` (the default): the vertex positions are fixed along a circle.  `TRUE`: some randomness is applied for positioning the vertices.
#'
#' @param thick.td A numeric value for the maximal thickness of edges in \code{taildependograph}. Default value is 5.
#'
#' @param thick.ec A numeric value for the maximal thickness of edges in \code{iecgraph}. Default value is 5.
#'
#' @details
#'
#' The tail dependence structure is set by a \code{ds} object. It thus corresponds to the stable tail dependence function \eqn{\ell}. The way to deduce the stable tail dependence function \eqn{\ell} from \code{ds} is explained in the Details section of \code{\link[satdad]{gen.ds}}.
#'
#'
#' @returns
#'
#' The function returns either the tail dependograph or the inverse extremal coefficients graph, or both, for the tail dependence structure `ds`.
#'
#' The tail dependograph displays pairwise tail superset importance coefficients, which measure the extent to which pairs of components (and their supersets) contribute to the overall variance of the stable tail dependence function.
#' We refer to  Mercadier, C. and Roustant, O. (2019) for more details. These coefficients are computed using the `tsic` function with the `"ind = 2"` option.
#'
#' The inverse extremal coefficients graph shows the inverse renormalized pairwise coefficients computed as \eqn{\theta_{ij}=1-\ell(1_i,1_j,\bold{0})/2}.
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{tsic}}, \code{\link[satdad]{ec}}, \code{\link[satdad]{ellMevlog}}
#'
#' @references
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#' Tiago de Oliveira, J. (1962/63)
#' Structure theory of bivariate extremes, extensions.
#' \emph{Estudos de Matematica, Estatistica, e Economicos}, 7:165--195.
#'
#' Smith, R. L. (1990)
#' Max-stable processes and spatial extremes.
#' \emph{Dept. of Math., Univ. of Surrey}, Guildford GU2 5XH, England.
#'
#' @examples
#'
#' ## Fix a 8-dimensional asymmetric tail dependence structure
#' ds8 <- gen.ds(d = 8)
#'
#' ## Plot the graphs that illustrate  characteristics of the tail dependence structure
#' graphs(ds = ds8, which = "both")
#'
#' @export graphs

graphs <- function(ds, names = NULL, n.MC = 1000, which = "taildependograph", random = FALSE, thick.td = 5, thick.ec = 5){

  if(!((which == "both")|(which == "taildependograph")|(which=="iecgraph"))){stop("The entry which should be a character string: 'both', 'taildependograph' or 'iecgraph'.")}

  d <- ds$d

  if(!is.null(names)){noms <- names}else{noms <- 1:d}

	E <- t(combn(d, 2))
	g <- igraph::graph(as.vector(t(E)), n = d, directed = FALSE)
	if(random){layout <- igraph::layout.fruchterman.reingold(g)
	}else{layout <- matrix(c(cos((2*pi)/d*(0:(d-1))), sin((2*pi)/d*(0:(d-1)))), ncol = 2)}

	if((which == "both")|(which == "taildependograph")){
	res <- tsic(ds = ds, ind = 2, n.MC = n.MC, sobol = TRUE, norm = TRUE)
	igraph::plot.igraph(g,layout = layout, edge.width = res$tsic/max(res$tsic) * thick.td, vertex.frame.color="darkgrey", vertex.color ="white",
												vertex.label.color = "darkgrey", vertex.label = noms, main="Tail Dependograph")
	mtext(paste("n.MC =",n.MC,"    ","thick.td = ",thick.td),side=1,col="gray")
	}
	if((which=="both")|(which=="iecgraph")){
	res_ec <- ec(ds=ds, ind = 2)
	TWOMINUSec <- 2-res_ec$ec
	igraph::plot.igraph(g,layout = layout, edge.width = TWOMINUSec * thick.ec, vertex.frame.color="darkgrey", vertex.color ="white",
	                    vertex.label.color = "darkgrey", vertex.label = noms, main="Inverse Extremal Coeff. Graph")
	mtext(paste("thick.ec = ",thick.ec), side = 1, col = "gray")
	}

}





#' Empirical stable tail dependence function.
#'
#' The stable tail dependence function of \code{sample} is estimated at each row of \code{x} and for all values of the threshold parameter \code{k}.
#'
#' @name ellEmp
#'
#' @param sample A \code{(n times d)} matrix.
#'
#' @param x A \code{(N.x times d)} matrix.
#'
#' @param k A vector of \code{N.k} integers smaller or equal to \code{n}.
#'
#'
#' @returns
#'
#' A \code{(N.k times N.x)} matrix is returned.
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{ellMevlog}}, \code{\link[satdad]{gen.ds}}
#'
#' @references
#'
#' Huang, X. (1992).
#' Statistics of bivariate extremes.
#' PhD Thesis, Erasmus University Rotterdam, Tinbergen Institute Research series No. 22.
#'
#' de Haan, L. and Resnick, S. I. (1993).
#' Estimating the limit distribution of multivariate extremes. Communications in Statistics.
#' Stochastic Models 9, 275--309.
#'
#' Fougeres, A.-L., de Haan, L. and  Mercadier, C.  (2015).
#' Bias correction in multivariate extremes.
#' Annals of Statistics 43 (2), 903--934.
#'
#' @examples
#'
#' ## Fix a 5-dimensional asymmetric tail dependence structure
#' ds5 <- gen.ds(d = 5)
#'
#' ## Construct a 1000-sample of Mevlog random vector associated with ds5
#' sample5 <- rMevlog(n = 1000, ds = ds5)
#'
#' ## Select 3 vectors in R^5
#' x5 <- matrix(runif(5*3), ncol = 5)
#'
#' ## Select 4 values for the threshold parameter
#' k5 <- (2:5)*10
#'
#' ## Estimation of the stable tail dependence function
#' # We thus get a 4 x 3 matrix
#' ellEmp(sample = sample5, x = x5, k = k5)
#'
#' ## Theoretical values of the stable tail dependence function inherited from ds5
#' ellMevlog(x = x5, ds = ds5)
#'
#' @export ellEmp

ellEmp <- function(sample, x, k){
    if(class(x)[1] != "matrix"){x <- matrix(x, ncol = ncol(sample), byrow = TRUE)}
    res <- .ellEmp_cpp(k, x, sample)
  rownames(res) <- paste("k", k, sep = ".")
  colnames(res) <- paste(paste("x[,", 1:nrow(x), sep = ""),"]",sep = "")
  if(length(k) == 1){res <- as.numeric(res)}
  return(res)
}





#' Empirical tail superset importance coefficients.
#'
#' Computes on a sample the tail superset importance coefficients (tsic) associated with threshold \code{k}. The value may be renormalized by the empirical global variance (Sobol version) and/or by its theoretical upper bound.
#'
#' @param sample A \code{(n times d)} matrix.
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = 2}, all pairwise coefficients are computed.
#' @param k An integer smaller or equal to \code{n}.
#' @param sobol A boolean. `FALSE` (the default). If `TRUE`:  the index is normalized by the empirical global variance.
#' @param norm A boolean. `FALSE` (the default). If `TRUE`: the index is normalized by its theoretical upper bound.
#'
#' @details
#' The theoretical functional decomposition of the variance of the stdf \eqn{\ell} consists in writing \eqn{D(\ell) = \sum_{I \subseteq \{1,...,d\}} D_I(\ell) } where \eqn{D_I(\ell)} measures the variance of \eqn{\ell_I(U_I)} the term associated with subset \eqn{I} in the Hoeffding-Sobol decomposition of \eqn{\ell}
#' ; note that \eqn{U_I} represents a random vector with independent standard uniform entries.
#'
#' Fixing a subset of components \eqn{I}, the theoretical tail superset importance coefficient is defined by \eqn{\Upsilon_I(\ell)=\sum_{J \supseteq I} D_J(\ell)}.
#' A theoretical upper bound for tsic \eqn{\Upsilon_I(\ell)} is given by Theorem 2 in Mercadier and Ressel (2021)
#' which states that \eqn{\Upsilon_I(\ell)\leq 2(|I|!)^2/((2|I|+2)!)}.
#'
#' Here, the function \code{tsicEmp} evaluates, on a \eqn{n}-sample and threshold \eqn{k},  the empirical tail superset  importance coefficient \eqn{\hat{\Upsilon}_{I,k,n}} the empirical counterpart of \eqn{\Upsilon_I(\ell)}.
#'
#' Under the option \code{sobol = TRUE}, the function \code{tsicEmp} returns  \eqn{\dfrac{\hat{\Upsilon}_{I,k,n}}{\hat{D}_{k,n}}} the empirical counterpart of \eqn{\dfrac{\Upsilon_I(\ell)}{D_I(\ell)}}.
#'
#' Under the option \code{norm = TRUE}, the quantities are multiplied by \eqn{\dfrac{(2|I|+2)!}{2(|I|!)^2}}.
#'
#' Proposition 1 and Theorem 2 of Mercadier and Roustant (2019) provide several rank-based expressions
#'
#' \eqn{\hat{\Upsilon}_{I,k,n}=\frac{1}{k^2}\sum_{s=1}^n\sum_{s^\prime=1}^n \prod_{t\in I}(\min(\overline{R}^{(t)}_s,\overline{R}^{(t)}_{s^\prime})-\overline{R}^{(t)}_{s}\overline{R}^{(t)}_{s^\prime}) \prod_{t\notin I} \min(\overline{R}^{(t)}_s,\overline{R}^{(t)}_{s^\prime})}
#'
#' \eqn{\hat{D}_{k,n}=\frac{1}{k^2}\sum_{s=1}^n\sum_{s^\prime=1}^n \prod_{t\in I}\min(\overline{R}^{(t)}_s,\overline{R}^{(t)}_{s^\prime})- \prod_{t\in I}\overline{R}^{(t)}_{s}\overline{R}^{(t)}_{s^\prime}}
#'
#' where
#'
#' \itemize{
#' \item{} \eqn{k} is the threshold parameter,
#' \item{} \eqn{n} is the sample size,
#' \item{} \eqn{X_1,...,X_n} describes the \code{sample}, each \eqn{X_s} is a d-dimensional vector \eqn{X_s^{(t)}} for \eqn{t=1,...,d},
#' \item{} \eqn{R^{(t)}_s} denotes the rank of \eqn{X^{(t)}_s} among \eqn{X^{(t)}_1, ..., X^{(t)}_n},
#' \item{} and  \eqn{\overline{R}^{(t)}_s = \min((n- R^{(t)}_s+1)/k,1)}.
#' }
#'
#'
#' @return
#' The function returns a list of two elements:
#'
#' \itemize{
#' \item{\code{subsets}} A list of subsets from  \eqn{\{1,...,d\}}.
#'
#' When \code{ind} is given as an integer, \code{subsets} is the list of subsets from  \eqn{\{1,...,d\}} with cardinality \code{ind}. When \code{ind} is the list, it corresponds to \code{subsets}.
#'
#' When \code{ind = "with.singletons"}  subsets is the list of all non empty subsets in \eqn{\{1,...,d\}}.
#'
#' When \code{ind = "all"}   subsets is the list of all subsets in \eqn{\{1,...,d\}} with cardinality larger or equal to 2.
#'
#' \item{\code{tsic}} A vector of empirical tail superset importance coefficients associated with the list \code{subsets}. When \code{norm = TRUE}, then tsic  are normalized in the sense that the original values are divided by corresponding upper bounds.
#' }
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{graphsEmp}}, \code{\link[satdad]{ellEmp}}
#'
#' @references
#'
#' Mercadier, C. and Ressel, P. (2021)
#' Hoeffding–Sobol decomposition of homogeneous co-survival functions: from Choquet representation to extreme value theory application.
#' Dependence Modeling, \bold{9(1)}, 179--198.
#'
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#' @export tsicEmp
#'
#' @examples
#'
#' ## Fix a 6-dimensional asymmetric tail dependence structure
#' ds <- gen.ds(d = 6, sub = list(1:4,5:6))
#'
#' ## Plot the  tail dependograph
#' graphs(ds)
#'
#' ## Generate a 1000-sample of Archimax Mevlog random vectors
#' ## associated with ds and underlying distribution exp
#' sample <- rArchimaxMevlog(n = 1000, ds = ds, dist = "exp", dist.param = 1.3)
#'
#' ## Compute tsic values associated with subsets
#' ## of cardinality 2 or more \code{ind = "all"}
#' res <- tsicEmp(sample = sample, ind = "all", k = 100, sobol = TRUE, norm = TRUE)
#'
#' ## Select the significative tsic
#' indices_nonzero <- which(res$tsic %in% boxplot.stats(res$tsic)$out == TRUE)
#'
#' ## Subsets associated with significative tsic reflecting the tail support
#' as.character(res$subsets[indices_nonzero])
#'
#' ## Pairwise tsic are obtained by
#' res_pairs <- tsicEmp(sample = sample, ind = 2, k = 100, sobol = TRUE, norm = TRUE)
#'
#' ## and plotted in the tail dependograph
#' graphsEmp(sample, k = 100)

tsicEmp <- function(sample, ind = 2, k, sobol = FALSE, norm = FALSE){
  if(!is.list(ind)){
    if(ind == "all"){
      ind <-  .subsets_cpp(ncol(sample))
      ind <- ind[-(1:ncol(sample))]
    }else if(ind == "with.singletons"){
      ind <-  .subsets_cpp(ncol(sample))
    }else{
      ind <- as.list(as.data.frame(combn(1:ncol(sample), ind)))
      names(ind) <- NULL
    }
  }
  Rbar <- .rankbar_cpp(sample, k)
  ratio <- .ratioall_cpp(Rbar)
  prodmin <- .prodmin_cpp(Rbar)
  if(sobol){
    prodprod <-  .prodprod_cpp(Rbar)
    D <- sum(prodmin-prodprod)/k^2
  }
  N = length(ind)
  upsi <- vector("list", N)
  for(m in 1:N){
    sub <- ind[[m]]
    aux <- 1
    for(i in sub){aux <- aux*ratio[[i]]}
    upsi[[m]] <- aux*prodmin
  }
  upsi <- as.numeric(lapply(upsi, sum))/k^2
  if(sobol){upsi <- upsi/D}
  if(norm){
    s <- unlist(lapply(ind, length))
    dims <- 2*(factorial(s))^2/factorial(2*s+2)
    b <- 1/dims
    upsi <- upsi*b}
  return(list(subsets = ind, tsic = upsi))
}

#' Empirical Extremal coefficients.
#'
#' Computes on a sample the extremal coefficients associated with threshold \code{k}.
#'
#'
#' @param sample A \code{(n times d)} matrix.
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = 2}, all pairwise coefficients are computed.
#' @param k An integer smaller or equal to \code{n}.
#' @param norm A boolean. `FALSE` (the default): empirical ec is computed. `TRUE`:  inverse normalized empirical ec is computed.
#'
#' @return
#' The function returns a list of two elements:
#' \itemize{
#' \item{\code{subsets}} A list of subsets from  \eqn{\{1,...,d\}}.
#'
#' When \code{ind} is given as an integer, \code{subsets} is the list of subsets from  \eqn{\{1,...,d\}} with cardinality \code{ind}. When \code{ind} is the list, it corresponds to \code{subsets}.
#'
#' When \code{ind = "with.singletons"}  subsets is the list of all non empty subsets in \eqn{\{1,...,d\}}.
#'
#' When \code{ind = "all"}   subsets is the list of all subsets in \eqn{\{1,...,d\}} with cardinality larger or equal to 2.
#'
#' \item{\code{ec}} A vector of empirical  extremal coefficients.
#'
#' An empirical extremal coefficient associated with the subset \eqn{I} is \eqn{\hat{\ell}_{k,n}(1_I,0_{I^c})}. Its value lies in \eqn{(1, |I|)}.
#'
#' When \code{norm = TRUE}, then inverse normalized empirical ec are computed by \eqn{1 - \dfrac{\hat{\ell}_{k,n}(1_I,0_{I^c})}{|I|}}.
#' }
#'
#' @export ecEmp
#'
#' @seealso
#' \code{\link[satdad]{ec}}, \code{\link[satdad]{ellEmp}},  \code{\link[satdad]{graphsEmp}}
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @examples
#' ## We produce below a figure on the dataset used in Mercadier  and Roustant (2019).
#'
#' data(France)
#'
#' ec_ymt <- ecEmp(sample = France$ymt, ind = 2, k = 25)
#'
#' ## The 9 largest inverse empirical pairwise extremal coefficients.
#'  graphsMapEmp(France$ymt, region='france', coord=France$coord, k=25, which="iecgraph", select=9)
#'
#' ## The 30 largest inverse empirical pairwise extremal coefficients.
#' graphsMapEmp(France$ymt, region='france', coord=France$coord, k=25, which="iecgraph", select=30)
#'
#' ## All the inverse empirical pairwise extremal coefficients.
#' graphsMapEmp(France$ymt, region='france', coord=France$coord, k=25, which="iecgraph")
#'
ecEmp <- function(sample, ind = 2, k, norm = TRUE){
  if(!is.list(ind)){
    if(ind == "all"){
      ind <-  .subsets_cpp(ncol(sample))
      ind <- ind[-(1:ncol(sample))]
    }else if(ind == "with.singletons"){
      ind <-  .subsets_cpp(ncol(sample))
    }else{
      ind <- as.list(as.data.frame(combn(1:ncol(sample), ind)))
      names(ind) <- NULL
    }
  }
  N = length(ind)
  empec <- rep(NA, N)
  for(m in 1:N){
    sub <- ind[[m]]
    x <- rep(0, ncol(sample))
    x[sub] <- 1
    empec[m] <- ellEmp(sample, x = x, k = k)
  }
  if(norm){u <- unlist(lapply(ind, length)); empec <- (u-empec)/(u-1)}
  return(list(subsets = ind, ec = empec))
}


#' Empirical graphs of the tail dependence structure.
#'
#' Empirical tail dependograph and  empirical inverse extremal coefficients graph of the tail dependence structure on a \code{sample} associated with threshold k.
#'
#' @name graphsEmp
#'
#' @param sample  A \code{(n times d)} matrix.
#'
#' @param layout The vertex coordinates as a \code{(d times 2)} matrix. The default is NULL. See also the parameter \code{random}.
#'
#' @param names A character vector of length \code{d} which replaces \code{as.character(1:d)} (the default ones).
#'
#' @param k An integer smaller or equal to \code{n}.
#'
#' @param which A character string: \code{taildependograph} (the default), \code{iecgraph} or  \code{both}.
#'
#' @param select If select = NULL (default) all edges are plotted. If select is an integer between 1 and the number of possible pairs of components of sample, then only the select largest edges are plotted.
#'
#' @param simplify If select is not NULL, and if a vertex is not associated with one of the selected edges, this vertex is not printed.
#'
#' @param random A boolean. `FALSE` (the default): the vertex positions are fixed along a circle when layout is NULL.  `TRUE`: some randomness is applied for positioning the vertices.
#'
#' @param thick.td A numeric value for the maximal thickness of edges in \code{taildependograph}. Default value is 5.
#'
#' @param thick.ec A numeric value for the maximal thickness of edges in \code{iecgraph}. Default value is 5.
#'
#' @returns
#'
#'  It returns both (or one among) the empirical tail dependograph and the empirical inverse extremal coefficients graph of the  \code{sample}.
#'
#'  The empirical tail dependograph represents the pairwise empirical  tail superset importance coefficients, see  Mercadier, C. and Roustant, O. (2019).
#'  These indices are computed by the function \code{tsicEmp}.
#'  It measures how much a pair of components (included supersets of this pair of components) is involved in the asymptotic dependence of the sample.
#'
#'  The empirical Inverse extremal coefficients graph represents empirical pairwise coefficients that estimate \eqn{1-\ell(1_i,1_j,\bold{0})/2}.
#'
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#' @seealso
#'
#' \code{\link[satdad]{ecEmp}}, \code{\link[satdad]{tsicEmp}}
#'
#' @references
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#' @examples
#'
#' ## Fix a 8-dimensional asymmetric tail dependence structure
#' ds8 <- gen.ds(d = 8)
#'
#' ## Generate a 200-sample of Frechet margins Mevlog model associated with ds8
#' sample8 <- rMevlog(n = 200 , ds = ds8)
#'
#' ## Plot the tail dependograph of ds8
#' graphs(ds = ds8)
#'
#' ## Its empirical version for k = 20
#' graphsEmp(sample = sample8, k = 20)
#'
#' ## Its empirical version for k = 20 restricted to the 3 largest edges
#' graphsEmp(sample = sample8, k = 20, select = 3)
#'
#' ## Plot the Inverse extremal coefficients graph of ds8
#' graphs(ds = ds8, which = "iecgraph")
#'
#' ## Its empirical version for k = 20
#' graphsEmp(sample = sample8, k = 20, which = "iecgraph")
#'
#' ## Its empirical version for k = 20 restricted to the 3 largest edges
#' graphsEmp(sample = sample8, k = 20, which = "iecgraph", select = 3)
#'
#' ## Plot the empirical tail dependograph
#' ## on river discharge data for tributaries
#' ## of the Danube extracted from
#' ##  Asadi P., Davison A.C., Engelke S. (2015).
#' ## “Extremes on river networks.”
#' ## The Annals of Applied Statistics, 9(4), 2023 – 2050.
#' #NOT RUN dan <- graphicalExtremes::danube$data_clustered
#' #NOT loc <- as.matrix(graphicalExtremes::danube$info[,c('PlotCoordX', 'PlotCoordY')])
#' #NOT graphsEmp(dan,  k=50, layout = loc)
#'
#' @export graphsEmp

graphsEmp <- function(sample, layout = NULL, names = NULL, k, which = "taildependograph", select = NULL, simplify = FALSE, random = FALSE, thick.td = 5, thick.ec = 5){
  if(!((which == "both")|(which == "taildependograph")|(which == "iecgraph"))){stop("The entry which should be a character string: 'both', 'taildependograph' or 'iecgraph'.")}

  d <- ncol(sample)

  if(!is.null(names)){noms <- names}else{noms <- 1:d}

  paires <- combn(1:d, 2)
  if(!is.null(select)){if(select > ncol(paires)){stop("The entry `select` should be smaller than the number of pairs.")}}

  E <- t(paires)
  g <- igraph::graph(as.vector(t(E)), n = d, directed = FALSE)
  if(is.null(layout)){
    if(random){layout <- igraph::layout.fruchterman.reingold(g)
    }else{layout <- matrix(c(cos((2*pi)/d*(0:(d-1))), sin((2*pi)/d*(0:(d-1)))), ncol = 2)
    }}

  if((which == "both")|(which == "taildependograph")){
    res <- tsicEmp(sample, ind = 2, k = k, sobol = TRUE, norm = TRUE)
    gp <- g
    present <- 1:d
    nomsp <- noms
    layoutp <- layout
    if(!is.null(select)){
      indices <- order(res$tsic, decreasing = TRUE)[1:select]
      res$tsic[-indices] <- 0
      if(simplify){present <- sort(unique(unlist(res$subsets[indices])))
      }else{present <- 1:d}
      dp <- length(present)
      samplep <- sample[,present]
      nomsp <- noms[present]
      pairesp <- combn(1:dp, 2)
      Ep <- t(pairesp)
      gp <- igraph::graph(as.vector(t(Ep)), n = dp, directed = FALSE)
      layoutp <- layoutp[present,]
      res <- tsicEmp(samplep, ind = 2, k = k, sobol = TRUE, norm = TRUE)
      indices <- order(res$tsic, decreasing = TRUE)[1:select]
      res$tsic[-indices] <- 0
    }
    igraph::plot.igraph(gp,layout = layoutp, edge.width = res$tsic/max(res$tsic)  * thick.td, vertex.frame.color = "darkgrey", vertex.color = "white",
                        vertex.label.color = "darkgrey", vertex.label = nomsp,  main = "Emp. Tail Dependograph")
    message <- paste("k =", k ,"    ","thick.td = ",thick.td)
    if(!is.null(select)){if(select != d*(d-1)/2){ message <- paste(message,"    ", "nb.selected.edges = ", select, "\n vertices.unprinted = ",deparse(noms[(1:d)[-present]],width.cutoff = 500))}}
    mtext(message, side = 1, col = "gray")
  }
  if((which == "both")|(which == "iecgraph")){
    res_ec <- ecEmp(sample, ind = 2, k = k, norm = TRUE)
    TWOMINUSec <- res_ec$ec
    gp <- g
    present <- 1:d
    nomsp <- noms
    layoutp <- layout
    if(!is.null(select)){
      indices <- order(TWOMINUSec, decreasing = TRUE)[1:select]
      TWOMINUSec[-indices] <- 0
      if(simplify){present <- sort(unique(unlist(res_ec$subsets[indices])))
      }else{present <- 1:d}
      dp <- length(present)
      samplep <- sample[,present]
      nomsp <- noms[present]
      pairesp <- combn(1:dp, 2)
      Ep <- t(pairesp)
      gp <- igraph::graph(as.vector(t(Ep)), n = dp, directed = FALSE)
      layoutp <- layoutp[present,]
      res_ec <- ecEmp(samplep, ind = 2, k = k, norm = TRUE)
      TWOMINUSec <- res_ec$ec
      indices <- order(TWOMINUSec, decreasing = TRUE)[1:select]
      TWOMINUSec[-indices] <- 0
    }
    igraph::plot.igraph(gp, layout = layoutp, edge.width = TWOMINUSec * thick.ec, vertex.frame.color = "darkgrey", vertex.color = "white",
                        vertex.label.color = "darkgrey", vertex.label = nomsp,  main = "Emp. Inverse Extremal Coefficients Graph")
    message <- paste("k =", k ,"    ","thick.ec = ", thick.ec)
    if(!is.null(select)){if(select != d*(d-1)/2){ message <- paste(message,"    ", "nb.selected.edges = ", select, "\n vertices.unprinted = ",deparse(noms[(1:d)[-present]],width.cutoff = 500))}}
    mtext(message, side = 1, col = "gray")
  }
}

#' Empirical graphs drawn on geographical maps of the tail dependence structure.
#'
#' Empirical tail dependograph and  Empirical inverse extremal coefficients graph drawn on geographical maps for the tail dependence structure of \code{sample} associated with threshold k.
#'
#' @name graphsMapEmp
#'
#' @param sample  A \code{(n times d)} matrix.
#'
#' @param k An integer smaller or equal to \code{n}.
#'
#' @param which A character string: \code{both}, \code{taildependograph} (the default) or \code{iecgraph}.
#'
#' @param names A character vector for \code{sample} columns which replaces \code{as.character(1:d)} (the default ones).
#'
#' @param coord Latitudes and Longitudes associated with \code{sample} columns associated to \code{region} map when \code{region} is furnished.
#'
#' @param region  A geographical region from \code{maps} package. The default value is NULL.
#'
#' @param select If select is NULL (the default) all edges are plotted. If select is an integer between 1 and the number of possible pairs of components of sample, then only the select largest edges are plotted.
#'
#' @param thick.td A numeric value for the maximal thickness of edges in \code{taildependograph}. Default value is 5.
#'
#' @param thick.ec A numeric value for the maximal thickness of edges in \code{iecgraph}. Default value is 5.
#'
#' @param eps A numerical graphical value fixing the distance between the plotted point and its names.  The default value is 0.03.
#'
#' @returns
#'
#'  It returns both (or one among) the empirical tail dependograph and the empirical inverse extremal coefficients graph on a geographical map of the  \code{sample}.
#'
#'  The empirical tail dependograph on a geographical map represents the pairwise empirical  tail superset importance coefficients of the locations associated with \code{sample} columns, see  Mercadier, C. and Roustant, O. (2019).
#'  These indices are computed by the function \code{tsicEmp}.
#'  It measures how much a pair of components (included supersets of this pair of components) is involved in the asymptotic dependence of the sample.
#'
#'  The empirical inverse extremal coefficients graph  on a geographical map represents empirical pairwise coefficients of the locations associated with \code{sample} columns that estimate \eqn{1-\ell(1_i,1_j,\bold{0})/2}.
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{graphsEmp}}
#'
#' @references
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#'
#' Becker, R. A.,  Wilks, A. R. (Original S code),  Brownrigg, R. (R version),  Minka, T. P. and  Deckmyn A. (Enhancements).  (2022)
#' maps : Draw Geographical Maps.
#' \emph{R package version 3.4.1.}
#'
#' @examples
#'
#' data(France)
#'
#' ## Figure 9 (a) of Mercadier  and Roustant (2019).
#' graphsMapEmp(France$ymt, k = 55,
#'  coord = France$coord,  region = 'France', select = 9)
#'
#' ## Figure 9 (b) of Mercadier  and Roustant (2019).
#' graphsMapEmp(France$ymt, k = 55,
#'  coord = France$coord,  region = 'France', select = 30)
#'
#' ## Figure 9 (c) of Mercadier  and Roustant (2019).
#' graphsMapEmp(France$ymt, k = 55,
#'  coord = France$coord,  region = 'France')
#'

#' @export graphsMapEmp

graphsMapEmp <- function(sample, k, which = "taildependograph", names = NULL,  coord, region = NULL, select = NULL, thick.td = 5, thick.ec = 5, eps = 0.03){
  if(!((which == "both")|(which == "taildependograph")|(which == "iecgraph"))){stop("The entry which should be a character string: 'both', 'taildependograph' or 'iecgraph'.")}

  d <- ncol(sample)

  if(!is.null(names)){noms <- names}else{noms <- 1:d}

  if((which == "both")|(which == "taildependograph")){
    res <- tsicEmp(sample, ind = 2, k = k, sobol = TRUE, norm = TRUE)
    if(!is.null(select)){
      indices <- order(res$tsic, decreasing = TRUE)[1:select]
      res$tsic[-indices] <- 0
    }
    norm_tsic <- res$tsic*thick.td
    if(is.null(select)){select = length(res$subsets)}
    if(!is.null(region)){
      maps::map(regions = region)
    }else{
      plot(coord$lon,coord$lat, pch = 20, col = 1,xaxt="n",yaxt="n",axes=FALSE,xlab="",ylab="")
    }
    for(a in 1:select){
      m <- order(res$tsic, decreasing = TRUE)[a]
      lines(coord$lon[res$subsets[[m]]], coord$lat[res$subsets[[m]]], lwd = norm_tsic[m], col = "gray")
    }
    points(coord$lon,coord$lat, pch = 20, col = 1)
    text(coord$lon,coord$lat-eps, labels = noms, cex = .8)
    title("Emp. Tail Dependograph")
    message <- paste("k =", k ,"    ","thick.td = ",thick.td)
    if(!is.null(select)){if(select != d*(d-1)/2){ message <- paste(message,"    ", "select = ", select)}}
    mtext(message, side = 1, col = "gray")
  }
  if((which == "both")|(which == "iecgraph")){
    res_ec <- ecEmp(sample, ind = 2, k = k)
    TWOMINUSec <- (2-res_ec$ec)/2
    if(!is.null(select)){
      indices <- order(TWOMINUSec, decreasing = TRUE)[1:select]
      TWOMINUSec[-indices] <- 0
    }
    norm_TWOMINUSec <- TWOMINUSec*thick.ec
    if(!is.null(region)){
      maps::map(regions = region)
    }else{
      plot(coord$lon,coord$lat, pch = 20, col = 1,xaxt="n",yaxt="n",axes=FALSE,xlab="",ylab="")
    }
    if(is.null(select)){select = length(res_ec$subsets)}
    for(a in 1:select){
      m <- order(TWOMINUSec, decreasing = TRUE)[a]
      lines(coord$lon[res_ec$subsets[[m]]], coord$lat[res_ec$subsets[[m]]], lwd = norm_TWOMINUSec[m], col = "gray")
    }
    points(coord$lon, coord$lat, pch = 20, col = 1)
    text(coord$lon, coord$lat-eps, labels = noms, cex = .8)
    title("Emp. Inverse Extremal Coeff. Graph")
    message <- paste("k =", k ,"    ","thick.ec = ", thick.ec)
    if(!is.null(select)){if(select != d*(d-1)/2){ message <- paste(message,"    ", "select = ", select)}}
    mtext(message, side = 1, col = "gray")
  }
}


#' r function for Archimax Mevlog models.
#'
#' Random vectors generation for some Archimax Mevlog models.
#'
#' @param n The number of observations.
#' @param ds An object of class \code{ds}.
#' @param dist The underlying distribution. A character string among \code{"exp"} (the default value), \code{"gamma"} and \code{"ext"}.
#'
#' @param dist.param The parameter associated with the choice \code{dist}.  If \code{dist} is \code{"exp"}, then \code{dist.param} is a postive real, the parameter of an exponential distribution. The default value is 1. If \code{dist} is \code{"gamma"}, then \code{dist.param} is a vector that concatenates the shape  and  scale parameters (in this order) of a gamma distribution.
#'
#' @returns returns a \code{n x d} matrix containing \code{n} realizations of a \code{d}-variate Archimax Mevlog random vector.
#'
#' @details
#'
#' We follow below  Algorithm 4.1 of p. 124 in Charpentier et al. (2014). Let \eqn{\psi} defined by \eqn{\psi(x)=\int_0^\infty \exp(-x t) dF(t)}, the Laplace transform of a positive random variable with  cumulative distribution function \eqn{F}.
#'
#' Define the random vector \eqn{(U_1,...,U_d)} as \eqn{U_i=\psi(-\log(Y_i)/V)} where
#' \itemize{
#' \item{} \eqn{Z} has a multivariate extreme value distribution with stable tail dependence function \eqn{\ell} ; here \eqn{Z} has standard Frechet margins,
#' \item{} \eqn{(Y_1,...,Y_d)=(\exp(-1/Z_1),...,\exp(-1/Z_d))} the margin transform of \eqn{Z} so that \eqn{Y} is sampled from the extreme value copula associated with \eqn{\ell},
#' \item{} \eqn{V} has the distribution function \eqn{F},
#' \item{} \eqn{Y} and \eqn{V} are independent.
#' }
#' Then,  \eqn{U} is sampled from the Archimax copula \eqn{C(u_1,...,u_d) = \psi(\ell(\psi^{-1}(u_1),...,\psi^{-1}(u_d)))}.
#'
#' We restrict here the function \eqn{\ell} to those associated with Mevlog models. See \code{\link[satdad]{ellMevlog}} and \code{\link[satdad]{gen.ds}}.
#'
#' We restrict also the distribution of \eqn{V} to
#' \itemize{
#' \item{} exponential ; For a positive \eqn{\lambda}, set \eqn{dF(t)=\lambda \exp(-\lambda t) 1_{t>0} dt},  then \eqn{\psi(x)=\frac{\lambda}{x+\lambda}} and \eqn{\psi^{-1}(x)=\lambda \frac{1-x}{x}}.
#'  \item{} gamma ; For positive scale \eqn{\sigma} and shape \eqn{a}, set \eqn{dF(t)= \frac{1}{\sigma^a \Gamma(a)}t^{a-1}\exp(-t/\sigma)1_{t>0}}, then \eqn{\psi(x)=\frac{1}{(x+\sigma)^a}} and \eqn{\psi^{-1}(x)=\frac{x^{-1/a}-1}{\sigma}}.
#'  }
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#' \code{\link[satdad]{rMevlog}},    \code{\link[satdad]{copArchimaxMevlog}}, \code{\link[satdad]{psiArchimaxMevlog}},  \code{\link[satdad]{psiinvArchimaxMevlog}}, \code{\link[satdad]{gen.ds}}
#'
#' @export rArchimaxMevlog
#'
#' @references
#'
#' Charpentier, A.,   Fougères,  A.-L.,  Genest, C. and  Nešlehová, J.G. (2014)
#' Multivariate Archimax copulas.
#' \emph{Journal of Multivariate Analysis}, \bold{126}, 118--136.
#'
#' @examples
#'
#' ## Fix a  5-dimensional asymmetric tail dependence structure
#' (ds5 <- gen.ds(d = 5))
#'
#' ## Generate a 1000-sample of Archimax Mevlog random vectors
#' ## associated with ds5 and underlying distribution gamma
#' (shape5 <- runif(1, 0.01, 5))
#' (scale5 <- runif(1, 0.01, 5))
#' sample5.gamma <- rArchimaxMevlog(n = 1000, ds = ds5, dist = "gamma", dist.param = c(shape5, scale5))
#'
#' ## Compare theoretical (left) and empirical (right) tail dependographs
#' oldpar <- par(mfrow = c(1,2))
#' graphs(ds = ds5)
#' graphsEmp(sample = sample5.gamma, k = 100)
#' par(oldpar)
#'
#' ## Generate a 1000-sample of Archimax Mevlog random vectors
#' ## associated with ds5 and underlying distribution exp
#' (lambda <- runif(1, 0.01, 5))
#' sample5.exp <- rArchimaxMevlog(n = 1000, ds = ds5, dist = "exp", dist.param = lambda)
#' ## Compare theoretical (left) and empirical (right) tail dependographs
#' graphs(ds = ds5)
#' graphsEmp(sample = sample5.exp, k = 100)
#'

rArchimaxMevlog <- function(n, ds,  dist = "exp", dist.param = 1){
  d <- ds$d
  z <- rMevlog(n, ds, mar = c(1,1,1))
  y <- exp(-1/z)
  if(dist=="exp"){
    if(length(dist.param) != 1){stop("The entry dist.param should be a positive real.")
      }else{
      v <- rexp(1,dist.param)
      samp <- psiArchimaxMevlog(-log(y)/v, dist = dist, dist.param = dist.param)
      }
  }
  if(dist=="gamma"){
    if(length(dist.param) != 2){stop("The entry dist.param should be a shape  and a scale, two positive parameters.")
      }else{
        v.shape <- dist.param[1]
        v.scale <- dist.param[2]
        v <- rgamma(1, shape = v.shape, scale = v.scale)
        samp <- psiArchimaxMevlog(-log(y)/v, dist = dist, dist.param = dist.param)
      }
  }
  if(dist=="ext"){
    samp <- y
  }
   return(samp)
}

#' cop-ell-psi-psiinv- functions for Archimax Mevlog models.
#'
#' Copula function, stable tail dependence function, psi function, psi inverse function for Archimax Mevlog models.
#'
#' @aliases copArchimaxMevlog
#' @aliases ellArchimaxMevlog
#' @aliases psiArchimaxMevlog
#' @aliases psiinvArchimaxMevlog
#'
#' @usage
#'
#' copArchimaxMevlog(x, ds,  dist = "exp", dist.param = 1)
#' ellArchimaxMevlog(x, ds)
#' psiArchimaxMevlog(t, dist = "exp", dist.param = 1)
#' psiinvArchimaxMevlog(t, dist = "exp", dist.param = 1)
#'
#' @param x A vector of size \code{d} or a \code{(N.x times d)} matrix.
#'
#' @param ds An object of class \code{ds}.
#'
#' @param dist The underlying distribution. A character string among \code{"exp"} (the default value), \code{"gamma"} and \code{"ext"}.
#'
#' @param dist.param The parameter associated with the choice \code{dist}.  If \code{dist} is \code{"exp"}, then \code{dist.param} is a postive real, the parameter of an exponential distribution. The default value is 1. If \code{dist} is \code{"gamma"}, then \code{dist.param} is a vector that concatenates the shape  and  scale parameters (in this order) of a gamma distribution.
#'
#' @param t A non negative scalar or vector.
#'
#'
#' @details
#'
#' The tail dependence structure is set by a \code{ds} object.  See Section \bold{Value} in \code{\link[satdad]{gen.ds}}.
#'
#' Turning to Archimax structures, we follow Charpentier et al. (2014). Their algorithm (4.1 of p. 124) has been applied in \code{\link[satdad]{rArchimaxMevlog}} to generate observations sampled from the copula
#'
#' \eqn{C(x_1,...,x_d) = \psi(\ell(\psi^{-1}(x_1),...,\psi^{-1}(x_d)))}
#'
#' when \eqn{\ell} is here the stable tail dependence function of a Mevlog model. In this package, the stdf function \eqn{\ell} is completely characterized by the \code{ds} object. See \code{\link[satdad]{ellMevlog}}.
#'
#' @returns
#'
#' When the underlying distribution \code{dist} is
#' \itemize{
#' \item{} "exp" ; For a positive \eqn{\lambda} given by \code{dist.param}, \eqn{\psi(t)=\frac{\lambda}{t+\lambda}} and \eqn{\psi^{-1}(t)=\lambda \frac{1-t}{t}}.
#'  \item{} "gamma" ; For positive scale \eqn{\sigma} and shape \eqn{a} given by \code{dist.param}, \eqn{\psi(t)=\frac{1}{(t+\sigma)^a}} and \eqn{\psi^{-1}(t)=\frac{t^{-1/a}-1}{\sigma}}.
#'  \item{} "ext" ; \eqn{\psi(t)=\exp(-t)} and \eqn{\psi^{-1}(t)=-\ln(t)}.
#' }
#'
#'  \code{copArchimaxMevlog} returns the copula function \eqn{C(x_1,...,x_d) = \psi(\ell(\psi^{-1}(x_1),...,\psi^{-1}(x_d)))}.
#'
#'  \code{ellArchimaxMevlog} returns the stable tail dependence function  \eqn{\ell(x_1,...,x_d)}.
#'
#'  \code{psiArchimaxMevlog} returns the psi function  \eqn{\psi(t)}.
#'
#'  \code{psiinvArchimaxMevlog} returns the psi inverse function  \eqn{\psi^{-1}(t)}.
#'
#'
#' @references
#'
#' Charpentier, A.,   Fougères,  A.-L.,  Genest, C. and  Nešlehová, J.G. (2014)
#' Multivariate Archimax copulas.
#' \emph{Journal of Multivariate Analysis}, \bold{126}, 118--136.
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{rArchimaxMevlog}},  \code{\link[satdad]{gen.ds}}, \code{\link[satdad]{ellMevlog}}
#'
#' @examples
#'
#' ## Fix a 7-dimensional tail dependence structure
#' ds7 <- gen.ds(d = 7)
#'
#' ## Fix the parameters for the underlying distribution
#' (lambda <- runif(1, 0.01, 5))
#' (shape <- runif(1, 0.01, 5))
#' (scale <- runif(1, 0.01, 5))
#'
#' ## Fix x and t
#' x <- c(0.8, 0.9, 0.5, 0.8, 0.4, 0.9, 0.9)
#' t <- 2
#'
#' ## Evaluate the functions under the underlying exponential construction
#' copArchimaxMevlog(x = x, ds = ds7, dist = "exp", dist.param = lambda)
#' ellArchimaxMevlog(x = x, ds = ds7)
#' psiArchimaxMevlog(t = t, dist = "exp", dist.param = lambda)
#' psiinvArchimaxMevlog(t = t, dist = "exp", dist.param = lambda)
#'
#' ## Evaluate the functions under the underlying gamma construction
#' copArchimaxMevlog(x = x, ds = ds7, dist = "gamma", dist.param = c(shape, scale))
#' ellArchimaxMevlog(x = x, ds = ds7)
#' psiArchimaxMevlog(t = t, dist = "gamma", dist.param = c(shape, scale))
#' psiinvArchimaxMevlog(t = t, dist = "gamma", dist.param = c(shape, scale))
#'
#' @export copArchimaxMevlog
#' @export ellArchimaxMevlog
#' @export psiArchimaxMevlog
#' @export psiinvArchimaxMevlog
#'
copArchimaxMevlog <- function(x, ds,  dist = "exp", dist.param = 1){
  copArchimaxMevlog_vec <- function(x.loc, ds.loc,  dist.loc, dist.param.loc){
    y.loc <- psiinvArchimaxMevlog(x, dist = dist.loc, dist.param = dist.param.loc)
    t.loc <- ellArchimaxMevlog(x = y.loc, ds = ds.loc)
    psiArchimaxMevlog(t.loc, dist = dist.loc, dist.param = dist.param.loc)
  }
  if(is.vector(x)){
    copArchimaxMevlog_vec(x.loc = x, ds.loc = ds,  dist.loc = dist, dist.param.loc = dist.param)
  }else{
    apply(x, 1, copArchimaxMevlog_vec, ds.loc = ds, dist.loc = dist, dist.param.loc = dist.param)
  }
}

ellArchimaxMevlog <- function(x, ds){
  ellMevlog(x, ds)
}

psiArchimaxMevlog <- function(t, dist = "exp", dist.param = 1){
  if(dist == "exp"){res <- dist.param/(t+dist.param)}
  if(dist == "gamma"){
    v.shape <- dist.param[1]
    v.scale <- dist.param[2]
    res <- 1/(t+v.scale)^v.shape
  }
  if(dist == "ext"){
    res <- exp(-t)
  }
  return(res)
}

psiinvArchimaxMevlog <- function(t, dist = "exp", dist.param = 1){
  if(dist == "exp"){res=dist.param*(1-t)/t}
  if(dist == "gamma"){
    v.shape <- dist.param[1]
    v.scale <- dist.param[2]
    res <- (t^{-1/v.shape}-1)/v.scale}
  if(dist == "ext"){
    res <- -log(t)
  }
  return(res)
}




#' Empirical tail importance coefficients.
#'
#' Computes on a sample the tail  importance coefficients (tic) associated with threshold \code{k}. The value may be renormalized by the empirical global variance (Sobol version).
#'
#' @param sample A \code{(n times d)} matrix.
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = 2}, all pairwise coefficients are computed.
#' @param k An integer smaller or equal to \code{n}.
#' @param sobol A boolean. `FALSE` (the default). If `TRUE`:  the index is normalized by the empirical global variance.
#'
#' @details
#' The theoretical functional decomposition of the variance of the stdf \eqn{\ell} consists in writing \eqn{D(\ell) = \sum_{I \subseteq \{1,...,d\}} D_I(\ell) } where \eqn{D_I(\ell)} measures the variance of \eqn{\ell_I(U_I)} the term associated with subset \eqn{I} in the Hoeffding-Sobol decomposition of \eqn{\ell}
#' ; note that \eqn{U_I} represents  a random vector with independent standard uniform entries. The theoretical tail variance contribution is thus \eqn{D_I(\ell)} and the theoretical tail sobol index is \eqn{S_I(\ell)=\dfrac{D_I(\ell)}{D(\ell)}}.
#'
#' Here, the function \code{ticEmp} evaluates  \eqn{\hat{D}_{I,k,n}} the empirical counterpart of \eqn{D_I(\ell)} under the option \code{sobol = FALSE}, and \eqn{\hat{S}_{I,k,n}} the empirical counterpart of \eqn{S_I(\ell)} under the option \code{sobol = TRUE}.
#'
#' Proposition 1 and Theorem 2 of Mercadier and Roustant (2019) furnish their rank-based expressions. For the subset of components \eqn{I},
#'
#' \eqn{\hat{D}_{I,k,n}=\frac{1}{k^2}\sum_{s=1}^n\sum_{s^\prime=1}^n \prod_{t\in I}(\min(\overline{R}^{(t)}_s,\overline{R}^{(t)}_{s^\prime})-\overline{R}^{(t)}_{s}\overline{R}^{(t)}_{s^\prime}) \prod_{t\notin I} \overline{R}^{(t)}_s\overline{R}^{(t)}_{s^\prime}}
#'
#' \eqn{\hat{D}_{k,n}=\frac{1}{k^2}\sum_{s=1}^n\sum_{s^\prime=1}^n \prod_{t\in I}\min(\overline{R}^{(t)}_s,\overline{R}^{(t)}_{s^\prime})- \prod_{t\in I}\overline{R}^{(t)}_{s}\overline{R}^{(t)}_{s^\prime}}
#'
#' and \eqn{\hat{S}_{I,k,n}=\dfrac{\hat{D}_{I,k,n}}{\hat{D}_{k,n}}}
#'
#' where
#'
#' \itemize{
#' \item{} \eqn{k} is the threshold parameter,
#' \item{} \eqn{n} is the sample size,
#' \item{} \eqn{X_1,...,X_n} describes the \code{sample}, each \eqn{X_s} is a d-dimensional vector \eqn{X_s^{(t)}} for \eqn{t=1,...,d},
#' \item{} \eqn{R^{(t)}_s} denotes the rank of \eqn{X^{(t)}_s} among \eqn{X^{(t)}_1, ..., X^{(t)}_n},
#' \item{} and  \eqn{\overline{R}^{(t)}_s = \min((n- R^{(t)}_s+1)/k,1)}.
#' }
#'
#' @return
#' The function returns a list of two elements:
#' \itemize{
#' \item{\code{subsets}} A list of subsets from  \eqn{\{1,...,d\}}.
#'
#' When \code{ind} is given as an integer, \code{subsets} is the list of subsets from  \eqn{\{1,...,d\}} with cardinality \code{ind}. When \code{ind} is the list, it corresponds to \code{subsets}.
#'
#' When \code{ind = "with.singletons"}  subsets is the list of all non empty subsets in \eqn{\{1,...,d\}}.
#'
#' When \code{ind = "all"}   subsets is the list of all subsets in \eqn{\{1,...,d\}} with cardinality larger or equal to 2.
#'
#' \item{\code{tic}} A vector of tail importance coefficients, or their sobol versions when \code{sobol = "TRUE"}.
#' }
#'
#' @references
#'
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#' @export ticEmp
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{tic}} and \code{\link[satdad]{tsicEmp}}
#'
#' @examples
#'
#' ## Fix a 5-dimensional asymmetric tail dependence structure
#' (ds5 <- gen.ds(d = 5))
#'
#' ## Generate a 1000-sample of Mevlog random vectors associated with ds5
#' sample5 <- rMevlog(n = 1000, ds = ds5)
#'
#' ## Compute empirical tic values according cardinality
#' res2 <- ticEmp(sample5, ind = 2, k = 100, sobol = TRUE)
#' res3 <- ticEmp(sample5, ind = 3, k = 100, sobol = TRUE)
#' res4 <- ticEmp(sample5, ind = 4, k = 100, sobol = TRUE)
#'
#' ## Represent the empirical indices associated with pairs
#' barplot(res2$tic ~ as.character(res2$subsets), las = 2,
#'      xlab = "", ylab = "", main = "Tail Sobol Indices (cardinality 2)")
#'
#' ## Represent the empirical indices associated with triplets
#' barplot(res3$tic ~ as.character(res3$subsets), las = 2,
#'      xlab = "", ylab = "", main = "Tail Sobol Indices (cardinality 3)")
#'
#' ## Represent the empirical indices associated with quadriplets
#' barplot(res4$tic ~ as.character(res4$subsets), las = 2,
#'      xlab = "", ylab ="", main = "Tail Sobol Indices (cardinality 4)")
#'
#' ## Check the sum-to-one constraint of empirical tail Sobol indices
#' sum(ticEmp(sample5, ind = "with.singletons", k = 100,  sobol = TRUE)$tic)
#'


ticEmp <- function(sample, ind = 2, k, sobol = FALSE){
  if(!is.list(ind)){
    if(ind == "all"){
      ind <-  .subsets_cpp(ncol(sample))
      ind <- ind[-(1:ncol(sample))]
    }else if(ind == "with.singletons"){
      ind <-  .subsets_cpp(ncol(sample))
    }else{
      ind <- as.list(as.data.frame(combn(1:ncol(sample), ind)))
      names(ind) <- NULL
    }
  }
  Rbar <- .rankbar_cpp(sample, k)
  ratio <- .ratioallsobol_cpp(Rbar)
  prodprod <- .prodprod_cpp(Rbar)
  prodmin <-  .prodmin_cpp(Rbar)
  N = length(ind)
  di <- vector("list", N)
  for(m in 1:N){
    sub <- ind[[m]]
    aux <- 1
    for(i in sub){aux <- aux*ratio[[i]]}
    di[[m]] <- aux*prodprod
  }
  di <- (as.numeric(lapply(di, sum))/k^2)
  if(sobol){
    D <- sum(prodmin-prodprod)/k^2
    di <- di/D}
  return(list(subsets = ind, tic = di))
}


#' Tail importance coefficients for Mevlog models.
#'
#' Computes  the tail  importance coefficients (tic) on a \code{Mevlog} model which is a multivariate extreme value (symmetric or asymmetric) logistic model, descibed here by its dependence structure.
#'
#' @param ds An object of class \code{ds}.
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = 2}, all pairwise coefficients are computed.
#' @param n.MC Monte Carlo sample size. Default value is 1000. See details in  \code{\link[satdad]{tsic}}.
#' @param sobol A boolean. `FALSE` (the default). If `TRUE`:  the index is normalized by the theoretical global variance.
#'
#' @details
#' The tail dependence structure is specified using a \code{ds} object, which corresponds to the stable tail dependence function  \eqn{\ell}.
#' The process for deducing  the stable tail dependence function \eqn{\ell} from \code{ds} is explained in the Details section of \code{\link[satdad]{gen.ds}}.
#'
#' The theoretical functional decomposition of the variance of the stdf \eqn{\ell} consists in writing \eqn{D(\ell) = \sum_{I \subseteq \{1,...,d\}} D_I(\ell) } where \eqn{D_I(\ell)} measures the variance of \eqn{\ell_I(U_I)} the term associated with subset \eqn{I} in the Hoeffding-Sobol decomposition of \eqn{\ell}
#' ; note that \eqn{U_I} represents  a random vector with independent standard uniform entries.
#' The theoretical tail importance coefficient (tic) is thus \eqn{D_I(\ell)} and  its sobol version is \eqn{S_I(\ell)=\dfrac{D_I(\ell)}{D(\ell)}}.
#'
#' The function \code{tic} uses the Mobius inversion formula, see Formula (8) in Liu and Owen (2006), to derive the tic from the tsic. The latter are the tail superset importance coefficients obtained by the function \code{\link[satdad]{tsic}}.
#'
#' @return
#' The function returns a list of two elements:
#' \itemize{
#' \item{\code{subsets}} A list of subsets from  \eqn{\{1,...,d\}}.
#'
#' When \code{ind} is given as an integer, \code{subsets} is the list of subsets from  \eqn{\{1,...,d\}} with cardinality \code{ind}.
#'
#' When \code{ind} is a list, it corresponds to \code{subsets}.
#'
#' When \code{ind = "with.singletons"}  subsets is the list of all non empty subsets in \eqn{\{1,...,d\}}.
#'
#' When \code{ind = "all"}   subsets is the list of all subsets in \eqn{\{1,...,d\}} with cardinality larger or equal to 2.
#'
#' \item{\code{tic}} A vector of tail importance coefficients, or their Sobol versions when \code{sobol = "TRUE"}.
#' }
#'
#' @seealso
#'
#'  \code{\link[satdad]{tsic}}, \code{\link[satdad]{ticEmp}}
#'
#' @references
#' Liu, R. and Owen, A. B. (2006)
#' Estimating mean dimensionality of analysis of variance decompositions.
#' \emph{J. Amer. Statist. Assoc.}, \bold{101(474)}:712--721.
#'
#' Mercadier, C. and Roustant, O. (2019)
#' The tail dependograph.
#' \emph{Extremes}, \bold{22}, 343--372.
#'
#' @export tic
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{ticEmp}} and \code{\link[satdad]{tsic}}
#'
#' @examples
#'
#' ## Fix a 4-dimensional asymmetric tail dependence structure
#' ds4 <- gen.ds(d = 4, sub = list(1:2,3:4,1:3))
#'
#' ## Compute all tic values
#' res4 <- tic(ds4, ind = "with.singletons", sobol = TRUE)
#'
#' ## Check the sum-to-one constraint of tail Sobol indices
#' sum(res4$tic)
#'


tic <- function(ds, ind = 2, n.MC = 1000, sobol = FALSE){
  res <- tsic(ds = ds, ind = "with.singletons", n.MC = n.MC, norm = FALSE)
  sig <- (-1)^unlist(lapply(res$subsets, length))
  tsic <- c(-sum(res$tsic*sig), res$tsic)

  mobiusInversion <- function(coefficients){
    d.loc <- as.integer(log2(length(coefficients)))
    subsets.loc <- vector("list",2^d.loc)
    subsets.loc[-1] <-  .subsets_cpp(d.loc)
    transformed_coefficients <- rep(NA,2^d.loc)
    for(i in 1:(2^d.loc)){
      sub.loc <-  subsets.loc[[i]]
      ind.loc <- lapply(subsets.loc,function(v){prod(sub.loc %in% v)})
      signe.loc <- (-1)^(unlist(lapply(subsets.loc[which(ind.loc==1)],length))-length(sub.loc))
      transformed_coefficients[i] <- sum(coefficients[which(ind.loc==1)]*signe.loc)
    }
    transformed_coefficients[transformed_coefficients<0]=0
    transformed_coefficients
  }

  tic <- mobiusInversion(tsic)
  if(sobol){tic <- tic/tsic[1]}

  if(ind == "all"){
    tic <- tic[-(1:(ds$d+1))]; subsets <-  res$subsets[-(1:ds$d)]}
  if(ind == "with.singletons"){tic <- tic[-1]; subsets <-  res$subsets}
  if((!is.list(ind))&(is.numeric(ind))){
    k <- sum(choose(ds$d,0:(ind-1)))
    nk <- choose(ds$d,ind)
    select <- k+0:(nk-1)
    tic <- tic[-1]
    tic <- tic[select]
    subsets <- res$subsets[select]
  }
  return(list(tic = tic, subsets =subsets))
}



#' Cleveland's Dot Plots of the tail dependence structure.
#'
#' Global comparison of the theoretical tail superset importance coefficients (tsic) via a Cleveland's Dot Plot.
#'
#' @param ds An object of class \code{ds}.
#'
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = "all"}.
#'
#' @param which A character string among "tsic" (normalized tsic plot), and "ec" (normalized ec plot).
#'
#' @param labels A boolean. `TRUE` the default indicates that the names of the subsets are printed. `FALSE` if only points are drawn.
#'
#' @return
#'
#' Draws a Cleveland dot plot of the normalized theoretical tsic when \code{superset = TRUE}, the default value.
#'
#' Otherwise theoretical normalized ec are drawn.
#'
#' @export plotClev
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{plotClevEmp}}
#'
#' @examples
#' ## Fix a 6-dimensional asymmetric tail dependence structure
#' ## Two blocks of components are specified
#' ds6 <- gen.ds(d = 6, sub = list(1:4,5:6))
#'
#' ## Plot the associated Cleveland dot plot
#' plotClev(ds6)
#'
plotClev <- function(ds, ind = "all", which = "tsic", labels = TRUE){
  if((which !="tsic")&(which !="ec")){stop("The entry `which` is incorrect.")}
  if(!is.list(ind)){
    if(ind == "all"){
      ind <-  .subsets_cpp(ds$d)
      ind <- ind[-(1:ds$d)]
    }else if(ind == "with.singletons"){
      ind <-  .subsets_cpp(ds$d)
    }else{
      ind <- as.list(as.data.frame(combn(1:ds$d, ind)))
      names(ind) <- NULL
    }
  }
  if(which == "tsic"){
    res <- tsic(ds, ind = ind, sobol = TRUE, norm = TRUE)
    if(labels){
      labels = as.character(res$subsets)
    }else{
      labels = NULL
    }
    dotchart(res$tsic, groups = unlist(lapply(res$subsets,length)), labels = labels)
    title("Cleveland's Dot Plot of tsic")
  }else{
    res <- ec(ds, ind = ind, norm = TRUE)
    if(labels){
      labels = as.character(res$subsets)
    }else{
      labels = NULL
    }
    dotchart(res$ec, groups = unlist(lapply(res$subsets,length)), labels = labels)
    title("Cleveland's Dot Plot of\n  inv. normalized ec")
  }
}



#' Empirical Cleveland's Dot Plots of the tail dependence structure.
#'
#' Global comparison of the empirical tail superset importance coefficients (tsicEmp) via a Cleveland's Dot Plot.
#'
#' @param sample A \code{(n times d)} matrix.
#' @param k An integer smaller or equal to \code{n}.
#' @param ind A character string among "with.singletons" and "all" (without singletons), or an integer in \eqn{\{2,...,d\}} or a list of subsets from  \eqn{\{1,...,d\}}. The default is \code{ind = "all"}.
#' @param which A character string among "tsic" (empirical normalized tsic plot), and "ec" (empirical normalized ec plot).
#' @param labels A boolean. `TRUE` the default indicates that the names of the subsets are printed. `FALSE` if only points are drawn.
#'
#' @return
#'
#' Draws a Cleveland dot plot of the normalized empirical tsic when \code{superset = TRUE}, the default value.
#'
#' Otherwise empirical normalized ec are drawn.
#'
#' @author
#' Cécile Mercadier (\code{mercadier@math.univ-lyon1.fr})
#'
#' @seealso
#'
#' \code{\link[satdad]{plotClev}}
#'
#'
#' @export plotClevEmp
#'
#' @examples
#'
#' ## Fix a 5-dimensional asymmetric tail dependence structure
#' (ds5 <- gen.ds(d = 5))
#'
#' ## Generate a 1000-sample of Mevlog random vectors associated with ds5
#' sample5 <- rMevlog(n = 1000, ds = ds5)
#'
#' ## Plot the empirical Cleveland dot plot (restricted to pairs)
#' plotClevEmp(sample5,  k = 100,  ind = 2)
#'


plotClevEmp <- function(sample, k, ind = "all", which = "tsic", labels = TRUE){
  if((which !="tsic")&(which !="ec")){stop("The entry `which` is incorrect.")}
  if(!is.list(ind)){
    if(ind == "all"){
      ind <-  .subsets_cpp(ncol(sample))
      ind <- ind[-(1:ncol(sample))]
    }else if(ind == "with.singletons"){
      ind <-  .subsets_cpp(ncol(sample))
    }else{
      ind <- as.list(as.data.frame(combn(1:ncol(sample), ind)))
      names(ind) <- NULL
    }
  }

  if(which == "tsic"){
    res <- tsicEmp(sample, k = k, ind = ind, sobol = TRUE, norm = TRUE)
    if(labels){
      labels = as.character(res$subsets)
    }else{
      labels = NULL
    }
    dotchart(res$tsic, groups = unlist(lapply(res$subsets,length)), labels = labels)
    title("Cleveland's Dot Plot of Emp. tsic")
  }else{
    res <- ecEmp(sample, k = k, ind = ind, norm = TRUE)
    if(labels){
      labels = as.character(res$subsets)
    }else{
      labels = NULL
    }
    dotchart(res$ec, groups = unlist(lapply(res$subsets,length)), labels = labels)
    title("Cleveland's Dot Plot of\n Emp. inv. normalized ec")
  }
}
