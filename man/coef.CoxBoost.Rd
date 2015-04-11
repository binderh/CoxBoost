\name{coef.CoxBoost}
\alias{coef.CoxBoost}
\title{Coeffients from CoxBoost fit}
\description{
Extracts the coefficients from the specified boosting steps of a CoxBoost object fitted by \code{\link{CoxBoost}}.  
}
\usage{
\method{coef}{CoxBoost}(object,at.step=NULL,scaled=TRUE,\ldots)
}
\arguments{
\item{object}{fitted CoxBoost object from a \code{\link{CoxBoost}} call.}
\item{at.step}{scalar or vector of boosting step(s) at which prediction is wanted. If no step is given, the final boosting step is used.}
\item{scaled}{logical value indicating whether coefficients should be returned scaled to be at the level of the original covariate values, or unscaled as used internally when \code{standardize=TRUE} is used in the \code{\link{CoxBoost}} call.}
\item{\dots}{miscellaneous arguments, none of which is used at the moment.}
}
\value{
For a vector of length \code{p} (number of covariates) (\code{at.step} being a scalar) or a \code{length(at.step) * p} matrix (\code{at.step} being a vector).  
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
