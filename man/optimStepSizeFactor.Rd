\name{optimStepSizeFactor}
\alias{optimStepSizeFactor}
\title{Coarse line search for optimum step-size modification factor}
\description{
This routine helps in finding an optimum step-size modification factor for \code{\link{CoxBoost}}, i.e., that results in an optimum in terms of cross-validated partial log-likelihood.  
}
\usage{
optimStepSizeFactor(time,status,x,
                    direction=c("down","up","both"),start.stepsize=0.1,
                    iter.max=10,constant.cv.res=NULL,
                    parallel=FALSE,trace=FALSE,...) 
}
\arguments{
\item{time}{vector of length \code{n} specifying the observed times.}
\item{status}{censoring indicator, i.e., vector of length \code{n} with entries \code{0} for censored observations and \code{1} for uncensored observations. If this vector contains elements not equal to \code{0} or \code{1}, these are taken to indicate events from a competing risk and a model for the subdistribution hazard with respect to event \code{1} is fitted (see e.g. Fine and Gray, 1999).}
\item{x}{\code{n * p} matrix of covariates.}
\item{direction}{direction of line search for an optimal step-size modification factor (starting from value 1).}
\item{start.stepsize}{step size used for the line search. A final step is performed using half this size.}
\item{iter.max}{maximum number of search iterations.}
\item{constant.cv.res}{result of \code{\link{cv.CoxBoost}} for \code{stepsize.factor=1}, that can be provided for saving computing time, if it already is available.}
\item{parallel}{logical value indicating whether computations in the cross-validation folds should be performed in parallel on a compute cluster. Parallelization is performed via the package \code{snowfall} and the initialization function of of this package, \code{sfInit}, should be called before calling \code{cv.CoxBoost}.}
\item{trace}{logical value indicating whether information on progress should be printed.}
\item{\dots}{miscellaneous parameters for \code{\link{cv.CoxBoost}}.}
}
\details{
A coarse line search is performed for finding the best parameter \code{stepsize.factor} for \code{\link{CoxBoost}}. If an \code{pendistmat} argument is provided (which is passed on to \code{\link{CoxBoost}}), a search for factors smaller than 1 is sensible (corresponding to \code{direction="down"}). If no connection information is provided, it is reasonable to employ \code{direction="both"}, for avoiding restrictions without subject matter knowledge.
}
\value{
List with the following components:
\item{factor.list}{array with the evaluated step-size modification factors.}
\item{critmat}{matrix with the mean partial log-likelihood for each step-size modification factor in the course of the boosting steps.}
\item{optimal.factor.index}{index of the optimal step-size modification factor.}
\item{optimal.factor}{optimal step-size modification factor.}
\item{optimal.step}{optimal boosting step number, i.e., with minimum mean partial log-likelihood, for step-size modification factor \code{optimal.factor}.}
}
\author{
Written by Harald Binder \email{binderh@uni-mainz.de}. 
}
\references{
Binder, H. and Schumacher, M. (2009). Incorporating pathway information into boosting estimation of high-dimensional risk prediction models. BMC Bioinformatics. 10:18.
}
\seealso{
\code{\link{CoxBoost}}, \code{\link{cv.CoxBoost}}
}
\examples{
\dontrun{
#   Generate some survival data with 10 informative covariates 
n <- 200; p <- 100
beta <- c(rep(1,10),rep(0,p-10))
x <- matrix(rnorm(n*p),n,p)
real.time <- -(log(runif(n)))/(10*exp(drop(x \%*\% beta)))
cens.time <- rexp(n,rate=1/10)
status <- ifelse(real.time <= cens.time,1,0)
obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)

#  Determine step-size modification factor. As there is no connection matrix,
#  perform search into both directions 

optim.res <- optimStepSizeFactor(direction="both",
                                time=obs.time,status=status,x=x,
                                trace=TRUE)

#   Fit with obtained step-size modification parameter and optimal number of boosting
#   steps obtained by cross-validation

cbfit <- CoxBoost(time=obs.time,status=status,x=x,
                  stepno=optim.res$optimal.step,
                  stepsize.factor=optim.res$optimal.factor) 
summary(cbfit)

}

}
\keyword{models} \keyword{smooth} \keyword{regression}
