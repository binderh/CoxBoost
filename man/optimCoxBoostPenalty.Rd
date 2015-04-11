\name{optimCoxBoostPenalty}
\alias{optimCoxBoostPenalty}
\title{Coarse line search for adequate penalty parameter}
\description{
This routine helps in finding a penalty value that leads to an ``optimal'' number of boosting steps for CoxBoost, determined by cross-validation, that is not too small/in a specified range.  
}
\usage{
optimCoxBoostPenalty(time,status,x,minstepno=50,maxstepno=200,
                     start.penalty=9*sum(status==1),iter.max=10,
                     upper.margin=0.05,parallel=FALSE,
                     trace=FALSE,...)
}
\arguments{
\item{time}{vector of length \code{n} specifying the observed times.}
\item{status}{censoring indicator, i.e., vector of length \code{n} with entries \code{0} for censored observations and \code{1} for uncensored observations. If this vector contains elements not equal to \code{0} or \code{1}, these are taken to indicate events from a competing risk and a model for the subdistribution hazard with respect to event \code{1} is fitted (see e.g. Fine and Gray, 1999).}
\item{x}{\code{n * p} matrix of covariates.}
\item{minstepno, maxstepno}{range of boosting steps in which the ``optimal'' number of boosting steps is wanted to be.}
\item{start.penalty}{start value for the search for the appropriate penalty.}
\item{iter.max}{maximum number of search iterations.}
\item{upper.margin}{specifies the fraction of \code{maxstepno} which is used as an upper margin in which a cross-validation minimum is not taken to be one. This is necessary because of random fluctuations of cross-validated partial log-likelihood.}
\item{parallel}{logical value indicating whether computations in the cross-validation folds should be performed in parallel on a compute cluster. Parallelization is performed via the package \code{snowfall} and the initialization function of of this package, \code{sfInit}, should be called before calling \code{cv.CoxBoost}.}
\item{trace}{logical value indicating whether information on progress should be printed.}
\item{\dots}{miscellaneous parameters for \code{\link{cv.CoxBoost}}.}
}
\details{
The penalty parameter for \code{\link{CoxBoost}} has to be chosen only very coarsely.  In Tutz and Binder (2006) it is suggested for likelihood based boosting just to make sure, that the optimal number of boosting steps, according to some criterion such as cross-validation, is larger or equal to 50.  With a smaller number of steps, boosting may become too ``greedy'' and show sub-optimal performance.  This procedure uses a very coarse line search and so one should specify a rather large range of boosting steps.  
}
\value{
List with element \code{penalty} containing the ``optimal'' penalty and \code{cv.res} containing the corresponding result of \code{cv.CoxBoost}.
}
\author{
Written by Harald Binder \email{binderh@uni-mainz.de}. 
}
\references{
Tutz, G. and Binder, H. (2006) Generalized additive modelling with implicit variable selection by likelihood based boosting. \emph{Biometrics}, 62:961-971.
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

#  determine penalty parameter

optim.res <- optimCoxBoostPenalty(time=obs.time,status=status,x=x,
                                  trace=TRUE,start.penalty=500)

#   Fit with obtained penalty parameter and optimal number of boosting
#   steps obtained by cross-validation

cbfit <- CoxBoost(time=obs.time,status=status,x=x,
                  stepno=optim.res$cv.res$optimal.step,
                  penalty=optim.res$penalty) 
summary(cbfit)

}

}
\keyword{models} \keyword{smooth} \keyword{regression}
