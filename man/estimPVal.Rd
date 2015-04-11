\name{estimPVal}
\alias{estimPVal}
\title{Estimate p-values for a model fitted by CoxBoost}
\description{
Performs permutation-based p-value estimation for the optional covariates in a fit from \code{\link{CoxBoost}}.}
\usage{
estimPVal(object,x,permute.n=10,per.covariate=FALSE,parallel=FALSE,
          multicore=FALSE,trace=FALSE,...)
}
\arguments{
\item{object}{fit object obtained from \code{\link{CoxBoost}}.}
\item{x}{\code{n * p} matrix of covariates. This has to be the same that was used in the call to \code{\link{CoxBoost}}.}
\item{permute.n}{number of permutations employed for obtaining a null distribution.}
\item{per.covariate}{logical value indicating whether a separate null distribution should be considered for each covariate. A larger number of permutations will be needed if this is wanted.}
\item{parallel}{logical value indicating whether computations for obtaining a null distribution via permutation should be performed in parallel on a compute cluster. Parallelization is performed via the package \code{snowfall} and the initialization function of of this package, \code{sfInit}, should be called before calling \code{estimPVal}.}
\item{multicore}{indicates whether computations in the permuted data sets should be performed in parallel, using package \code{parallel}. If \code{TRUE}, package \code{parallel} is employed using the default number of cores. A value larger than \code{1} is taken to be the number of cores that should be employed.}
\item{trace}{logical value indicating whether progress in estimation should be indicated by printing the number of the permutation that is currently being evaluated.}
\item{\dots}{miscellaneous parameters for the calls to \code{\link{CoxBoost}}}
}
\details{
As p-value estimates are based on permutations, random numbers are drawn for determining permutation indices. Therfore, the results depend on the state of the random number generator. This can be used to explore the variability due to random variation and help to determine an adequate value for \code{permute.n}. A value of 100 should be sufficient, but this can be quite slow. If there is a considerable number of covariates, e.g., larger than 100, a much smaller number of permutations, e.g., 10, might already work well. The estimates might also be negatively affected, if only a small number of boosting steps (say <50) was employed for the original fit. 
}
\value{
Vector with p-value estimates, one value for each optional covariate specificed in the original call to \code{\link{CoxBoost}}.
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
\references{
Binder, H., Porzelius, C. and Schumacher, M. (2009). Rank-based p-values for sparse high-dimensional risk prediction models fitted by componentwise boosting. FDM-Preprint Nr. 101, University of Freiburg, Germany.
}
\seealso{
\code{\link{CoxBoost}}
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

#   Fit a Cox proportional hazards model by CoxBoost

cbfit <- CoxBoost(time=obs.time,status=status,x=x,stepno=100,
                  penalty=100) 

#   estimate p-values

p1 <- estimPVal(cbfit,x,permute.n=10)

#   get a second vector of estimates for checking how large
#   random variation is

p2 <- estimPVal(cbfit,x,permute.n=10)

plot(p1,p2,xlim=c(0,1),ylim=c(0,1),xlab="permute 1",ylab="permute 2")
}

}
\keyword{models} \keyword{regression} \keyword{survial}
