\name{predict.iCoxBoost}
\alias{predict.iCoxBoost}
\title{Predict method for iCoxBoost fits}
\description{
Obtains predictions at specified boosting steps from a iCoxBoost object fitted by \code{\link{iCoxBoost}}.  
}
\usage{
\method{predict}{iCoxBoost}(object,newdata=NULL,
        subset=NULL,at.step=NULL,times=NULL,
        type=c("lp","logplik","risk","CIF"),\ldots)
}
\arguments{
\item{object}{fitted CoxBoost object from a \code{\link{CoxBoost}} call.}
\item{newdata}{data frame with new covariate values (for \code{n.new} observations). If just prediction for the training data is wanted, it can be omitted. If the predictive partial log-likelihood is wanted (\code{type=logplik}), this frame also has to contain the response information.}
\item{subset}{an optional vector specifying a subset of observations to be used for evaluation.}
\item{at.step}{scalar or vector of boosting step(s) at which prediction is wanted. If \code{type="risk"} is used, only one step is admissible. If no step is given, the final boosting step is used.}
\item{times}{vector with \code{T} time points where prediction is wanted. Only needed for \code{type="risk"}}
\item{type}{type of prediction to be returned: \code{"lp"} gives the linear predictor, \code{"logplik"} the partial log-likelihood, \code{"risk"} the predicted probability of not yet having had the event at the time points given in \code{times}, and \code{"CIF"} the predicted cumulative incidence function, i.e., the predicted probability of having had the event of interest.}
\item{\dots}{miscellaneous arguments, none of which is used at the moment.}
}
\value{
For \code{type="lp"} and \code{type="logplik"} a vector of length \code{n.new} (\code{at.step} being a scalar) or a \code{n.new * length(at.step)} matrix (\code{at.step} being a vector) with predictions is returned.  
For \code{type="risk"} or \code{type="CIF"} a \code{n.new * T} matrix with predicted probabilities at the specific time points is returned.  
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
\examples{
n <- 200; p <- 100
beta <- c(rep(1,2),rep(0,p-2))
x <- matrix(rnorm(n*p),n,p)
actual.data <- as.data.frame(x)
real.time <- -(log(runif(n)))/(10*exp(drop(x \%*\% beta)))
cens.time <- rexp(n,rate=1/10)
actual.data$status <- ifelse(real.time <= cens.time,1,0)
actual.data$time <- ifelse(real.time <= cens.time,real.time,cens.time)

#   define training and test set

train.index <- 1:100
test.index <- 101:200

#   Fit a Cox proportional hazards model by iCoxBoost

\donttest{cbfit <- iCoxBoost(Surv(time,status) ~ .,data=actual.data[train.index,],
				   stepno=300,cv=FALSE)}

#   mean partial log-likelihood for test set in every boosting step

\donttest{step.logplik <- predict(cbfit,newdata=actual.data[test.index,],
                        at.step=0:300,type="logplik")

plot(step.logplik)}

#   names of covariates with non-zero coefficients at boosting step
#   with maximal test set partial log-likelihood

\donttest{print(coef(cbfit,at.step=which.max(step.logplik)-1))}

}
\keyword{models} \keyword{regression} \keyword{survial}
