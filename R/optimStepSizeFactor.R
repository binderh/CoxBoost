#' Coarse line search for optimum step-size modification factor
#'
#' This routine helps in finding an optimum step-size modification factor for
#' \code{\link{CoxBoost}}, i.e., that results in an optimum in terms of
#' cross-validated partial log-likelihood.
#'
#' A coarse line search is performed for finding the best parameter
#' \code{stepsize.factor} for \code{\link{CoxBoost}}. If an \code{pendistmat}
#' argument is provided (which is passed on to \code{\link{CoxBoost}}), a
#' search for factors smaller than 1 is sensible (corresponding to
#' \code{direction="down"}). If no connection information is provided, it is
#' reasonable to employ \code{direction="both"}, for avoiding restrictions
#' without subject matter knowledge.
#'
#' @param time vector of length \code{n} specifying the observed times.
#' @param status censoring indicator, i.e., vector of length \code{n} with
#' entries \code{0} for censored observations and \code{1} for uncensored
#' observations. If this vector contains elements not equal to \code{0} or
#' \code{1}, these are taken to indicate events from a competing risk and a
#' model for the subdistribution hazard with respect to event \code{1} is
#' fitted (see e.g. Fine and Gray, 1999).
#' @param x \code{n * p} matrix of covariates.
#' @param direction direction of line search for an optimal step-size
#' modification factor (starting from value 1).
#' @param start.stepsize step size used for the line search. A final step is
#' performed using half this size.
#' @param iter.max maximum number of search iterations.
#' @param constant.cv.res result of \code{\link{cv.CoxBoost}} for
#' \code{stepsize.factor=1}, that can be provided for saving computing time, if
#' it already is available.
#' @param parallel logical value indicating whether computations in the
#' cross-validation folds should be performed in parallel on a compute cluster.
#' Parallelization is performed via the package \code{snowfall} and the
#' initialization function of of this package, \code{sfInit}, should be called
#' before calling \code{cv.CoxBoost}.
#' @param trace logical value indicating whether information on progress should
#' be printed.
#' @param \dots miscellaneous parameters for \code{\link{cv.CoxBoost}}.
#' @return List with the following components: \item{factor.list}{array with
#' the evaluated step-size modification factors.} \item{critmat}{matrix with
#' the mean partial log-likelihood for each step-size modification factor in
#' the course of the boosting steps.} \item{optimal.factor.index}{index of the
#' optimal step-size modification factor.} \item{optimal.factor}{optimal
#' step-size modification factor.} \item{optimal.step}{optimal boosting step
#' number, i.e., with minimum mean partial log-likelihood, for step-size
#' modification factor \code{optimal.factor}.}
#' @author Written by Harald Binder \email{binderh@@uni-mainz.de}.
#' @seealso \code{\link{CoxBoost}}, \code{\link{cv.CoxBoost}}
#' @references Binder, H. and Schumacher, M. (2009). Incorporating pathway
#' information into boosting estimation of high-dimensional risk prediction
#' models. BMC Bioinformatics. 10:18.
#' @keywords models smooth regression
#' @examples
#'
#' \dontrun{
#' #   Generate some survival data with 10 informative covariates
#' n <- 200; p <- 100
#' beta <- c(rep(1,10),rep(0,p-10))
#' x <- matrix(rnorm(n*p),n,p)
#' real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
#' cens.time <- rexp(n,rate=1/10)
#' status <- ifelse(real.time <= cens.time,1,0)
#' obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)
#'
#' #  Determine step-size modification factor. As there is no connection matrix,
#' #  perform search into both directions
#'
#' optim.res <- optimStepSizeFactor(direction="both",
#'                                 time=obs.time,status=status,x=x,
#'                                 trace=TRUE)
#'
#' #   Fit with obtained step-size modification parameter and optimal number of boosting
#' #   steps obtained by cross-validation
#'
#' cbfit <- CoxBoost(time=obs.time,status=status,x=x,
#'                   stepno=optim.res$optimal.step,
#'                   stepsize.factor=optim.res$optimal.factor)
#' summary(cbfit)
#'
#' }
#'
#' @export
optimStepSizeFactor <- function(time,status,x,direction=c("down","up","both"),start.stepsize=0.1,iter.max=10,
                                constant.cv.res=NULL,parallel=FALSE,trace=FALSE,...)
{
    if (parallel) {
        if (!requireNamespace("snowfall")) {
            parallel <- FALSE
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            snowfall::sfExport("x")
        }
    }

    direction <- match.arg(direction)

    factor.list <- switch(direction,down=c(1,1-start.stepsize,1-2*start.stepsize),
                                    up=c(1,1+start.stepsize,1+2*start.stepsize),
                                    both=c(1-start.stepsize,1,1+start.stepsize))

    critmat <- NULL

    folds <- NULL
    if (!is.null(constant.cv.res)) folds <- constant.cv.res$folds

    i <- 1
    step.size <- start.stepsize
    reduction.done <- FALSE
    iter.count <- 0

    while (i <= length(factor.list) && iter.count < iter.max) {
        iter.count <- iter.count + 1
        if (trace) cat("iteration ",iter.count,": evaluating factor ",factor.list[i],"\n",sep="")

        if (factor.list[i] == 1 && !is.null(constant.cv.res)) {
            cv.res.act.path <- constant.cv.res
        } else {
            if (is.null(folds)) {
                cv.res.act.path <- cv.CoxBoost(time=time,status=status,x=x,stepsize.factor=factor.list[i],parallel=parallel,upload.x=FALSE,trace=trace,...)
            } else {
                cv.res.act.path <- cv.CoxBoost(time=time,status=status,x=x,stepsize.factor=factor.list[i],parallel=parallel,upload.x=FALSE,folds=folds,trace=trace,...)
            }
        }

        critmat <- rbind(critmat,cv.res.act.path$mean.logplik)

        if (is.null(folds)) folds <- cv.res.act.path$folds

        i <- i + 1
        if (i > length(factor.list)) {
            if (reduction.done) break

            actual.max.val <- apply(critmat,1,max)
            actual.max <- which.max(actual.max.val)

            do.reduction <- TRUE

            if (direction %in% c("down","both") && factor.list[actual.max] == min(factor.list)) {
                if (min(factor.list) < step.size*2) break
                factor.list <- c(factor.list,min(factor.list)-step.size)
                do.reduction <- FALSE
            }

            if (direction %in% c("up","both") && factor.list[actual.max] == max(factor.list)) {
                factor.list <- c(factor.list,max(factor.list)+step.size)
                do.reduction <- FALSE
            }

            if (do.reduction) {
                sort.factor.list <- sort(factor.list)
                max.pos <- which(sort.factor.list == factor.list[actual.max])
                max.factor <- sort.factor.list[max.pos]

                if (max.pos == 1) {
                    second.max.pos <- 2
                } else {
                    if (max.pos == length(sort.factor.list)) {
                        second.max.pos <- length(sort.factor.list) - 1
                    } else {
                        candidate1.pos <- max.pos - 1
                        candidate2.pos <- max.pos + 1

                        if (actual.max.val[which(sort.factor.list[candidate1.pos] == factor.list)] > actual.max.val[which(sort.factor.list[candidate2.pos] == factor.list)]) {
                            second.max.pos <- candidate1.pos
                        } else {
                            second.max.pos <- candidate2.pos
                        }
                    }
                }

                second.max.factor <- sort.factor.list[second.max.pos]

                factor.list <- c(factor.list,mean(c(max.factor,second.max.factor)))

                reduction.done <- TRUE
            }
        }
    }

    optimal.factor.index <- which.max(apply(critmat,1,max))

    list(factor.list=factor.list,
         critmat=critmat,
         optimal.factor.index=optimal.factor.index,
         optimal.factor=factor.list[optimal.factor.index],
         optimal.step=which.max(critmat[optimal.factor.index,])-1)
}
