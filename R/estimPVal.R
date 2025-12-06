#' Estimate p-values for a model fitted by CoxBoost
#'
#' Performs permutation-based p-value estimation for the optional covariates in
#' a fit from \code{\link{CoxBoost}}.
#'
#' As p-value estimates are based on permutations, random numbers are drawn for
#' determining permutation indices. Therfore, the results depend on the state
#' of the random number generator. This can be used to explore the variability
#' due to random variation and help to determine an adequate value for
#' \code{permute.n}. A value of 100 should be sufficient, but this can be quite
#' slow. If there is a considerable number of covariates, e.g., larger than
#' 100, a much smaller number of permutations, e.g., 10, might already work
#' well. The estimates might also be negatively affected, if only a small
#' number of boosting steps (say <50) was employed for the original fit.
#'
#' @param object fit object obtained from \code{\link{CoxBoost}}.
#' @param x \code{n * p} matrix of covariates. This has to be the same that was
#' used in the call to \code{\link{CoxBoost}}.
#' @param permute.n number of permutations employed for obtaining a null
#' distribution.
#' @param per.covariate logical value indicating whether a separate null
#' distribution should be considered for each covariate. A larger number of
#' permutations will be needed if this is wanted.
#' @param parallel logical value indicating whether computations for obtaining
#' a null distribution via permutation should be performed in parallel on a
#' compute cluster. Parallelization is performed via the package
#' \code{snowfall} and the initialization function of of this package,
#' \code{sfInit}, should be called before calling \code{estimPVal}.
#' @param multicore indicates whether computations in the permuted data sets
#' should be performed in parallel, using package \code{parallel}. If
#' \code{TRUE}, package \code{parallel} is employed using the default number of
#' cores. A value larger than \code{1} is taken to be the number of cores that
#' should be employed.
#' @param trace logical value indicating whether progress in estimation should
#' be indicated by printing the number of the permutation that is currently
#' being evaluated.
#' @param \dots miscellaneous parameters for the calls to
#' \code{\link{CoxBoost}}
#' @return Vector with p-value estimates, one value for each optional covariate
#' specificed in the original call to \code{\link{CoxBoost}}.
#' @author Harald Binder \email{binderh@@uni-mainz.de}
#' @seealso \code{\link{CoxBoost}}
#' @references Binder, H., Porzelius, C. and Schumacher, M. (2009). Rank-based
#' p-values for sparse high-dimensional risk prediction models fitted by
#' componentwise boosting. FDM-Preprint Nr. 101, University of Freiburg,
#' Germany.
#' @keywords models regression survial
#' @examples
#'
#' \donttest{
#' #   Generate some survival data with 10 informative covariates
#' n <- 200; p <- 100
#' beta <- c(rep(1,10),rep(0,p-10))
#' x <- matrix(rnorm(n*p),n,p)
#' real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
#' cens.time <- rexp(n,rate=1/10)
#' status <- ifelse(real.time <= cens.time,1,0)
#' obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)
#'
#' #   Fit a Cox proportional hazards model by CoxBoost
#'
#' cbfit <- CoxBoost(time=obs.time,status=status,x=x,stepno=100,
#'                   penalty=100)
#'
#' #   estimate p-values
#'
#' p1 <- estimPVal(cbfit,x,permute.n=10)
#'
#' #   get a second vector of estimates for checking how large
#' #   random variation is
#'
#' p2 <- estimPVal(cbfit,x,permute.n=10)
#'
#' plot(p1,p2,xlim=c(0,1),ylim=c(0,1),xlab="permute 1",ylab="permute 2")
#' }
#'
#' @export
estimPVal <- function(object,x,permute.n=10,per.covariate=FALSE,parallel=FALSE,multicore=FALSE,trace=FALSE,...) {
    if (is.matrix(object$penalty)) {
        stop("Uncertainty cannot be estimated with penalty updates")
    }

    if (object$standardize) x <- scale(x,center=object$meanx,scale=object$sdx)

    permute.index <- matrix(NA,permute.n,nrow(x))
    for (actual.permute in 1:permute.n) permute.index[actual.permute,] <- sample(nrow(x))

    time <- object$time
    status <- object$status
    unpen.index <- object$unpen.index
    if (length(object$unpen.index) > 0) {
        pen.index <- (1:ncol(x))[-object$unpen.index]
    } else {
        pen.index <- 1:ncol(x)
    }
    stepno <- object$stepno
    penalty <- rep(object$penalty,2)

    if (stepno == 0) return(rep(NA,length(pen.index)))

    eval.permute <- function(actual.permute,...) {
        if (trace) cat("permutation",actual.permute,"\n")

        actual.x <- cbind(x,x[permute.index[actual.permute,],pen.index])

        permute.res <- CoxBoost(time=time,status=status,x=actual.x,unpen.index=unpen.index,
                                standardize=FALSE,stepno=stepno,penalty=penalty,trace=FALSE,...)

        actual.score <- colMeans(permute.res$scoremat[,1:length(pen.index),drop=FALSE], na.rm = TRUE)
        null.score <- colMeans(permute.res$scoremat[,(length(pen.index)+1):ncol(permute.res$scoremat),drop=FALSE], na.rm = TRUE)

        if (per.covariate) {
            return(actual.score <= null.score)
        } else {
            return(unlist(lapply(actual.score,function(arg) mean(arg <= null.score))))
        }
    }

    done.parallel <- FALSE

    if (parallel) {
        if (!requireNamespace("snowfall")) {
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            snowfall::sfLibrary(CoxBoost)
            snowfall::sfExport("time","status","x","unpen.index","pen.index","stepno","penalty","trace","permute.index")
            permute.mat <- matrix(unlist(snowfall::sfClusterApplyLB(1:permute.n,eval.permute,...)),length(pen.index),permute.n)
            done.parallel <- TRUE
        }
    }

    if (!done.parallel & multicore) {
        if (!requireNamespace("parallel")) {
            warning("package 'parallel' not found, i.e., parallelization cannot be performed using this package")
        } else {
            if (multicore > 1) {
                permute.mat <- matrix(unlist(parallel::mclapply(1:permute.n,eval.permute,mc.preschedule=FALSE,mc.cores=multicore,...)),length(pen.index),permute.n)
            } else {
                permute.mat <- matrix(unlist(parallel::mclapply(1:permute.n,eval.permute,mc.preschedule=FALSE,...)),length(pen.index),permute.n)
            }
            done.parallel <- TRUE
        }
    }

    if (!done.parallel) {
        permute.mat <- matrix(unlist(lapply(1:permute.n,eval.permute)),length(pen.index),permute.n)
    }

    rowMeans(permute.mat)
}
