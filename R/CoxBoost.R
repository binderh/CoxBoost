weight.at.times <- function(weights,weight.times,times) {
    res <- matrix(0,nrow(weights),length(times))
    for (i in seq(along=times)) {
        res[,i] <- weights[,max(which(weight.times <= times[i]))]
    }

    res
}

efron.weightmat <- function(time,status,cause.code,weights=NULL,prune.times=FALSE,weight.times=NULL) {
    n <- length(time)
    uncens <- which(status == cause.code)
    weightmat <- matrix(0,n,length(uncens))

    rept <- rep(0,n)
    for (i in 1:n) rept[i] <- sum(time[i:n]==time[i] & status[i:n]== cause.code)

    for (i in 1:length(uncens)) {
        weightmat[time >= time[uncens][i] | (status != 0 & status != cause.code), i] <- 1
        tie <- time == time[uncens[i]] & status==cause.code
        di <- max(rept[tie])
        weightmat[tie, i] <- weightmat[tie, i] - (di - rept[uncens[i]])/di
    }

    #   check for competing risks scenario
    if (any(status != 0 & status != cause.code)) {
        cens.ind <- ifelse(status == 0,1,0)
        surv.res <- summary(survival::survfit(Surv(time,cens.ind) ~ 1),times=sort(time))$surv
        invcensprob <- rep(surv.res[length(surv.res)],length(time))
        invcensprob[order(time)[1:length(surv.res)]] <- surv.res

        for (i in 1:length(uncens)) {
            current.invcens <- invcensprob[uncens][i]
            weightmat[,i] <- weightmat[,i] *
                current.invcens/ifelse(time < time[uncens][i],invcensprob,current.invcens)
        }
    }

    if (!is.null(weights)) {
        if (is.matrix(weights) && prune.times) {
            weights <- weight.at.times(weights,weight.times,time[uncens])
        }

        weightmat <- weightmat * weights
    }

    weightmat
}

find.best <- function(x.double.vec,n,p,uncens.C,uncens,
                      actual.beta,actual.risk.score,actual.linear.predictor,
                      weight.double.vec,max.nz.vec,max.1.vec,
                      weightmat.times.risk,weightmat.times.risk.sum,weight.time.dependent,
                      penalty,criterion,actual.step,
                      x.is.01,presel.index=NULL,first.score=NULL,heuristic=TRUE)
{
    if (actual.step == 1 || is.null(presel.index)) {
        actual.presel <- 0:(p-1)
    } else {
        actual.presel <- as.integer(presel.index - 1)
    }

    res <- .C("find_best_candidate",
            x.double.vec,
            as.integer(n),
            as.integer(p),
            uncens.C,
            as.integer(length(uncens)),
            as.double(actual.beta),
            as.double(actual.risk.score),
            as.double(actual.linear.predictor),
            weight.double.vec,
            max.nz.vec,
            max.1.vec,
            as.double(weightmat.times.risk),
            as.double(weightmat.times.risk.sum),
            as.integer(weight.time.dependent),
            as.double(penalty),
            as.integer(criterion == "pscore" || criterion == "hpscore"),
            actual.presel,
            as.integer(length(actual.presel)),
            as.integer(x.is.01),
            warncount=integer(1),
            min.index=integer(1),
            min.beta.delta=double(1),
            score.vec=double(p),
            beta.delta.vec=double(p),
            U.vec=double(p),
            I.vec=double(p),
            NAOK=TRUE,
            PACKAGE="CoxBoost"
            )

    if (heuristic && !is.null(presel.index) && actual.step > 1) {
        min.presel.score <- min(res$score.vec[presel.index])
        if (length(presel.index) < length(first.score) &&
            min.presel.score < max(first.score[-presel.index]))
        {
            new.candidates <- sort(union(which(first.score > min.presel.score),presel.index))

            res <- .C("find_best_candidate",
                    x.double.vec,
                    as.integer(n),
                    as.integer(p),
                    uncens.C,
                    as.integer(length(uncens)),
                    as.double(actual.beta),
                    as.double(actual.risk.score),
                    as.double(actual.linear.predictor),
                    weight.double.vec,
                    max.nz.vec,
                    max.1.vec,
                    as.double(weightmat.times.risk),
                    as.double(weightmat.times.risk.sum),
                    as.integer(weight.time.dependent),
                    as.double(penalty),
                    as.integer(criterion == "pscore"),
                    as.integer(new.candidates - 1),
                    as.integer(length(new.candidates)),
                    as.integer(x.is.01),
                    warncount=integer(1),
                    min.index=integer(1),
                    min.beta.delta=double(1),
                    score.vec=double(p),
                    beta.delta.vec=double(p),
                    U.vec=double(p),
                    I.vec=double(p),
                    NAOK=TRUE,
                    PACKAGE="CoxBoost"
                    )

        }
    }

    res
}

calc.Lambda <- function(event.times,time,uncens,weightmat.times.risk.sum,weights)
{
    actual.Lambda <- rep(NA,length(event.times))
    for (i in seq(along=event.times)) {
        actual.mask <- time[uncens] <= event.times[i]
        if (is.null(weights)) {
            actual.Lambda[i] <- sum(1/weightmat.times.risk.sum[actual.mask])
        } else {
            if (is.matrix(weights)) {
                actual.Lambda[i] <- sum(weights[uncens,i][actual.mask]/weightmat.times.risk.sum[actual.mask])
            } else {
                actual.Lambda[i] <- sum(weights[uncens][actual.mask]/weightmat.times.risk.sum[actual.mask])
            }
        }
    }

    actual.Lambda
}

update.ml.fraction <- function(ml.fraction,weightmat,actual.risk.score,x,subset.time.order,pen.index,min.index,
                               n.uncens,penalty)
{
    actual.x.bar <- apply(weightmat*actual.risk.score*x[subset.time.order,pen.index[min.index]],2,sum)/apply(weightmat*actual.risk.score,2,sum)
    I <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[subset.time.order,pen.index[min.index]],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))
    nu <- I / (I + penalty[min.index])

    ml.fraction + (1-ml.fraction)*nu
}

update.penalty <- function(penalty,sf.scheme,actual.stepsize.factor,
                           min.index,ml.fraction,pendistmat,connected.index,penpos,
                           weightmat,actual.risk.score,x,subset.time.order,pen.index,
                           n.uncens,uncens,n,weightmat.times.risk,weightmat.times.risk.sum,xnames,trace)
{
    if (is.null(pendistmat) || (min.index %in% connected.index && any(pendistmat[penpos[min.index],] != 0))) {
        I.index <- min.index

        actual.x.bar <- apply(weightmat*actual.risk.score*x[subset.time.order,pen.index[min.index]],2,sum)/apply(weightmat*actual.risk.score,2,sum)
        I <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[subset.time.order,pen.index[min.index]],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))

        if (!is.null(pendistmat)) {
            connected <- connected.index[which(pendistmat[penpos[min.index],] != 0)]
            if (length(connected) > 0) {
                change.index <- connected[ml.fraction[connected] < 1]
                I.index <- c(I.index,change.index)
            }
        }

        I.vec <- .C("get_I_vec",
                  as.double(x[subset.time.order,pen.index[I.index]]),
                  as.integer(n),
                  as.integer(length(I.index)),
                  as.integer(length(uncens)),
                  as.double(weightmat.times.risk),
                  as.double(weightmat.times.risk.sum),
                  I.vec=double(length(I.index)),
                  PACKAGE="CoxBoost"
                  )$I.vec

        old.penalty <- penalty[min.index]
        if (sf.scheme == "sigmoid") {
            new.nu <- max(1 - (1-(I.vec[1]/(I.vec[1]+penalty[min.index])))^actual.stepsize.factor,0.00001)  # prevent penalty -> Inf
            penalty[min.index] <- (1/new.nu - 1)*I
        } else {
            penalty[min.index] <- (1/actual.stepsize.factor - 1)*I.vec[1] + penalty[min.index]/actual.stepsize.factor
        }
        if (penalty[min.index] < 0) penalty[min.index] <- 0

        if (length(I.vec) > 1) {
            if (trace) {
                cat("\npenalty update for ",xnames[pen.index][min.index]," (mlf: ",round(ml.fraction[min.index],3),"): ",old.penalty," -> ",penalty[min.index],"\n",sep="")
            }

            change.I <- I.vec[2:length(I.vec)]
            if (trace) cat("adjusting penalty for covariates",paste(xnames[pen.index][change.index],collapse=", ",sep=""),"\n")

            new.target.penalty <- pendistmat[penpos[min.index],penpos[change.index]]*
                                  (1 - ml.fraction[change.index])*change.I/
                                  ((1-actual.stepsize.factor)*pendistmat[penpos[min.index],penpos[change.index]]*
                                                              (1-ml.fraction[min.index])*I.vec[1]/(I.vec[1]+old.penalty) +
                                   (1-ml.fraction[change.index])*change.I/(change.I+penalty[change.index])) -
                                  change.I

            penalty[change.index] <- ifelse(new.target.penalty > 0,new.target.penalty,penalty[change.index])
        }
    }

    penalty
}




#' Fit a Cox model by likelihood based boosting
#'
#' \code{CoxBoost} is used to fit a Cox proportional hazards model by
#' componentwise likelihood based boosting.  It is especially suited for models
#' with a large number of predictors and allows for mandatory covariates with
#' unpenalized parameter estimates.
#'
#' In contrast to gradient boosting (implemented e.g. in the \code{glmboost}
#' routine in the R package \code{mboost}, using the \code{CoxPH} loss
#' function), \code{CoxBoost} is not based on gradients of loss functions, but
#' adapts the offset-based boosting approach from Tutz and Binder (2007) for
#' estimating Cox proportional hazards models. In each boosting step the
#' previous boosting steps are incorporated as an offset in penalized partial
#' likelihood estimation, which is employed for obtain an update for one single
#' parameter, i.e., one covariate, in every boosting step. This results in
#' sparse fits similar to Lasso-like approaches, with many estimated
#' coefficients being zero. The main model complexity parameter, which has to
#' be selected (e.g. by cross-validation using \code{\link{cv.CoxBoost}}), is
#' the number of boosting steps \code{stepno}. The penalty parameter
#' \code{penalty} can be chosen rather coarsely, either by hand or using
#' \code{\link{optimCoxBoostPenalty}}.
#'
#' The advantage of the offset-based approach compared to gradient boosting is
#' that the penalty structure is very flexible. In the present implementation
#' this is used for allowing for unpenalized mandatory covariates, which
#' receive a very fast coefficient build-up in the course of the boosting
#' steps, while the other (optional) covariates are subjected to penalization.
#' For example in a microarray setting, the (many) microarray features would be
#' taken to be optional covariates, and the (few) potential clinical covariates
#' would be taken to be mandatory, by including their indices in
#' \code{unpen.index}.
#'
#' If a group of correlated covariates has influence on the response, e.g.
#' genes from the same pathway, componentwise boosting will often result in a
#' non-zero estimate for only one member of this group. To avoid this,
#' information on the connection between covariates can be provided in
#' \code{pendistmat}. If then, in addition, a penalty updating scheme with
#' \code{stepsize.factor} < 1 is chosen, connected covariates are more likely
#' to be chosen in future boosting steps, if a directly connected covariate has
#' been chosen in an earlier boosting step (see Binder and Schumacher, 2009b).
#'
#' @param time vector of length \code{n} specifying the observed times.
#' @param status censoring indicator, i.e., vector of length \code{n} with
#' entries \code{0} for censored observations and \code{1} for uncensored
#' observations. If this vector contains elements not equal to \code{0} or
#' \code{1}, these are taken to indicate events from a competing risk and a
#' model for the subdistribution hazard with respect to event \code{1} is
#' fitted (see e.g. Fine and Gray, 1999; Binder et al. 2009a).
#' @param x \code{n * p} matrix of covariates.
#' @param unpen.index vector of length \code{p.unpen} with indices of mandatory
#' covariates, where parameter estimation should be performed unpenalized.
#' @param standardize logical value indicating whether covariates should be
#' standardized for estimation. This does not apply for mandatory covariates,
#' i.e., these are not standardized.
#' @param subset a vector specifying a subset of observations to be used in the
#' fitting process.
#' @param weights optional vector of length \code{n}, for specifying weights
#' for the individual observations.
#' @param penalty penalty value for the update of an individual element of the
#' parameter vector in each boosting step.
#' @param criterion indicates the criterion to be used for selection in each
#' boosting step. \code{"pscore"} corresponds to the penalized score
#' statistics, \code{"score"} to the un-penalized score statistics. Different
#' results will only be seen for un-standardized covariates (\code{"pscore"}
#' will result in preferential selection of covariates with larger covariance),
#' or if different penalties are used for different covariates.
#' \code{"hpscore"} and \code{"hscore"} correspond to \code{"pscore"} and
#' \code{"score"}. However, a heuristic is used for evaluating only a subset of
#' covariates in each boosting step, as described in Binder et al. (2011). This
#' can considerably speed up computation, but may lead to different results.
#' @param stepsize.factor determines the step-size modification factor by which
#' the natural step size of boosting steps should be changed after a covariate
#' has been selected in a boosting step. The default (value \code{1}) implies
#' constant penalties, for a value < 1 the penalty for a covariate is increased
#' after it has been selected in a boosting step, and for a value > 1 the
#' penalty it is decreased. If \code{pendistmat} is given, penalty updates are
#' only performed for covariates that have at least one connection to another
#' covariate.
#' @param sf.scheme scheme for changing step sizes (via
#' \code{stepsize.factor}). \code{"linear"} corresponds to the scheme described
#' in Binder and Schumacher (2009b), \code{"sigmoid"} employs a sigmoid shape.
#' @param pendistmat connection matrix with entries ranging between 0 and 1,
#' with entry \code{(i,j)} indicating the certainty of the connection between
#' covariates \code{i} and \code{j}. According to this information penalty
#' changes due to \code{stepsize.factor} < 1 are propagated, i.e., if entry
#' \code{(i,j)} is non-zero, the penalty for covariate \code{j} is decreased
#' after it has been increased for covariate \code{i}, after it has been
#' selected in a boosting step. This matrix either has to have dimension
#' \code{(p - p.unpen) * (p - p.unpen)} or the indicices of the
#' \code{p.connected} connected covariates have to be given in
#' \code{connected.index}, in which case the matrix has to have dimension
#' \code{p.connected * p.connected}. Generally, sparse matices from package
#' \code{Matrix} can be used to save memory.
#' @param connected.index indices of the \code{p.connected} connected
#' covariates, for which \code{pendistmat} provides the connection information
#' for distributing changes in penalties. No overlap with \code{unpen.index} is
#' allowed. If \code{NULL}, and a connection matrix is given, all covariates
#' are assumed to be connected.
#' @param stepno number of boosting steps (\code{m}).
#' @param x.is.01 logical value indicating whether (the non-mandatory part of)
#' \code{x} contains just values 0 and 1, i.e., binary covariates. If this is
#' the case and indicated by this argument, computations are much faster.
#' @param return.score logical value indicating whether the value of the score
#' statistic (or penalized score statistic, depending on \code{criterion}), as
#' evaluated in each boosting step for every covariate, should be returned. The
#' corresponding element \code{scoremat} can become very large (and needs much
#' memory) when the number of covariates and boosting steps is large.
#' @param trace logical value indicating whether progress in estimation should
#' be indicated by printing the name of the covariate updated.
#' @param cmprsk type of competing risk, specific hazards or cause-specific
#' @param coupled.strata logical value indicating whether strata should be
#' coupled during variable selection in each boosting step. If \code{TRUE}
#' (default), the same covariate is selected and updated simultaneously across
#' all strata, assuming proportional effects between strata. If \code{FALSE},
#' strata are treated independently during selection, allowing stratum-specific
#' updates (only relevant when multiple strata are defined and
#' \code{criterion} is set to \code{"hscore"} or \code{"hpscore"}).
#' @param stratum vector specifying different groups of individuals for a
#' stratified Cox regression. In \code{CoxBoost} fit each group gets its own
#' baseline hazard.
#'
#' @return \code{CoxBoost} returns an object of class \code{CoxBoost}.
#'
#' \item{n, p}{number of observations and number of covariates.}
#' \item{stepno}{number of boosting steps.} \item{xnames}{vector of length
#' \code{p} containing the names of the covariates. This information is
#' extracted from \code{x} or names following the scheme \code{V1, V2, ...}}
#' are used. \item{coefficients}{\code{(stepno+1) * p} matrix containing the
#' coefficient estimates for the (standardized) optional covariates for
#' boosting steps \code{0} to \code{stepno}. This will typically be a sparse
#' matrix, built using package \code{\link[Matrix]{Matrix}}.}
#' \item{scoremat}{\code{stepno * p} matrix containing the value of the score
#' statistic for each of the optional covariates before each boosting step.}
#' \item{meanx, sdx}{vector of mean values and standard deviations used for
#' standardizing the covariates.} \item{unpen.index}{indices of the mandatory
#' covariates in the original covariate matrix \code{x}.} \item{penalty}{If
#' \code{stepsize.factor != 1}, \code{stepno * (p - p.unpen)} matrix containing
#' the penalties used for every boosting step and every penalized covariate,
#' otherwise a vector containing the unchanged values of the penalty employed
#' in each boosting step.} \item{time}{observed times given in the
#' \code{CoxBoost} call.} \item{status}{censoring indicator given in the
#' \code{CoxBoost} call.} \item{event.times}{vector with event times from the
#' data given in the \code{CoxBoost} call.}
#' \item{linear.predictors}{\code{(stepno+1) * n} matrix giving the linear
#' predictor for boosting steps \code{0} to \code{stepno} and every
#' observation.} \item{Lambda}{matrix with the Breslow estimate for the
#' cumulative baseline hazard for boosting steps \code{0} to \code{stepno} for
#' every event time.} \item{logplik}{partial log-likelihood of the fitted model
#' in the final boosting step.}
#' @author Written by Harald Binder \email{binderh@@uni-mainz.de}.
#' @seealso \code{\link{predict.CoxBoost}}, \code{\link{cv.CoxBoost}}.
#' @references Binder, H., Benner, A., Bullinger, L., and Schumacher, M.
#' (2013). Tailoring sparse multivariable regression techniques for prognostic
#' single-nucleotide polymorphism signatures. Statistics in Medicine, doi:
#' 10.1002/sim.5490.
#'
#' Binder, H., Allignol, A., Schumacher, M., and Beyersmann, J. (2009).
#' Boosting for high-dimensional time-to-event data with competing risks.
#' Bioinformatics, 25:890-896.
#'
#' Binder, H. and Schumacher, M. (2009). Incorporating pathway information into
#' boosting estimation of high-dimensional risk prediction models. BMC
#' Bioinformatics. 10:18.
#'
#' Binder, H. and Schumacher, M. (2008). Allowing for mandatory covariates in
#' boosting estimation of sparse high-dimensional survival models. BMC
#' Bioinformatics. 9:14.
#'
#' Tutz, G. and Binder, H. (2007) Boosting ridge regression. Computational
#' Statistics & Data Analysis, 51(12):6044-6059.
#'
#' Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the
#' subdistribution of a competing risk. Journal of the American Statistical
#' Association. 94:496-509.
#' @keywords models regression survial
#' @examples
#'
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
#' cbfit <- CoxBoost(time=obs.time,status=status,x=x,stepno=100,penalty=100)
#' summary(cbfit)
#'
#' #   ... with covariates 1 and 2 being mandatory
#'
#' cbfit.mand <- CoxBoost(time=obs.time,status=status,x=x,unpen.index=c(1,2),
#'                        stepno=100,penalty=100)
#' summary(cbfit.mand)
#'
#' @export
CoxBoost <- function(time,status,x,unpen.index=NULL,standardize=TRUE,subset=1:length(time),
                     weights=NULL,stratum=NULL,stepno=100,penalty=9*sum(status[subset]==1),
                     criterion=c("pscore","score","hpscore","hscore"),
                     cmprsk=c("sh","csh","ccsh"),coupled.strata=TRUE,
                     stepsize.factor=1,sf.scheme=c("sigmoid","linear"),pendistmat=NULL,connected.index=NULL,
                     x.is.01=FALSE,return.score=TRUE,trace=FALSE)
{
    sf.scheme <- match.arg(sf.scheme)
    criterion <- match.arg(criterion)
    cmprsk <- match.arg(cmprsk)

    ori.stepno <- stepno
    if (is.list(stepno)) stepno <- max(unlist(stepno))

    if (any(is.na(x))) {
        stop("'x' may not contain missing values")
    }

    if (length(unpen.index) >= ncol(x)) {
        stop("All covariates are indicated as mandatory. At least one non-mandatory covariate is needed.")
    }

    if (!is.null(weights) && is.matrix(weights)) {
        if (is.null(attr(weights,"times"))) {
            unique.times <- sort(unique(time[status != 0]))
            if (ncol(weights) != length(unique.times)) {
                stop("if 'attr(,\"times\")' is not specified for matrix 'weights', the number of columns has to be the number of unique event times.")
            }
            attr(weights,"times") <- unique.times
        }

        weight.times <- attr(weights,"times")
    }

    object <- list()

    #   reduce response to subset
    time <- time[subset]
    status <- status[subset]

    #   reorder observations according to time to speed up computations
    object$time <- time
    object$status <- status
    object$event.times <- sort(unique(object$time[object$status != 0]))

    time.order <- order(time,decreasing=TRUE)
    subset.time.order <- (1:nrow(x))[subset][time.order]
    status <- status[time.order]
    time <- time[time.order]

    if (!is.null(weights)) {
        if (is.matrix(weights)) {
            weights <- weights[subset,]
            attr(weights,"times") <- weight.times
        } else {
            weights <- weights[subset]
        }
    } else {
        weights <- rep(1,length(time))
    }

    object$weights <- weights
    if (is.matrix(weights)) attr(object$weights,"times") <- weight.times

    if (is.matrix(weights)) {
        weights <- weights[time.order,]
        attr(weights,"times") <- weight.times
    } else {
        weights <- weights[time.order]
    }

    if (is.null(stratum)) {
        object$strata <- "1"
        object$stratum <- stratum <- rep(1,length(time))

    } else {
        object$stratum <- stratum <- stratum[subset]
        object$strata <- names(table(stratum))
    }

    stratum <- stratum[time.order]

    object$stepno <- stepno
    object$unpen.index <- unpen.index
    pen.index <- 1:ncol(x)
    if (!is.null(unpen.index)) pen.index <- pen.index[-unpen.index]

    if (is.null(colnames(x))) {
        object$xnames <- paste("V",1:ncol(x),sep="")
    } else {
        object$xnames <- colnames(x)
    }

    if (!is.null(connected.index) && any(connected.index %in% unpen.index)) {
        stop("sets of unpenalized and connected covariates may not overlap")
    }

    if (!is.null(pendistmat) && is.null(connected.index)) {
        if (ncol(pendistmat) == ncol(x) - length(unpen.index)) {
            if (!is.null(unpen.index)) {
                connected.index <- (1:ncol(x))[-unpen.index]
            } else {
                connected.index <- 1:ncol(x)
            }
        } else {
            stop("'connected.index' is missing and cannot be guessed")
        }
    }

    if (!is.null(unpen.index)) {
        if (!is.null(connected.index)) connected.index <- match(connected.index,(1:ncol(x))[-unpen.index])
        unpen.x <- x[subset.time.order,unpen.index,drop=FALSE]
    }

    n <- length(status)
    p <- length(pen.index)

    penpos <- match(1:p,connected.index)

    object$n <- n
    object$p <- p

    object$meanx <- rep(0,length(object$xnames))
    object$sdx <- rep(1,length(object$xnames))
    if (standardize) {
        pen.sdx <- apply(x[subset,pen.index,drop=FALSE],2,sd)
        pen.sdx <- ifelse(pen.sdx == 0,1,pen.sdx)
        pen.meanx <- apply(x[subset,pen.index,drop=FALSE],2,mean)
        x[subset,pen.index] <- scale(x[subset,pen.index],center=pen.meanx,scale=pen.sdx)

        object$meanx[pen.index] <- pen.meanx
        object$sdx[pen.index] <- pen.sdx
        object$standardize <- TRUE
    } else {
        object$standardize <- FALSE
    }

    #   data structures for different causes and different strata

    if (!is.null(unpen.index)) {
        unpen.x.ori <- unpen.x
        unpen.x <- list()
        for (i in seq(along=object$strata)) {
            unpen.x[[i]] <- unpen.x.ori[stratum == object$strata[i],,drop=FALSE]
        }
        unpen.x.ori <- NULL
    }

    x.double.vec <- list()
    for (i in seq(along=object$strata)) {
        x.double.vec[[i]] <- as.double(x[subset.time.order,pen.index,drop=FALSE][stratum == object$strata[i],,drop=FALSE])
    }

    object$causes <- names(table(status[status != 0]))
    if (length(object$causes) > 1 && cmprsk == "sh") {
        object$causes <- object$causes[1]
        if (is.matrix(weights)) {
            warning("time-dependent weights might conflict with IPC weights")
        }
    }

    model <- list()
    for (i in seq(along=object$causes)) {
        for (j in seq(along=object$strata)) {
            model.index <- (i - 1)*length(object$strata)+j

            model[[model.index]] <- list(code=object$causes[i],stratum=object$strata[j])

            actual.smask <- stratum == model[[model.index]]$stratum
            actual.n <- sum(actual.smask)

            model[[model.index]]$n <- actual.n

            if (is.list(penalty)) {
                if (is.null(names(penalty))) {
                    actual.penalty <- penalty[[i]]
                } else {
                    actual.penalty <- penalty[[model[[model.index]]$code]]
                }
            } else {
                actual.penalty <- penalty
            }

            if (is.list(ori.stepno)) {
                if (is.null(names(ori.stepno))) {
                    model[[model.index]]$stepno <- ori.stepno[[i]]
                } else {
                    model[[model.index]]$stepno <- ori.stepno[[model[[model.index]]$code]]
                }
            } else {
                model[[model.index]]$stepno <- ori.stepno
            }

            if (length(actual.penalty) < length(pen.index)) model[[model.index]]$penalty <- rep(actual.penalty[1],length(pen.index))

            if (any(stepsize.factor != 1)) {
                model[[model.index]]$penaltymat <- matrix(NA,stepno,length(model[[model.index]]$penalty))
            }

            model[[model.index]]$reverse.time.order <- match(seq(along=object$time)[object$stratum == object$strata[j]],time.order[actual.smask])

            model[[model.index]]$uncens <- which(status[actual.smask] == model[[model.index]]$code)
            model[[model.index]]$n.uncens <- length(model[[model.index]]$uncens)
            model[[model.index]]$event.times <- sort(unique(time[actual.smask][model[[model.index]]$uncens]))

            if (is.matrix(weights)) {
                model[[model.index]]$weight.at.times <- weight.at.times(weights[actual.smask,,drop=FALSE],weight.times,time[actual.smask][model[[model.index]]$uncens])

                model[[model.index]]$weight.at.event <- diag(model[[model.index]]$weight.at.times[model[[model.index]]$uncens,,drop=FALSE])
            } else {
                if (!is.null(weights)) {
                    model[[model.index]]$weight.at.event <- weights[actual.smask][model[[model.index]]$uncens]
                } else {
                    model[[model.index]]$weight.at.event <- rep(1.0,length(model[[model.index]]$uncens))
                }
            }

            model[[model.index]]$coefficients <- Matrix::Matrix(0,stepno+1,p)
            if (!is.null(unpen.index)) {
                model[[model.index]]$unpen.coefficients <- matrix(NA,stepno+1,ncol(unpen.x[[1]]))
            } else {
                model[[model.index]]$unpen.coefficients <- NULL
            }
            model[[model.index]]$linear.predictor <- matrix(NA,stepno+1,actual.n)

            model[[model.index]]$Lambda <- matrix(NA,stepno+1,length(model[[model.index]]$event.times))
            if (return.score) model[[model.index]]$scoremat <- matrix(NA,max(1,stepno),object$p)

            #   Efron handling of ties
            if (cmprsk == "sh") {
                if (is.matrix(weights)) {
                    model[[model.index]]$weightmat <- efron.weightmat(time[actual.smask],status[actual.smask],model[[model.index]]$code,model[[model.index]]$weight.at.times)
                } else {
                    model[[model.index]]$weightmat <- efron.weightmat(time[actual.smask],status[actual.smask],model[[model.index]]$code,weights[actual.smask])
                }
            } else {
                if (is.matrix(weights)) {
                    model[[model.index]]$weightmat <- efron.weightmat(time[actual.smask],ifelse(status == model[[model.index]]$code,status,0)[actual.smask],model[[model.index]]$code,model[[model.index]]$weight.at.times)
                } else {
                    model[[model.index]]$weightmat <- efron.weightmat(time[actual.smask],ifelse(status == model[[model.index]]$code,status,0)[actual.smask],model[[model.index]]$code,weights[actual.smask])
                }
            }
            model[[model.index]]$weightmat <- model[[model.index]]$weightmat

            model[[model.index]]$actual.beta <- rep(0,p)
            if (!is.null(unpen.index)) model[[model.index]]$actual.unpen.beta <- rep(0,ncol(unpen.x[[1]]))
            model[[model.index]]$actual.linear.predictor <- rep(0,actual.n)
            model[[model.index]]$actual.risk.score <- rep(1,actual.n)
            model[[model.index]]$ml.fraction <- rep(0,p)

            model[[model.index]]$weight.double.vec <- as.double(model[[model.index]]$weightmat)
            model[[model.index]]$max.nz.vec <- as.integer(apply(model[[model.index]]$weightmat,2,function(arg) max(which(arg != 0))))
            model[[model.index]]$max.1.vec <- as.integer(c(0,rev(cummin(rev(apply(model[[model.index]]$weightmat,2,function(arg) ifelse(!any(arg != 1),length(arg),min(which(arg != 1)-1))))))))
            model[[model.index]]$uncens.C <- as.integer(model[[model.index]]$uncens - 1)

            model[[model.index]]$warnstep <- NULL
            model[[model.index]]$unpen.warn <- NULL

            model[[model.index]]$first.score <- NULL
            model[[model.index]]$presel.index <- c()
        }
    }

    #   boosting iterations

    for (actual.step in 0:stepno) {
        model.score <- NULL
        model.beta.delta <- NULL
        model.U <- NULL
        model.I <- NULL

        for (cause.index in seq(along=object$causes)) {
            for (stratum.index in seq(along=object$strata)) {
                model.index <- (cause.index - 1)*length(object$strata)+stratum.index
                actual.smask <- stratum == object$strata[stratum.index]

                if (actual.step > 0 && any(stepsize.factor != 1)) {
                    model[[model.index]]$penaltymat[stepno,] <- model[[model.index]]$penalty
                }

                weightmat.times.risk <- model[[model.index]]$weightmat*model[[model.index]]$actual.risk.score
                weightmat.times.risk.sum <- colSums(weightmat.times.risk)

                #   update unpenalized covariates by one estimation step

                if (!is.null(unpen.index)) {
                    if (actual.step == 1 || !is.null(model[[model.index]]$unpen.warn)) {
                        model[[model.index]]$unpen.coefficients[actual.step+1,] <- model[[model.index]]$actual.unpen.beta
                    } else {
                        x.bar <- (t(weightmat.times.risk) %*% unpen.x[[stratum.index]]) / weightmat.times.risk.sum
                        U <- colSums((unpen.x[[stratum.index]][model[[model.index]]$uncens,] - x.bar)*model[[model.index]]$weight.at.event)

                        I <- matrix(0,ncol(unpen.x[[stratum.index]]),ncol(unpen.x[[stratum.index]]))
                        for (i in 1:model[[model.index]]$n.uncens) {
                            x.minus.bar <- t(t(unpen.x[[stratum.index]]) - x.bar[i,])
                            I <- I + (t(x.minus.bar*(weightmat.times.risk[,i])) %*% x.minus.bar)/weightmat.times.risk.sum[i]*model[[model.index]]$weight.at.event[i]
                        }

                        try.res <- try(unpen.beta.delta <- drop(solve(I) %*% U),silent=TRUE)
                        if (inherits(try.res, "try-error")) {
                            model[[model.index]]$unpen.warn <- actual.step
                            if (actual.step == 0) {
                                model[[model.index]]$unpen.coefficients[actual.step+1,] <- 0
                                model[[model.index]]$actual.unpen.beta <- model[[model.index]]$unpen.coefficients[actual.step+1,]
                            }
                        } else {
                            model[[model.index]]$actual.unpen.beta <- model[[model.index]]$actual.unpen.beta + unpen.beta.delta
                            model[[model.index]]$unpen.coefficients[actual.step+1,] <- model[[model.index]]$actual.unpen.beta

                            model[[model.index]]$actual.linear.predictor <- model[[model.index]]$actual.linear.predictor + drop(unpen.x[[stratum.index]] %*% unpen.beta.delta)
                            model[[model.index]]$actual.risk.score <- exp(drop(model[[model.index]]$actual.linear.predictor))
                            weightmat.times.risk <- model[[model.index]]$weightmat*model[[model.index]]$actual.risk.score
                            weightmat.times.risk.sum <- colSums(weightmat.times.risk)
                        }
                    }
                }

                if (actual.step == 0) {
                    model[[model.index]]$coefficients[actual.step+1,] <- model[[model.index]]$actual.beta
                    model[[model.index]]$linear.predictor[actual.step+1,] <- model[[model.index]]$actual.linear.predictor[model[[model.index]]$reverse.time.order]
                    if (is.matrix(weights)) {
                        model[[model.index]]$Lambda[actual.step+1,] <- calc.Lambda(model[[model.index]]$event.times,time[actual.smask],model[[model.index]]$uncens,weightmat.times.risk.sum,model[[model.index]]$weight.at.times)
                    } else {
                        model[[model.index]]$Lambda[actual.step+1,] <- calc.Lambda(model[[model.index]]$event.times,time[actual.smask],model[[model.index]]$uncens,weightmat.times.risk.sum,
                                                         weights[actual.smask])
                    }
                    next
                }

                res <- find.best(x.double.vec[[stratum.index]],model[[model.index]]$n,p,model[[model.index]]$uncens.C,model[[model.index]]$uncens,
                          model[[model.index]]$actual.beta,model[[model.index]]$actual.risk.score,model[[model.index]]$actual.linear.predictor,
                          model[[model.index]]$weight.at.event,model[[model.index]]$max.nz.vec,model[[model.index]]$max.1.vec,
                          weightmat.times.risk,weightmat.times.risk.sum,is.matrix(weights),
                          model[[model.index]]$penalty,criterion,actual.step,
                          x.is.01,model[[model.index]]$presel.index,model[[model.index]]$first.score,
                          heuristic=!(length(object$strata) > 1 && !coupled.strata && (criterion == "hscore" || criterion == "hpscore")))

                if (is.null(model[[model.index]]$warnstep) && res$warncount > 0) model[[model.index]]$warnstep <- actual.step

                if (return.score) model[[model.index]]$scoremat[actual.step,] <- res$score.vec
                if ((criterion == "hscore" || criterion == "hpscore") && actual.step == 1) {
                    model[[model.index]]$first.score <- res$score.vec
                }

                model.score <- rbind(model.score,res$score.vec)
                model.beta.delta <- rbind(model.beta.delta,res$beta.delta.vec)
                model.U <- rbind(model.U,res$U.vec)
                model.I <- rbind(model.I,res$I.vec)
            }
        }

        if (actual.step == 0) next

        if (length(object$strata) > 1 || cmprsk == "ccsh") {
            if (cmprsk == "ccsh") {
                cause.min.index <- rep(which.max(apply(model.score,2,function(arg) -2*sum(log(1-pchisq(arg,df=1))))),length(object$causes))
            } else {
                cause.min.index <- integer(length(object$causes))
                for (i in seq(along=object$causes)) {
                    if (coupled.strata) {
                        cause.min.index[i] <- which.max(apply(model.score[((i-1)*length(object$strata)+1):(i*length(object$strata)),,drop=FALSE],2,function(arg) -2*sum(log(1-pchisq(arg,df=1)))))
                    } else {
                        actual.first <- ((i-1)*length(object$strata)+1)
                        actual.last <- (i*length(object$strata))

                        stratum.U <- colSums(model.U[actual.first:actual.last,,drop=FALSE])
                        stratum.I <- colSums(model.I[actual.first:actual.last,,drop=FALSE])

                        if (criterion == "hpscore" || criterion == "pscore") {
                            stratum.score <- stratum.U*stratum.U/(stratum.I + model[[actual.first]]$penalty)
                        } else {
                            stratum.score <- stratum.U*stratum.U/(stratum.I + 1.0/0.1 - 1)
                        }

                        if (return.score) {
                            for (j in actual.first:actual.last) model[[j]]$scoremat[actual.step,] <- stratum.score
                        }

                        if (criterion == "hscore" || criterion == "hpscore") {

                            if (actual.step == 1) {
                                for (j in actual.first:actual.last) {
                                    model[[j]]$first.score <- stratum.score
                                }
                            }

                            if (actual.step > 1) {
                                min.presel.score <- min(stratum.score[model[[actual.first]]$presel.index])

                                if (length(model[[actual.first]]$presel.index) < length(model[[actual.first]]$first.score) &&
                                    min.presel.score < max(model[[actual.first]]$first.score[-model[[actual.first]]$presel.index])) {

                                    new.candidates <- sort(union(which(model[[actual.first]]$first.score > min.presel.score),model[[actual.first]]$presel.index))

                                    for (stratum.index in 1:length(object$strata)) {
                                        model.index <- ((i-1)*length(object$strata)+stratum.index)

                                        weightmat.times.risk <- model[[model.index]]$weightmat*model[[model.index]]$actual.risk.score
                                        weightmat.times.risk.sum <- colSums(weightmat.times.risk)

                                        res <- find.best(x.double.vec[[stratum.index]],model[[model.index]]$n,p,
                                                         model[[model.index]]$uncens.C,model[[model.index]]$uncens,
                                                         model[[model.index]]$actual.beta,
                                                         model[[model.index]]$actual.risk.score,
                                                         model[[model.index]]$actual.linear.predictor,
                                                         model[[model.index]]$weight.at.event,
                                                         model[[model.index]]$max.nz.vec,model[[model.index]]$max.1.vec,
                                                         weightmat.times.risk,weightmat.times.risk.sum,is.matrix(weights),
                                                         model[[model.index]]$penalty,criterion,actual.step,
                                                         x.is.01,new.candidates,
                                                         model[[model.index]]$first.score,heuristic=FALSE)


                                        if (is.null(model[[model.index]]$warnstep) && res$warncount > 0) {
                                            model[[model.index]]$warnstep <- actual.step
                                        }


                                        model.score[model.index,] <- res$score.vec
                                        model.beta.delta[model.index,] <- res$beta.delta.vec
                                        model.U[model.index,] <- res$U.vec
                                        model.I[model.index,] <- res$I.vec
                                    }

                                    stratum.U <- colSums(model.U[actual.first:actual.last,,drop=FALSE])
                                    stratum.I <- colSums(model.I[actual.first:actual.last,,drop=FALSE])

                                    if (criterion == "hpscore" || criterion == "pscore") {
                                        stratum.score <- stratum.U*stratum.U/(stratum.I + model[[actual.first]]$penalty)
                                    } else {
                                        stratum.score <- stratum.U*stratum.U/(stratum.I + 1.0/0.1 - 1)
                                    }

                                    if (return.score) {
                                        for (j in actual.first:actual.last) model[[j]]$scoremat[actual.step,] <- stratum.score
                                    }
                                }
                            }
                        }

                        candidate.indices <- (1:length(stratum.I))[stratum.I != 0]
                        candidate.penalty <- model[[actual.first]]$penalty[stratum.I != 0]
                        stratum.U <- stratum.U[stratum.I != 0]
                        stratum.I <- stratum.I[stratum.I != 0]

                        if (criterion == "pscore" || criterion == "hpscore") {
                            candidate.min <- which.max(stratum.U*stratum.U/(stratum.I + candidate.penalty))
                            cause.min.index[i] <- candidate.indices[candidate.min]
                        } else {
                            candidate.min <- which.max(stratum.U*stratum.U/(stratum.I + 1.0/0.1 - 1))
                            cause.min.index[i] <- candidate.indices[candidate.min]
                        }

                        model.beta.delta[actual.first:actual.last,cause.min.index[i]] <- stratum.U[candidate.min]/(stratum.I[candidate.min] + model[[actual.first]]$penalty[cause.min.index[i]])
                    }
                }
            }
        } else {
            cause.min.index <- apply(model.score,1,which.max)
        }

        for (cause.index in seq(along=object$causes)) {
            min.index <- cause.min.index[cause.index]
            if (trace) cat(object$xnames[pen.index][min.index]," ",sep="")

            for (stratum.index in seq(along=object$strata)) {
                model.index <- (cause.index - 1)*length(object$strata)+stratum.index
                actual.smask <- stratum == object$strata[stratum.index]
                min.beta.delta <- model.beta.delta[model.index,min.index]

                if (criterion == "hscore" || criterion == "hpscore") {
                    model[[model.index]]$presel.index <- sort(union(model[[model.index]]$presel.index,min.index))
                }

                if (!is.null(pendistmat)) {
                    model[[model.index]]$ml.fraction[min.index] <- update.ml.fraction(model[[model.index]]$ml.fraction[min.index],
                                                         model[[model.index]]$weightmat,model[[model.index]]$actual.risk.score,x,
                                                         subset.time.order[stratum == object$strata[stratum.index]],pen.index,min.index,
                                                         model[[model.index]]$n.uncens,model[[model.index]]$penalty)
                }

                model[[model.index]]$actual.beta[min.index] <- model[[model.index]]$actual.beta[min.index] + min.beta.delta

                if (length(object$strata) > 1) {
                    model[[model.index]]$actual.linear.predictor <- model[[model.index]]$actual.linear.predictor + x[subset.time.order[stratum == object$strata[stratum.index]],pen.index[min.index],drop=FALSE]*min.beta.delta
                } else {
                    model[[model.index]]$actual.linear.predictor <- model[[model.index]]$actual.linear.predictor + x[subset.time.order,pen.index[min.index]]*min.beta.delta
                }

                model[[model.index]]$actual.risk.score <- exp(drop(model[[model.index]]$actual.linear.predictor))
                weightmat.times.risk <- model[[model.index]]$weightmat*model[[model.index]]$actual.risk.score
                weightmat.times.risk.sum <- colSums(weightmat.times.risk)

                model[[model.index]]$coefficients[actual.step+1,] <- model[[model.index]]$actual.beta
                model[[model.index]]$linear.predictor[actual.step+1,] <- model[[model.index]]$actual.linear.predictor[model[[model.index]]$reverse.time.order]
                if (is.matrix(weights)) {
                    model[[model.index]]$Lambda[actual.step+1,] <- calc.Lambda(model[[model.index]]$event.times,time[actual.smask],model[[model.index]]$uncens,weightmat.times.risk.sum,model[[model.index]]$weight.at.times)
                } else {
                    model[[model.index]]$Lambda[actual.step+1,] <- calc.Lambda(model[[model.index]]$event.times,time[actual.smask],model[[model.index]]$uncens,weightmat.times.risk.sum,
                                                     weights[actual.smask])
                }

                #   update the penalty if the user has chosen any other value than the default

                actual.stepsize.factor <- ifelse(length(stepsize.factor) >= min.index,stepsize.factor[min.index],stepsize.factor[1])

                if (actual.stepsize.factor != 1 && model[[cause.index]]$ml.fraction[min.index] < 1) {
                    model[[model.index]]$penalty <- update.penalty(model[[model.index]]$penalty,sf.scheme,actual.stepsize.factor,
                           min.index,model[[model.index]]$ml.fraction,pendistmat,connected.index,penpos,
                           model[[model.index]]$weightmat,model[[model.index]]$actual.risk.score,x,subset.time.order[stratum == object$strata[stratum.index]],pen.index,
                           model[[model.index]]$n.uncens,model[[model.index]]$uncens,n,weightmat.times.risk,weightmat.times.risk.sum,object$xnames,trace)
                }
            }
        }
    }
    if (trace) cat("\n")

    warnsteps <- unlist(lapply(model,function(arg) arg$warnstep))
    if (length(warnsteps) > 0) warning(paste("potentially attempted to move towards a nonexisting maximum likelihood solution around step",min(warnsteps)))

    unpen.warns <- unlist(lapply(model,function(arg) arg$unpen.warn))
    if (length(unpen.warns) > 0) warning(paste("estimation for unpenalized covariates did not converge starting at step ",min(unpen.warns),". Values were kept fixed and might be unreliable",sep=""))

    #   combine penalized and unpenalized covariates
    if (!is.null(object$unpen.index)) {
        object$p <- object$p + length(object$unpen.index)

        for (i in seq(along=model)) {
            combined.coefficients <- Matrix::Matrix(0,nrow(model[[i]]$coefficients),object$p)
            combined.coefficients[,pen.index] <- model[[i]]$coefficients
            combined.coefficients[,object$unpen.index] <- model[[i]]$unpen.coefficients
            model[[i]]$coefficients <- combined.coefficients
        }
    }

    #   final assembly of CoxBoost object

    class(object) <- "CoxBoost"

    #   add results directly to the CoxBoost objects if there is only a single cause

    if (length(model) == 1) {
        if (stepsize.factor != 0) {
            object$penalty <- model[[1]]$penaltymat
        } else {
            object$penalty <- model[[1]]$penalty
        }

        object$event.times <- model[[1]]$event.times
        object$coefficients <- model[[1]]$coefficients
        object$linear.predictor <- model[[1]]$linear.predictor
        object$Lambda <- model[[1]]$Lambda
        if (return.score) object$scoremat <- model[[1]]$scoremat

        object$logplik <- predict(object,type="logplik")
    } else {
        object$model <- list()
        for (i in seq(along=model)) {
            object$model[[i]] <- list(cause.code=model[[i]]$code,stratum=model[[i]]$stratum)
            object$model[[i]]$stepno <- model[[i]]$stepno

            if (stepsize.factor != 0) {
                object$model[[i]]$penalty <- model[[i]]$penaltymat
            } else {
                object$model[[i]]$penalty <- model[[i]]$penalty
            }

            object$model[[i]]$event.times <- model[[i]]$event.times
            object$model[[i]]$weight.at.times <- model[[i]]$weight.at.times
            object$model[[i]]$coefficients <- model[[i]]$coefficients
            object$model[[i]]$linear.predictor <- model[[i]]$linear.predictor
            object$model[[i]]$Lambda <- model[[i]]$Lambda
            if (return.score) object$model[[i]]$scoremat <- model[[i]]$scoremat
        }

        names(object$model) <- object$causes

        actual.logplik <- predict(object,type="logplik")
        for (cause.index in seq(along=object$causes)) {
            for (stratum.index in seq(along=object$strata)) {
                object$model[[(cause.index-1)*length(object$strata)+stratum.index]]$logplik <- actual.logplik[[cause.index]]
            }
        }
    }

    object
}

joint.print <- function(x,long=FALSE) {
    for (i in 1:(length(x$causes)*length(x$strata))) {
        if (length(x$causes) > 1 || length(x$strata) > 1) {
            if (length(x$causes) > 1) cat("cause '",x$model[[i]]$cause.code,"\'",sep="")
            if (length(x$causes) > 1 && length(x$strata) > 1) cat(" ")
            if (length(x$strata) > 1) cat("stratum '",x$model[[i]]$stratum,"\'",sep="")
            cat(":\n")
            actual.stepno <- x$model[[i]]$stepno
            actual.coef <- x$model[[i]]$coefficients[x$model[[i]]$stepno+1,]
            actual.logplik <- x$model[[i]]$logplik
        } else {
            actual.stepno <- x$stepno
            actual.coef <- x$coefficients[x$stepno+1,]
            actual.logplik <- x$logplik
        }

        cat(actual.stepno,"boosting steps resulting in",sum(actual.coef != 0),"non-zero coefficients",
            ifelse(is.null(x$unpen.index),"",paste("(with",length(x$unpen.index),"being mandatory)")),
        "\n")
        cat("partial log-likelihood:",actual.logplik,"\n")

        if (long) {
            cat("\n")
            cat("parameter estimates > 0:\n",paste(x$xnames[actual.coef > 0],collapse=", "),"\n")
            cat("parameter estimates < 0:\n",paste(x$xnames[actual.coef < 0],collapse=", "),"\n")

            if (!is.null(x$unpen.index)) {
                cat("Parameter estimates for mandatory covariates:\n",sep="")
                print(matrix(signif(actual.coef[x$unpen.index],4),length(x$unpen.index),1,dimnames=list(x$xnames[x$unpen.index],c("Estimate"))))
            }
        }

        if (i < (length(x$causes)*length(x$strata))) cat("\n")
    }
}

#' Print a CoxBoost type
#' @param x a CoxBoost object
#' @param ... not used
#' @method print CoxBoost
#' @return
#' Invisibly returns \code{NULL}. The function is called for its
#' side-effects, printing for each cause–stratum combination:
#' \itemize{
#'   \item the number of boosting steps,
#'   \item the number of non-zero coefficients,
#'   \item the partial log-likelihood.
#' }
#' @export
print.CoxBoost <- function(x,...) {
    joint.print(x,long=FALSE)
}

#' Summary of a CoxBoost type
#' @param object a CoxBoost object
#' @param ... not used
#' @method summary CoxBoost
#' @return
#' Invisibly returns \code{NULL}. The function is called for its
#' side-effects, printing for each cause–stratum combination:
#' \itemize{
#'   \item the number of boosting steps and non-zero coefficients,
#'   \item the partial log-likelihood,
#'   \item the covariates with positive coefficients,
#'   \item the covariates with negative coefficients,
#'   \item and, if applicable, the coefficient estimates of mandatory
#'         (unpenalized) covariates.
#' }
#' @export
summary.CoxBoost <- function(object,...) {
    joint.print(object,long=TRUE)
}



#' Plot coefficient paths from CoxBoost fit
#'
#' Plots coefficient paths, i.e. the parameter estimates plotted against the
#' boosting steps as obtained from a CoxBoost object fitted by
#' \code{\link{CoxBoost}}.
#'
#'
#' @param x fitted CoxBoost object from a \code{\link{CoxBoost}} call.
#' @param line.col color of the lines of the coefficient path
#' @param label.cex scaling factor for the variable labels.
#' @param xlab label for the x axis, default label when \code{NULL}.
#' @param ylab label for the y axis, default label when \code{NULL}.
#' @param xlim,ylim plotting ranges, default label when \code{NULL}.
#' @param main a main title for the plot
#' @param \dots miscellaneous arguments, passed to the plot routine.
#' @return No value is returned, but a plot is generated.
#' @author Harald Binder \email{binderh@@uni-mainz.de}
#' @export
plot.CoxBoost <- function(x,line.col="dark grey",label.cex=0.6,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,main=NULL,...) {
    if (is.null(xlab)) xlab <- "boosting step"
    if (is.null(ylab)) ylab <- "estimated coefficients"

    plotmat <- list()
    plot.names <- list()
    nz.index <- list()

    for (i in 1:(length(x$causes)*length(x$strata))) {
        if (((length(x$causes) > 1 || length(x$strata) > 1) && x$model[[i]]$stepno == 0) ||
            (length(x$causes) == 1 && x$stepno == 0))
        {
            plotmat[[i]] <- plot.names[[i]] <- nz.index[[i]] <- c()
            next
        }

        if (length(x$causes) > 1 || length(x$strata) > 1) {
            nz.index[[i]] <- which(Matrix::colSums(abs(x$model[[i]]$coefficients)) > 0)
        } else {
            nz.index[[i]] <- which(Matrix::colSums(abs(x$coefficients)) > 0)
        }
        nz.index[[i]] <- nz.index[[i]][!(nz.index[[i]] %in% x$unpen.index)]
        plot.names[[i]] <- x$xnames[nz.index[[i]]]

        if (length(x$causes) > 1 || length(x$strata) > 1) {
            plotmat[[i]] <- as.matrix(x$model[[i]]$coefficients[,nz.index[[i]],drop=FALSE])[1:(x$model[[i]]$stepno+1),,drop=FALSE]
            ncoef <- ncol(x$model[[i]]$coefficients)
        } else {
            plotmat[[i]] <- as.matrix(x$coefficients[,nz.index[[i]],drop=FALSE])
            ncoef <- ncol(x$coefficients)
        }
    }

    if (is.null(ylim)) ylim <- c(min(unlist(lapply(plotmat,min))),max(unlist(lapply(plotmat,max))))

    for (cause.index in seq(along=x$causes)) {
        for (stratum.index in seq(along=x$strata)) {
            i <- (cause.index-1)*length(x$strata)+stratum.index

            actual.main <- main
            if (is.null(actual.main)) {
                if (length(x$causes) > 1) actual.main <- paste("cause '",x$causes[[cause.index]],"'",sep="")
                if (length(x$causes) > 1 && length(x$strata) > 1) actual.main <- paste(actual.main," ",sep="")
                if (length(x$strata) > 1) actual.main <- paste(actual.main,"stratum '",x$strata[[stratum.index]],"'",sep="")
            }

            if (length(nz.index[[i]]) == 0) {
                plot(1,type="n",xlab=xlab,ylab=ylab,ylim=ylim,main=actual.main,...)
                next
            }

            actual.xlim <- xlim

            if (is.null(actual.xlim)) actual.xlim <- c(0,(nrow(plotmat[[i]])-1)*(1+0.017*max(nchar(plot.names[[i]]))))

            plot(1,type="n",xlim=actual.xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=actual.main,...)

            if (length(nz.index[[i]]) < ncoef - length(x$unpen.index)) graphics::lines(c(0,nrow(plotmat[[i]])-1),c(0,0),col=line.col)

            for (coef.index in 1:ncol(plotmat[[i]])) {
              graphics::lines(0:(nrow(plotmat[[i]])-1),plotmat[[i]][,coef.index],col=line.col)
              graphics::text(actual.xlim[2],plotmat[[i]][nrow(plotmat[[i]]),coef.index],plot.names[[i]][coef.index],pos=2,cex=label.cex)
            }
        }
    }
}



#' Coeffients from CoxBoost fit
#'
#' Extracts the coefficients from the specified boosting steps of a CoxBoost
#' object fitted by \code{\link{CoxBoost}}.
#'
#'
#' @param object fitted CoxBoost object from a \code{\link{CoxBoost}} call.
#' @param at.step scalar or vector of boosting step(s) at which prediction is
#' wanted. If no step is given, the final boosting step is used.
#' @param scaled logical value indicating whether coefficients should be
#' returned scaled to be at the level of the original covariate values, or
#' unscaled as used internally when \code{standardize=TRUE} is used in the
#' \code{\link{CoxBoost}} call.
#' @param \dots miscellaneous arguments, none of which is used at the moment.
#' @return For a vector of length \code{p} (number of covariates)
#' (\code{at.step} being a scalar) or a \code{length(at.step) * p} matrix
#' (\code{at.step} being a vector).
#' @author Harald Binder \email{binderh@@uni-mainz.de}
#' @export
coef.CoxBoost <- function(object,at.step=NULL,scaled=TRUE,...) {
    beta <- list()

    for (model.index in 1:(length(object$causes)*length(object$strata))) {
        if (is.null(at.step) || length(at.step) == 1) {
            if (length(object$causes) > 1 || length(object$strata) > 1) {
                if (is.null(at.step)) {
                    actual.beta <- object$model[[model.index]]$coefficients[nrow(object$model[[model.index]]$coefficients),]
                } else {
                    actual.beta <- object$model[[model.index]]$coefficients[at.step+1,]
                }
            } else {
                if (is.null(at.step)) {
                    actual.beta <- object$coefficients[nrow(object$coefficients),]
                } else {
                    actual.beta <- object$coefficients[at.step+1,]
                }
            }

            if (scaled) actual.beta <- actual.beta * object$sdx
            names(actual.beta) <- object$xnames

        } else {
            if (length(object$causes) > 1 || length(object$strata) > 1) {
                            actual.beta <- as.matrix(object$model[[model.index]]$coefficients[at.step+1,])
            } else {
                actual.beta <- as.matrix(object$coefficients[at.step+1,])
            }

            if (scaled) actual.beta <- t(t(actual.beta) * object$sdx)
            colnames(actual.beta) <- object$xnames
        }

        if (length(object$causes) > 1 || length(object$strata) > 1) {
            beta[[model.index]] <- actual.beta
        } else {
            return(actual.beta)
        }
    }

    if (length(object$causes) > 1 && length(object$strata) > 1) {
            names(beta) <- paste(rep(object$causes,rep(length(object$strata),length(object$causes))),rep(object$strata,length(object$causes)))
    } else {
        if (length(object$causes) > 1) {
            names(beta) <- object$causes
        } else {
            names(beta) <- object$strata
        }
    }

    beta
}




#' Predict method for CoxBoost fits
#'
#' Obtains predictions at specified boosting steps from a CoxBoost object
#' fitted by \code{\link{CoxBoost}}.
#'
#'
#' @param object fitted CoxBoost object from a \code{\link{CoxBoost}} call.
#' @param newdata \code{n.new * p} matrix with new covariate values. If just
#' prediction for the training data is wanted, it can be omitted.
#' @param newtime,newstatus vectors with observed time and censoring indicator
#' (0 for censoring, 1 for no censoring, and any other values for competing
#' events in a competing risks setting) for new observations, where prediction
#' is wanted. Only required if predicted partial log-likelihood is wanted,
#' i.e., if \code{type="logplik"}. This can also be omitted when prediction is
#' only wanted for the training data, i.e., \code{newdata=NULL}.
#' @param subset an optional vector specifying a subset of observations to be
#' used for evaluation.
#' @param at.step scalar or vector of boosting step(s) at which prediction is
#' wanted. If \code{type="risk"} is used, only one step is admissible. If no
#' step is given, the final boosting step is used.
#' @param times vector with \code{T} time points where prediction is wanted.
#' Only needed for \code{type="risk"}
#' @param type type of prediction to be returned: \code{"lp"} gives the linear
#' predictor, \code{"logplik"} the partial log-likelihood, \code{"risk"} the
#' predicted probability of not yet having had the event at the time points
#' given in \code{times}, and \code{"CIF"} the predicted cumulative incidence
#' function, i.e., the predicted probability of having had the event of
#' interest.
#' @param weights weights for each time.
#' @param stratum vector specifying different groups of individuals for a
#' stratified Cox regression. In \code{CoxBoost} fit each group gets its own
#' baseline hazard.
#' @param \dots miscellaneous arguments, none of which is used at the moment.
#' @return For \code{type="lp"} and \code{type="logplik"} a vector of length
#' \code{n.new} (\code{at.step} being a scalar) or a \code{n.new *
#' length(at.step)} matrix (\code{at.step} being a vector) with predictions is
#' returned.  For \code{type="risk"} or \code{type="CIF"} a \code{n.new * T}
#' matrix with predicted probabilities at the specific time points is returned.
#' @author Harald Binder \email{binderh@@uni-mainz.de}
#' @keywords models regression survial
#' @examples
#'
#' #   Generate some survival data with 10 informative covariates
#' n <- 200; p <- 100
#' beta <- c(rep(1,10),rep(0,p-10))
#' x <- matrix(rnorm(n*p),n,p)
#' real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
#' cens.time <- rexp(n,rate=1/10)
#' status <- ifelse(real.time <= cens.time,1,0)
#' obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)
#'
#' #   define training and test set
#'
#' train.index <- 1:100
#' test.index <- 101:200
#'
#' #   Fit CoxBoost to the training data
#'
#' cbfit <- CoxBoost(time=obs.time[train.index],status=status[train.index],
#'                   x=x[train.index,],stepno=300,penalty=100)
#'
#' #   mean partial log-likelihood for test set in every boosting step
#'
#' step.logplik <- predict(cbfit,newdata=x[test.index,],
#'                         newtime=obs.time[test.index],
#'                         newstatus=status[test.index],
#'                         at.step=0:300,type="logplik")
#'
#' plot(step.logplik)
#'
#' #   names of covariates with non-zero coefficients at boosting step
#' #   with maximal test set partial log-likelihood
#'
#' print(cbfit$xnames[cbfit$coefficients[which.max(step.logplik),] != 0])
#'
#'
#'
#' @export
predict.CoxBoost <- function(object,newdata=NULL,newtime=NULL,newstatus=NULL,subset=NULL,weights=NULL,stratum=NULL,at.step=NULL,times=NULL,type=c("lp","logplik","risk","CIF"),...)
{
    type <- match.arg(type)

    if (!is.null(weights) && is.matrix(weights)) weight.times <- attr(weights,"times")

    if (is.null(at.step)) {
        if (length(object$causes) == 1) {
            at.step <- list(object$stepno)
        } else {
            at.step <- lapply(object$model,function(arg) arg$stepno)[seq(from=1,by=length(object$strata),length.out=length(object$causes))]
        }
        names(at.step) <- object$causes
    } else {
        if (is.list(at.step)) {
            if (is.null(names(at.step))) names(at.step) <- object$causes
        } else {
            at.step.ori <- at.step
            at.step <- list()
            for (i in seq(along=object$causes)) at.step[[i]] <- at.step.ori
            names(at.step) <- object$causes
        }
    }

    if (!is.null(newdata)) {
        if (is.null(subset)) {
            subset.index <- 1:nrow(newdata)
        } else {
            subset.index <- (1:nrow(newdata))[subset]
            if (!is.null(newtime)) newtime <- newtime[subset]
            if (!is.null(newstatus)) newstatus <- newstatus[subset]
        }

        if (is.null(stratum)) {
            stratum <- rep(object$strata[1],length(subset.index))
        } else {
            stratum <- stratum[subset.index]
            if (any(!(stratum %in% object$strata))) {
                warning(paste("'stratum' contains strata no model has been fitted for. Replacing them by stratum '",object$strata[1],"'",sep=""))
                stratum[!(stratum %in% object$strata)] <- as.numeric(object$strata[1])
            }
        }

        if (is.null(weights)) {
            weights <- rep(1,length(subset.index))
        } else {
            if (is.matrix(weights)) {
                weights <- weights[subset.index,,drop=FALSE]
            } else {
                weights <- weights[subset.index]
            }
        }
    } else {
        stratum <- object$stratum
        weights <- object$weights

        if (is.matrix(weights)) {
            weight.times <- attr(object$weights,"times")
        }
    }

    res <- list()

    for (cause.index in seq(along=object$causes)) {
        actual.cause <- object$causes[cause.index]

        if (is.null(newdata)) {
            if (length(object$causes) > 1 || length(object$strata) > 1) {
                linear.predictor <- matrix(0,length(at.step[[actual.cause]]),length(object$time))
                for (stratum.index in seq(along=object$strata)) {
                    linear.predictor[object$stratum == object$strata[stratum.index]] <- object$model[[(cause.index-1)*length(object$strata)+stratum.index]]$linear.predictor[at.step[[actual.cause]]+1,,drop=FALSE]
                }
            } else {
                linear.predictor <- object$linear.predictor[at.step[[actual.cause]]+1,,drop=FALSE]
            }
        } else {
            linear.predictor <- matrix(0,length(at.step[[actual.cause]]),length(subset.index))

            for (stratum.index in seq(along=object$strata)) {
                model.index <- (cause.index-1)*length(object$strata)+stratum.index

                if (length(object$causes) > 1 || length(object$strata) > 1) {
                    nz.index <- which(Matrix::colSums(abs(object$model[[model.index]]$coefficients[at.step[[actual.cause]]+1,,drop=FALSE])) > 0)
                } else {
                    nz.index <- which(Matrix::colSums(abs(object$coefficients[at.step[[actual.cause]]+1,,drop=FALSE])) > 0)
                }

                if (length(nz.index) > 0) {
                    if (length(object$causes) > 1 || length(object$strata) > 1) {
                        nz.coef <- object$model[[model.index]]$coefficients[at.step[[actual.cause]]+1,nz.index,drop=FALSE]
                    } else {
                        nz.coef <- object$coefficients[at.step[[actual.cause]]+1,nz.index,drop=FALSE]
                    }

                    if (object$standardize) {
                        linear.predictor[,stratum == object$strata[stratum.index]] <- as.matrix(Matrix::tcrossprod(nz.coef,scale(newdata[subset.index[stratum == object$strata[stratum.index]],nz.index,drop=FALSE],center=object$meanx[nz.index],scale=object$sdx[nz.index])))
                    } else {
                        linear.predictor[,stratum == object$strata[stratum.index]] <- as.matrix(Matrix::tcrossprod(nz.coef,newdata[subset.index[stratum == object$strata[stratum.index]],nz.index,drop=FALSE]))
                    }
                }
            }
        }

        if (type == "lp") {
            if (length(object$causes) > 1) {
                res[[cause.index]] <- linear.predictor
                next
            } else {
                return(linear.predictor)
            }
        }

        if (type == "logplik") {
            if (is.null(newdata)) {
                newtime <- object$time
                newstatus <- object$status
            } else {
                if (is.null(newtime) || is.null(newstatus)) stop("'newtime' and 'newstatus' required for prediction on new data")
            }

            logplik <- numeric(length(at.step[[actual.cause]]))

            for (stratum.index in seq(along=object$strata)) {
                model.index <- (cause.index-1)*length(object$strata)+stratum.index
                actual.smask <- (stratum == object$strata[stratum.index])

                if (length(object$causes) > 1) {
                    actual.status <- ifelse(newstatus[actual.smask] == actual.cause,newstatus[actual.smask],0)
                } else {
                    actual.status <- newstatus[actual.smask]
                }

                if (is.matrix(weights)) {
                    if (is.null(newdata)) {
                        weightmat <- efron.weightmat(newtime[actual.smask],actual.status,actual.cause,object$model[[model.index]]$weight.at.times)
                    } else {
                        weightmat <- efron.weightmat(newtime[actual.smask],actual.status,actual.cause,weights[actual.smask,,drop=FALSE],prune.times=TRUE,weight.times=weight.times)
                    }
                } else {
                    weightmat <- efron.weightmat(newtime[actual.smask],actual.status,actual.cause,weights[actual.smask])
                }

                uncens <- which(newstatus == actual.cause & actual.smask)

                if (is.matrix(weights)) {
                    weights.at.event <- c()
                    for (i in seq(along=uncens)) {
                        weights.at.event <- c(weights.at.event,weights[uncens,,drop=FALSE][i,which(weight.times == newtime[uncens][i])[1]])
                    }
                } else {
                    weights.at.event <- weights[uncens]
                }

                for (i in seq(along=at.step[[actual.cause]])) {
                    logplik[i] <- logplik[i] + sum(weights.at.event*(linear.predictor[i,uncens] - log(apply(weightmat*exp(linear.predictor[i,actual.smask]),2,sum))))
                }
            }

            if (length(object$causes) > 1) {
                res[[cause.index]] <- logplik
                next
            } else {
                return(logplik)
            }
        }

        if (type == "risk" || type == "CIF") {
            if (max(unlist(lapply(at.step[[actual.cause]],length))) > 1) warning("predicted risk is only calculated for a single step (the first in 'at.step')")
            if (is.null(times)) times <- unique(object$time)

            if (length(object$causes) > 1) {
                res[[cause.index]] <- linear.predictor
            } else {
                pred.risk <- matrix(0,length(linear.predictor),length(times))

                for (stratum.index in seq(along=object$strata)) {
                    actual.smask <- (stratum == object$strata[stratum.index])

                    if (length(object$strata) > 1) {
                        breslow.Lambda <- unlist(lapply(times,function(x) ifelse(x < object$model[[stratum.index]]$event.times[1],0,object$model[[stratum.index]]$Lambda[at.step[[actual.cause]][1]+1,rev(which(object$model[[stratum.index]]$event.times <= x))[1]])))
                    } else {
                        breslow.Lambda <- unlist(lapply(times,function(x) ifelse(x < object$event.times[1],0,object$Lambda[at.step[[actual.cause]][1]+1,rev(which(object$event.times <= x))[1]])))
                    }
                    pred.risk[actual.smask,] <- exp(exp(linear.predictor[1,actual.smask]) %*% -t(breslow.Lambda))
                }

                if (type == "risk") {
                    return(pred.risk)
                } else {
                    return(1-pred.risk)
                }
            }
        }
    }

    if (type == "risk" || type == "CIF") {
        all.event.times <- object$event.times
        if (all.event.times[1] > 0) all.event.times <- c(0,all.event.times)
        if (max(all.event.times) < max(times)) {
            all.event.times <- c(all.event.times,max(times))
        } else {
            all.event.times <- all.event.times[all.event.times <= max(times)]
        }

        cum.haz <- matrix(0,length(res[[1]][1,]),length(all.event.times))

        for (cause.index in seq(along=object$causes)) {
            add.haz <- matrix(0,length(res[[1]][1,]),length(all.event.times))

            for (stratum.index in seq(along=object$strata)) {
                actual.smask <- (stratum == object$strata[stratum.index])
                model.index <- (cause.index-1)*length(object$strata)+stratum.index

                actual.bas.haz <- unlist(lapply(all.event.times,function(x) ifelse(x < object$model[[model.index]]$event.times[1],0,object$model[[model.index]]$Lambda[at.step[[cause.index]][1]+1,rev(which(object$model[[model.index]]$event.times <= x))[1]])))
                add.haz[actual.smask,] <- (exp(res[[cause.index]][1,actual.smask]) %*% t(actual.bas.haz))
            }

            cum.haz <- cum.haz + add.haz
        }

        all.surv <- exp(-cum.haz)
        t.all.surv <- t((all.surv[,1:(ncol(all.surv)-1)] + all.surv[,2:ncol(all.surv)])/2)
        #t.all.surv <- t(all.surv[,2:ncol(all.surv)])

        all.cif <- matrix(0,nrow(cum.haz),length(all.event.times))

        for (cause.index in seq(along=object$causes)) {
            actual.cif <- matrix(0,nrow(cum.haz),length(all.event.times))

            for (stratum.index in seq(along=object$strata)) {
                actual.smask <- (stratum == object$strata[stratum.index])
                model.index <- (cause.index-1)*length(object$strata)+stratum.index

                actual.event.times <- object$model[[model.index]]$event.times
                actual.cum.haz <- exp(res[[cause.index]][1,actual.smask]) %*% t(object$model[[model.index]]$Lambda[at.step[[cause.index]][1]+1,])

                if (actual.event.times[1] != 0) {
                    actual.cum.haz <- cbind(rep(0,nrow(actual.cum.haz)),actual.cum.haz)
                    actual.event.times <- c(0,actual.event.times)
                }
                actual.haz <- apply(actual.cum.haz,1,diff)/diff(actual.event.times)
                all.event.index <- unlist(lapply(all.event.times[-length(all.event.times)],function(x) rev(which(actual.event.times <= x))[1]))
                actual.haz <- rbind(actual.haz,rep(0,ncol(actual.haz)))[all.event.index,]

                actual.cif[actual.smask,] <- cbind(rep(0,length(res[[cause.index]][1,actual.smask])),t(apply(actual.haz * diff(all.event.times) * t.all.surv,2,cumsum)))
            }

            res[[cause.index]] <- actual.cif

            all.cif <- all.cif + actual.cif
        }

        times.index <- unlist(lapply(times,function(x) rev(which(all.event.times <= x))[1]))
        scale.cif <- (1 - all.surv)/all.cif
        scale.cif[,1] <- 0

        for (i in seq(along=object$causes)) {
            res[[i]] <- (res[[i]]*scale.cif)[,times.index]
            if (type == "risk") res[[i]] <- 1 - res[[i]]
        }
    }

    names(res) <- object$causes

    res
}



#' Determines the optimal number of boosting steps by cross-validation
#'
#' Performs a K-fold cross-validation for \code{\link{CoxBoost}} in search for
#' the optimal number of boosting steps.
#'
#'
#' @param time vector of length \code{n} specifying the observed times.
#' @param status censoring indicator, i.e., vector of length \code{n} with
#' entries \code{0} for censored observations and \code{1} for uncensored
#' observations. If this vector contains elements not equal to \code{0} or
#' \code{1}, these are taken to indicate events from a competing risk and a
#' model for the subdistribution hazard with respect to event \code{1} is
#' fitted (see e.g. Fine and Gray, 1999).
#' @param x \code{n * p} matrix of covariates.
#' @param subset a vector specifying a subset of observations to be used in the
#' fitting process.
#' @param maxstepno maximum number of boosting steps to evaluate, i.e, the
#' returned ``optimal'' number of boosting steps will be in the range
#' \code{[0,maxstepno]}.
#' @param K number of folds to be used for cross-validation. If \code{K} is
#' larger or equal to the number of non-zero elements in \code{status},
#' leave-one-out cross-validation is performed.
#' @param type way of calculating the partial likelihood contribution of the
#' observation in the hold-out folds: \code{"verweij"} uses the more
#' appropriate method described in Verweij and van Houwelingen (1996),
#' \code{"naive"} uses the approach where the observations that are not in the
#' hold-out folds are ignored (often found in other R packages).
#' @param parallel logical value indicating whether computations in the
#' cross-validation folds should be performed in parallel on a compute cluster,
#' using package \code{snowfall}. Parallelization is performed via the package
#' \code{snowfall} and the initialization function of of this package,
#' \code{sfInit}, should be called before calling \code{cv.CoxBoost}.
#' @param multicore indicates whether computations in the cross-validation
#' folds should be performed in parallel, using package \code{parallel}. If
#' \code{TRUE}, package \code{parallel} is employed using the default number of
#' cores. A value larger than \code{1} is taken to be the number of cores that
#' should be employed.
#' @param upload.x logical value indicating whether \code{x} should/has to be
#' uploaded to the compute cluster for parallel computation. Uploading this
#' only once (using \code{sfExport(x)} from library \code{snowfall}) can save
#' much time for large data sets.
#' @param folds if not \code{NULL}, this has to be a list of length \code{K},
#' each element being a vector of indices of fold elements. Useful for
#' employing the same folds for repeated runs.
#' @param trace logical value indicating whether progress in estimation should
#' be indicated by printing the number of the cross-validation fold and the
#' index of the covariate updated.
#' @param cmprsk type of competing risk, specific hazards or cause-specific
#' @param weights weights to be passed to \code{\link{predict}}
#' @param \dots miscellaneous parameters for the calls to
#' \code{\link{CoxBoost}}
#' @param stratum vector specifying different groups of individuals for a
#' stratified Cox regression. In \code{CoxBoost} fit each group gets its own
#' baseline hazard.
#' @return List with the following components: \item{mean.logplik}{vector of
#' length \code{maxstepno+1} with the mean partial log-likelihood for boosting
#' steps \code{0} to \code{maxstepno}} \item{se.logplik}{vector with standard
#' error estimates for the mean partial log-likelihood criterion for each
#' boosting step.} \item{optimal.step}{optimal boosting step number, i.e., with
#' minimum mean partial log-likelihood.} \item{folds}{list of length \code{K},
#' where the elements are vectors of the indices of observations in the
#' respective folds.}
#' @author Harald Binder \email{binderh@@uni-mainz.de}
#' @seealso \code{\link{CoxBoost}}, \code{\link{optimCoxBoostPenalty}}
#' @references Verweij, P. J. M. and van Houwelingen, H. C. (1993).
#' Cross-validation in survival analysis. Statistics in Medicine,
#' 12(24):2305-2314.
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
#'
#' #  10-fold cross-validation
#'
#' cv.res <- cv.CoxBoost(time=obs.time,status=status,x=x,maxstepno=500,
#'                       K=10,type="verweij",penalty=100)
#'
#' #   examine mean partial log-likelihood in the course of the boosting steps
#' plot(cv.res$mean.logplik)
#'
#' #   Fit with optimal number of boosting steps
#'
#' cbfit <- CoxBoost(time=obs.time,status=status,x=x,stepno=cv.res$optimal.step,
#'                   penalty=100)
#' summary(cbfit)
#'
#' }
#'
#' @export
cv.CoxBoost <- function(time,status,x,subset=1:length(time),weights=NULL,stratum=NULL,maxstepno=100,K=10,type=c("verweij","naive"),cmprsk=c("sh","csh","ccsh"),parallel=FALSE,upload.x=TRUE,
                        multicore=FALSE,folds=NULL,
                        trace=FALSE,...) {
    type <- match.arg(type)
    cmprsk <- match.arg(cmprsk)

    subset <- (1:length(time))[subset]

    if (is.null(stratum)) stratum <- rep(1,length(time))

    if (!is.null(folds) && length(folds) != K) stop("'folds' has to be of length 'K'")

    if (is.null(folds)) {
        if (K >= sum(status != 0)) {    #   leave-one-out cross-validation
            if (type == "verweij") {
                folds <- as.list(1:length(time))
            } else {
                folds <- as.list(which(status != 0))
            }
        } else {
            event.table <- table(status[subset][status[subset] != 0])

            while(TRUE) {
                folds <- split(sample(1:length(subset)), rep(1:K, length = length(subset)))

                #   make sure there is at least one event in training and test folds respectively
                #   Note: the Verweij approach actually could deal with folds that contain only
                #   censored observations, but this is expected to considerably increase variability
                #   and therefore alse prevented

                if (length(event.table) > 1 && cmprsk != "sh") {
                    if (!any(unlist(lapply(folds,function(fold) length(table(status[subset][fold][status[subset][fold] != 0])) != length(event.table)))) &&
                        !any(unlist(lapply(folds,function(fold) length(table(status[subset][-fold][status[subset][-fold] != 0])) != length(event.table)))))
                    {
                        break
                    }
                } else {
                    if (!any(unlist(lapply(folds,function(fold) sum(status[subset][fold] == 1))) == 0) &&
                        !any(unlist(lapply(folds,function(fold) sum(status[subset][-fold] == 1))) == 0))
                    {
                        break
                    }
                }
            }
        }
    }

    criterion <- NULL

    eval.fold <- function(actual.fold,...) {
        if (trace) cat("cv fold ",actual.fold,": ",sep="")
        cv.fit <- CoxBoost(time=time,status=status,x=x,subset=subset[-folds[[actual.fold]]],weights=weights,stratum=stratum,
                           stepno=maxstepno,return.score=FALSE,trace=trace,cmprsk=cmprsk,...)

        if (type == "verweij") {
            full.ploglik <- predict(cv.fit,newdata=x,newtime=time,newstatus=status,subset=subset,weights=weights,stratum=stratum,type="logplik",at.step=0:maxstepno)
            fold.ploglik <- predict(cv.fit,newdata=x,newtime=time,newstatus=status,subset=subset[-folds[[actual.fold]]],weights=weights,stratum=stratum,
                                    type="logplik",at.step=0:maxstepno)

            if (is.list(full.ploglik)) {
                res <- lapply(1:length(full.ploglik),function(arg) full.ploglik[[arg]] - fold.ploglik[[arg]])
                names(res) <- names(full.ploglik)
                return(res)
            } else {
                return(full.ploglik - fold.ploglik)
            }
        } else {
            return(predict(cv.fit,newdata=x,newtime=time,newstatus=status,subset=subset[folds[[actual.fold]]],
                                                 type="logplik",at.step=0:maxstepno))
        }
    }

    eval.success <- FALSE

    if (parallel) {
        if (!requireNamespace("snowfall")) {
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed using this package")
        } else {
            snowfall::sfLibrary(CoxBoost)
            if (upload.x) {
                snowfall::sfExport("time","status","x","maxstepno","trace","type","folds","cmprsk","weights","stratum")
            } else {
                snowfall::sfExport("time","status","maxstepno","trace","type","folds","cmprsk","weights","stratum")
            }
            criterion <- snowfall::sfClusterApplyLB(1:length(folds),eval.fold,...)
            eval.success <- TRUE
        }
    }

    if (!eval.success & multicore) {
        if (!requireNamespace("parallel")) {
            warning("package 'parallel' not found, i.e., parallelization cannot be performed using this package")
        } else {
            if (multicore > 1) {
                criterion <- parallel::mclapply(1:length(folds),eval.fold,mc.preschedule=FALSE,mc.cores=multicore,...)
            } else {
                criterion <- parallel::mclapply(1:length(folds),eval.fold,mc.preschedule=FALSE,...)
            }
            eval.success <- TRUE
        }
    }

    if (!eval.success) {
        criterion <- lapply(1:length(folds),eval.fold,...)
    }

    if (is.list(criterion[[1]])) {
        mean.criterion <- list()
        se.criterion <- list()
        optimal.step <- list()
        for (i in seq(along=criterion[[1]])) {
            actual.criterion <- matrix(unlist(lapply(criterion,function(arg) arg[[i]])),nrow=length(folds),byrow=TRUE)
            mean.criterion[[i]] <- apply(actual.criterion,2,mean)
            se.criterion[[i]] <- apply(actual.criterion,2,sd)/sqrt(nrow(actual.criterion))
            optimal.step[[i]] <- which.max(mean.criterion[[i]])-1
        }
        names(mean.criterion) <- names(se.criterion) <- names(optimal.step) <- names(criterion[[1]])
    } else {
        criterion <- matrix(unlist(criterion),nrow=length(folds),byrow=TRUE)
        mean.criterion <- apply(criterion,2,mean)
        se.criterion <- apply(criterion,2,sd)/sqrt(nrow(criterion))
        optimal.step <- which.max(mean.criterion)-1
    }

    list(mean.logplik=mean.criterion,se.logplik=se.criterion,optimal.step=optimal.step,
         folds=folds)
}




#' Coarse line search for adequate penalty parameter
#'
#' This routine helps in finding a penalty value that leads to an ``optimal''
#' number of boosting steps for CoxBoost, determined by cross-validation, that
#' is not too small/in a specified range.
#'
#' The penalty parameter for \code{\link{CoxBoost}} has to be chosen only very
#' coarsely.  In Tutz and Binder (2006) it is suggested for likelihood based
#' boosting just to make sure, that the optimal number of boosting steps,
#' according to some criterion such as cross-validation, is larger or equal to
#' 50.  With a smaller number of steps, boosting may become too ``greedy'' and
#' show sub-optimal performance.  This procedure uses a very coarse line search
#' and so one should specify a rather large range of boosting steps.
#'
#' @param time vector of length \code{n} specifying the observed times.
#' @param status censoring indicator, i.e., vector of length \code{n} with
#' entries \code{0} for censored observations and \code{1} for uncensored
#' observations. If this vector contains elements not equal to \code{0} or
#' \code{1}, these are taken to indicate events from a competing risk and a
#' model for the subdistribution hazard with respect to event \code{1} is
#' fitted (see e.g. Fine and Gray, 1999).
#' @param x \code{n * p} matrix of covariates.
#' @param minstepno,maxstepno range of boosting steps in which the ``optimal''
#' number of boosting steps is wanted to be.
#' @param start.penalty start value for the search for the appropriate penalty.
#' @param iter.max maximum number of search iterations.
#' @param upper.margin specifies the fraction of \code{maxstepno} which is used
#' as an upper margin in which a cross-validation minimum is not taken to be
#' one. This is necessary because of random fluctuations of cross-validated
#' partial log-likelihood.
#' @param parallel logical value indicating whether computations in the
#' cross-validation folds should be performed in parallel on a compute cluster.
#' Parallelization is performed via the package \code{snowfall} and the
#' initialization function of of this package, \code{sfInit}, should be called
#' before calling \code{cv.CoxBoost}.
#' @param trace logical value indicating whether information on progress should
#' be printed.
#' @param \dots miscellaneous parameters for \code{\link{cv.CoxBoost}}.
#' @return List with element \code{penalty} containing the ``optimal'' penalty
#' and \code{cv.res} containing the corresponding result of \code{cv.CoxBoost}.
#' @author Written by Harald Binder \email{binderh@@uni-mainz.de}.
#' @seealso \code{\link{CoxBoost}}, \code{\link{cv.CoxBoost}}
#' @references Tutz, G. and Binder, H. (2006) Generalized additive modelling
#' with implicit variable selection by likelihood based boosting.
#' \emph{Biometrics}, 62:961-971.
#' @keywords models smooth regression
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
#' #  determine penalty parameter
#'
#' optim.res <- optimCoxBoostPenalty(time=obs.time,status=status,x=x,
#'                                   trace=TRUE,start.penalty=500)
#'
#' #   Fit with obtained penalty parameter and optimal number of boosting
#' #   steps obtained by cross-validation
#'
#' cbfit <- CoxBoost(time=obs.time,status=status,x=x,
#'                   stepno=optim.res$cv.res$optimal.step,
#'                   penalty=optim.res$penalty)
#' summary(cbfit)
#'
#' }
#'
#' @export
optimCoxBoostPenalty <- function(time,status,x,minstepno=50,maxstepno=200,start.penalty=9*sum(status==1),
                                 iter.max=10,upper.margin=0.05,parallel=FALSE,trace=FALSE,...)
{
    if (parallel) {
        if (!requireNamespace("snowfall")) {
            parallel <- FALSE
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            snowfall::sfExport("x")
        }
    }

    actual.penalty <- start.penalty

    #   default: start from a large penalty and go down, when gone to far use small steps up
    step.up <- 1.2
    step.down <- 0.5

    actual.res <- NULL

    for (i in 1:iter.max) {
        if (trace) cat("iteration",i,": evaluating penalty",actual.penalty,"\n")

        actual.res <- cv.CoxBoost(time=time,status=status,x=x,maxstepno=maxstepno,penalty=actual.penalty,
                                  parallel=parallel,upload.x=FALSE,trace=trace,...)
        actual.max <- max(unlist(actual.res$optimal.step))

        if (trace) cat("maximum partial log-likelihood at boosting step",actual.max,"\n")
        if (actual.max >= minstepno && actual.max < maxstepno*(1-upper.margin)) break

        #   check whether we are in a scenario where penalty is far to low to start with
        if (i == 1 && actual.max < minstepno) {
            step.up <- 2
            step.down <- 0.8
        }

        if (actual.max < minstepno) {
            actual.penalty <- actual.penalty * step.up
        } else {
            actual.penalty <- actual.penalty * step.down
        }

        if (i == iter.max) warning("Exceeded iter.max in search for penalty parameter")
    }

    list(penalty=actual.penalty,cv.res=actual.res)
}

