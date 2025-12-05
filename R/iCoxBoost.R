#' Control parameters for cross-validation in \code{iCoxBoost}
#'
#' This function allows to set the control parameters for cross-validation to
#' be passed into a call to \code{\link{iCoxBoost}}.
#'
#'
#' @param K number of folds to be used for cross-validation. If \code{K} is
#' larger or equal to the number of events in the data to be analyzed,
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
#' \code{sfInit}, should be called before calling \code{iCoxBoost}.
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
#' @return List with elements corresponding to the call arguments.
#' @author Written by Harald Binder \email{binderh@@uni-mainz.de}.
#' @seealso \code{\link{iCoxBoost}}, \code{\link{cv.CoxBoost}}
#' @references Verweij, P. J. M. and van Houwelingen, H. C. (1993).
#' Cross-validation in survival analysis. Statistics in Medicine,
#' 12(24):2305-2314.
#' @keywords models smooth regression
#' @export
cvcb.control <- function(K=10,type=c("verweij","naive"),parallel=FALSE,upload.x=TRUE,multicore=FALSE,
                 folds=NULL)
{
	type <- match.arg(type)
	list(K=K,type=type,parallel=parallel,upload.x=upload.x,multicore=multicore,folds=folds)
}



#' Interface for cross-validation and model fitting using a formula description
#'
#' Formula interface for fitting a Cox proportional hazards model by
#' componentwise likelihood based boosting (via a call to
#' \code{\link{CoxBoost}}), where cross-validation can be performed
#' automatically for determining the number of boosting steps (via a call to
#' \code{\link{cv.CoxBoost}}).
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
#' coefficients being zero. The main model complexity parameter, the number of
#' boosting steps, is automatically selected by cross-validation using a call
#' to \code{\link{cv.CoxBoost}}). Note that this will introduce random
#' variation when repeatedly calling \code{iCoxBoost}, i.e. it is advised to
#' set/save the random number generator state for reproducible results.
#'
#' The advantage of the offset-based approach compared to gradient boosting is
#' that the penalty structure is very flexible. In the present implementation
#' this is used for allowing for unpenalized mandatory covariates, which
#' receive a very fast coefficient build-up in the course of the boosting
#' steps, while the other (optional) covariates are subjected to penalization.
#' For example in a microarray setting, the (many) microarray features would be
#' taken to be optional covariates, and the (few) potential clinical covariates
#' would be taken to be mandatory, by including their names in
#' \code{mandatory}.
#'
#' If a group of correlated covariates has influence on the response, e.g.
#' genes from the same pathway, componentwise boosting will often result in a
#' non-zero estimate for only one member of this group. To avoid this,
#' information on the connection between covariates can be provided in
#' \code{varlink}. If then, in addition, a penalty updating scheme with
#' \code{stepsize.factor} < 1 is chosen, connected covariates are more likely
#' to be chosen in future boosting steps, if a directly connected covariate has
#' been chosen in an earlier boosting step (see Binder and Schumacher, 2009b).
#'
#' @param formula A formula describing the model to be fitted, similar to a
#' call to \code{coxph}. The response must be a survival object, either as
#' returned by \code{Surv} or \code{Hist} (in a competing risks application).
#' @param data data frame containing the variables described in the formula.
#' @param weights optional vector, for specifying weights for the individual
#' observations.
#' @param subset a vector specifying a subset of observations to be used in the
#' fitting process.
#' @param mandatory vector containing the names of the covariates whose effect
#' is to be estimated un-regularized.
#' @param cause cause of interest in a competing risks setting, when the
#' response is specified by \code{Hist} (see e.g. Fine and Gray, 1999; Binder
#' et al. 2009a).
#' @param standardize logical value indicating whether covariates should be
#' standardized for estimation. This does not apply for mandatory covariates,
#' i.e., these are not standardized.
#' @param stepno maximum number of boosting steps to be evaluated when
#' determining the number of boosting steps by cross-validation, otherwise the
#' number of boosting seps itself.
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
#' @param nu (roughly) the fraction of the partial maximum likelihood estimate
#' used for the update in each boosting step. This is converted into a penalty
#' for the call to \code{CoxBoost}. Use smaller values, e.g., 0.01 when there
#' is little information in the data, and larger values, such as 0.1, with much
#' information or when the number of events is larger than the number of
#' covariates. Note that the default for direct calls to \code{CoxBoost}
#' corresponds to \code{nu=0.1}.
#' @param stepsize.factor determines the step-size modification factor by which
#' the natural step size of boosting steps should be changed after a covariate
#' has been selected in a boosting step. The default (value \code{1}) implies
#' constant \code{nu}, for a value < 1 the value \code{nu} for a covariate is
#' decreased after it has been selected in a boosting step, and for a value > 1
#' the value \code{nu} is increased. If \code{pendistmat} is given, updates of
#' \code{nu} are only performed for covariates that have at least one
#' connection to another covariate.
#' @param varlink list for specifying links between covariates, used to
#' re-distribute step sizes when \code{stepsize.factor != 1}. The list needs to
#' contain at least two vectors, the first containing the name of the source
#' covariates, the second containing the names of the corresponding target
#' covariates, and a third (optional) vector containing weights between 0 and 1
#' (defaulting to 1). If \code{nu} is increased/descreased for one of the
#' source covariates according to \code{stepsize.factor}, the \code{nu} for the
#' corresponding target covariate is descreased/increased accordingly
#' (multiplied by the weight). If \code{formula} contains interaction terms,
#' als rules for these can be set up, using variable names such as \code{V1:V2}
#' for the interaction term between covariates \code{V1} and \code{V2}.
#' @param cv \code{TRUE}, for performing cross-validation, with default
#' parameters, \code{FALSE} for not performing cross-validation, or list
#' containing the parameters for cross-validation, as obtained from a call to
#' \code{\link{cvcb.control}}.
#' @param trace logical value indicating whether progress in estimation should
#' be indicated by printing the name of the covariate updated.
#' @param ... miscellaneous arguments, passed to the call to
#' \code{\link{cv.CoxBoost}}.
#' @return \code{iCoxBoost} returns an object of class \code{iCoxBoost}, which
#' also has class \code{CoxBoost}. In addition to the elements from
#' \code{\link{CoxBoost}} it has the following elements: \item{call, formula,
#' terms}{call, formula and terms from the formula interface.}
#' \item{cause}{cause of interest.} \item{cv.res}{result from
#' \code{\link{cv.CoxBoost}}, if cross-validation has been performed.}
#' @author Written by Harald Binder \email{binderh@@uni-mainz.de}.
#' @seealso \code{\link{predict.iCoxBoost}}, \code{\link{CoxBoost}},
#' \code{\link{cv.CoxBoost}}.
#' @param cmprsk type of competing risk, specific hazards or cause-specific
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
#' beta <- c(rep(1,2),rep(0,p-2))
#' x <- matrix(rnorm(n*p),n,p)
#' actual.data <- as.data.frame(x)
#' real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
#' cens.time <- rexp(n,rate=1/10)
#' actual.data$status <- ifelse(real.time <= cens.time,1,0)
#' actual.data$time <- ifelse(real.time <= cens.time,real.time,cens.time)
#'
#' #   Fit a Cox proportional hazards model by iCoxBoost
#'
#' \donttest{cbfit <- iCoxBoost(Surv(time,status) ~ .,data=actual.data)
#' summary(cbfit)
#' plot(cbfit)}
#'
#' #   ... with covariates 1 and 2 being mandatory
#'
#' \donttest{cbfit.mand <- iCoxBoost(Surv(time,status) ~ .,data=actual.data,mandatory=c("V1"))
#' summary(cbfit.mand)
#' plot(cbfit.mand)}
#'
#'
#' @export
iCoxBoost <- function(formula,data=NULL,weights=NULL,subset=NULL,
					  mandatory=NULL,cause=1,cmprsk=c("sh","csh","ccsh"),standardize=TRUE,
					  stepno=200,criterion=c("pscore","score","hpscore","hscore"),
					  nu=0.1,stepsize.factor=1,varlink=NULL,
					  cv=cvcb.control(),trace=FALSE,...)
{
	criterion <- match.arg(criterion)
	cmprsk <- match.arg(cmprsk)

    call <- match.call(expand.dots=FALSE)
    formula <- eval(call$formula)

    if (is.null(data)) {
    	actual.terms <- terms(formula)
    	response <- model.response(model.frame(formula,data=environment(formula)))
	    x <- model.matrix(actual.terms)[,-c(1),drop=FALSE]
    } else {
    	actual.terms <- terms(formula,data=data)
    	response <- model.response(model.frame(formula,data))
	    x <- model.matrix(actual.terms,data=data)[,-c(1),drop=FALSE]
    }

    term.labels <- attr(actual.terms,"term.labels")
    stratum <- NULL
    if (any(substr(term.labels,1,7) == "strata(")) {
    	stratum.var <- term.labels[substr(term.labels,1,7) == "strata("][1]
    	stratum.var <- substr(stratum.var,8,nchar(stratum.var)-1)
    	if (stratum.var %in% colnames(x)) {
    	    stratum <- x[,stratum.var]
    	} else {
    	    if (is.null(data)) {
    	        get(stratum.var,envir=environment(formula))
    	    } else {
    	        stratum <- data[,stratum.var]
    	    }
    	}
    	x <- x[,substr(colnames(x),1,7) != "strata("]
    	x <- x[,colnames(x) != stratum.var]
    }

    actual.time <- as.numeric(response[,"time"])
    if (requireNamespace("prodlim", quietly = TRUE) && inherits(response, "Hist")) {
        if (cmprsk == "sh") {
            actual.status <- rep(2,NROW(response))
            actual.status[as.numeric(response[,"event"] == cause) == 1] <- 1
            actual.status[as.numeric(response[,"status"]) == 0] <- 0
        } else {
            actual.status <- response[,"event"]
            actual.status[as.numeric(response[,"status"]) == 0] <- 0
            actual.status[actual.status != 0] <- as.numeric(attr(response,"states"))[actual.status[actual.status != 0]]
        }
    } else {
    	actual.status <- as.numeric(response[,"status"])
    }

    if (is.null(subset)) subset <- 1:length(actual.time)
    if (cmprsk == "sh") {
        penalty <- sum(actual.status[subset]==cause)*(1/nu-1)
    } else {
        event.no <- table(actual.status[actual.status != 0])
        if (length(event.no) > 1) {
            penalty <- as.list(event.no*(1/nu-1))
        } else {
            penalty <- sum(actual.status[subset]==cause)*(1/nu-1)
        }
    }

    actual.names <- colnames(x)
    unpen.index <- NULL
    if (!is.null(mandatory)) {
    	if (sum(mandatory %in% actual.names) < length(mandatory)) {
    		stop("Some variables in 'mandatory' are not part of the model.")
    	}
    	unpen.index <- which(actual.names %in% mandatory)
    }

    pendistmat <- NULL
    connected.index <- NULL
    if (!is.null(varlink)) {
    	if (!is.list(varlink) || length(varlink) < 2) {
    		stop("'varlink' has to be a list with at least two elements")
    	}

    	actual.source <- varlink[[1]]
    	actual.target <- varlink[[2]]

    	if (length(varlink) == 2) {
    		actual.factor <- rep(1,length(actual.source))
    	} else {
			actual.factor <- varlink[[3]]
    	}

    	if (length(actual.source) != length(actual.target) || length(actual.source) != length(actual.factor)) {
    		stop("Source and target (and factor) vectors of 'varlink' have to be of the same length.")
    	}

    	actual.connected <- union(actual.source,actual.target)
    	if (!all(actual.connected %in% actual.names)) {
    		stop("Some elements of 'varlink' do not appear to be part of the model.")
    	}

    	connected.index <- which(actual.names %in% actual.connected)
    	actual.source.match <- match(actual.source,actual.names[connected.index])
    	actual.target.match <- match(actual.target,actual.names[connected.index])

    	pendistmat <- Matrix(0,length(connected.index),length(connected.index))

    	for (i in seq(along=actual.source.match)) {
    		pendistmat[actual.source.match[i],actual.target.match[i]] <- actual.factor[i]
    	}
    }

    use.stepno <- stepno
    cv.res <- NULL
    if (is.list(cv) || cv == TRUE) {
    	if (!is.list(cv)) cv <- cvcb.control()
	    cv.res <- cv.CoxBoost(time=actual.time,status=actual.status,x=x,
	    					  maxstepno=stepno,K=cv$K,type=cv$type,parallel=cv$parallel,
	    					  upload.x=cv$upload.x,multicore=cv$multicore,folds=cv$folds,
	    					  weights=weights,subset=subset,stratum=stratum,standardize=standardize,
	    					  penalty=penalty,criterion=criterion,cmprsk=cmprsk,
	    					  stepsize.factor=stepsize.factor,
	    					  unpen.index=unpen.index,pendistmat=pendistmat,
	    					  connected.index=connected.index,
	    					  trace=trace,...)
	    use.stepno <- cv.res$optimal.step
    }

    res <- CoxBoost(time=actual.time,status=actual.status,x=x,
    				weights=weights,subset=subset,stratum=stratum,standardize=standardize,
    				penalty=penalty,criterion=criterion,cmprsk=cmprsk,
    				stepno=use.stepno,unpen.index=unpen.index,pendistmat=pendistmat,
	    			stepsize.factor=stepsize.factor,
    				connected.index=connected.index,
    				trace=trace,...)

    res$call <- call
    res$formula <- formula
    res$terms <- actual.terms
    res$cause.code <- cause
    res$cv.res <- cv.res

    class(res) <- c("iCoxBoost",class(res))
    res
}



#' Predict method for iCoxBoost fits
#'
#' Obtains predictions at specified boosting steps from a iCoxBoost object
#' fitted by \code{\link{iCoxBoost}}.
#'
#'
#' @param object fitted CoxBoost object from a \code{\link{CoxBoost}} call.
#' @param newdata data frame with new covariate values (for \code{n.new}
#' observations). If just prediction for the training data is wanted, it can be
#' omitted. If the predictive partial log-likelihood is wanted
#' (\code{type=logplik}), this frame also has to contain the response
#' information.
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
#' n <- 200; p <- 100
#' beta <- c(rep(1,2),rep(0,p-2))
#' x <- matrix(rnorm(n*p),n,p)
#' actual.data <- as.data.frame(x)
#' real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
#' cens.time <- rexp(n,rate=1/10)
#' actual.data$status <- ifelse(real.time <= cens.time,1,0)
#' actual.data$time <- ifelse(real.time <= cens.time,real.time,cens.time)
#'
#' #   define training and test set
#'
#' train.index <- 1:100
#' test.index <- 101:200
#'
#' #   Fit a Cox proportional hazards model by iCoxBoost
#'
#' \donttest{cbfit <- iCoxBoost(Surv(time,status) ~ .,data=actual.data[train.index,],
#' 				   stepno=300,cv=FALSE)}
#'
#' #   mean partial log-likelihood for test set in every boosting step
#'
#' \donttest{step.logplik <- predict(cbfit,newdata=actual.data[test.index,],
#'                         at.step=0:300,type="logplik")
#'
#' plot(step.logplik)}
#'
#' #   names of covariates with non-zero coefficients at boosting step
#' #   with maximal test set partial log-likelihood
#'
#' \donttest{print(coef(cbfit,at.step=which.max(step.logplik)-1))}
#'
#'
#' @export
predict.iCoxBoost <- function(object,newdata=NULL,subset=NULL,at.step=NULL,times=NULL,
							  type=c("lp","logplik","risk","CIF"),...)
{
    type <- match.arg(type)

	newtime <- NULL
	newstatus <- NULL

	if (is.null(newdata)) {
		newdata <- NULL
	} else {
		if (type == "logplik") {
			new.response <- model.response(model.frame(object$formula,newdata))
			newtime <- as.numeric(new.response[,"time"])

			if (inherits(new.response, "Hist")) {
			    if (length(object$causes) > 1) {
                    newstatus <- new.response[,"event"]
                    newstatus[as.numeric(new.response[,"status"]) == 0] <- 0
                    newstatus[newstatus != 0] <- as.numeric(attr(new.response,"states"))[newstatus[newstatus != 0]]
			    } else {
                    newstatus <- rep(2,NROW(new.response))
                    newstatus[as.numeric(new.response[,"event"] == object$cause.code) == 1] <- 1
                    newstatus[as.numeric(new.response[,"status"]) == 0] <- 0
                }

			} else {
				newstatus <- as.numeric(new.response[,"status"])
			}
		}

        term.labels <- attr(object$terms,"term.labels")
        stratum <- NULL
        stratum.var <- NULL
        if (any(substr(term.labels,1,7) == "strata(")) {
    	    stratum.var <- term.labels[substr(term.labels,1,7) == "strata("][1]
    	    stratum.var <- substr(stratum.var,8,nchar(stratum.var)-1)
    	    stratum <- newdata[,stratum.var]
        }

		newdata <- model.matrix(object$terms,
				   data=model.frame(object$formula,data=newdata))[,-c(1),drop=FALSE]

        if (!is.null(stratum.var)) {
        	newdata <- newdata[,substr(colnames(newdata),1,7) != "strata("]
        	newdata <- newdata[,colnames(newdata) != stratum.var]
        }
	}

	predict.CoxBoost(object,newdata=newdata,newtime=newtime,newstatus=newstatus,subset=subset,at.step=at.step,times=times,type=type,stratum=stratum,...)
}

