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
            DUP=FALSE,NAOK=TRUE
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
                    DUP=FALSE,NAOK=TRUE
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
                  DUP=FALSE
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

            model[[model.index]]$coefficients <- Matrix(0,stepno+1,p)
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
                        if (class(try.res) == "try-error") {
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
            combined.coefficients <- Matrix(0,nrow(model[[i]]$coefficients),object$p)
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

print.CoxBoost <- function(x,...) {
    joint.print(x,long=FALSE)
}

summary.CoxBoost <- function(object,...) {
    joint.print(object,long=TRUE)
}

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

            if (length(nz.index[[i]]) < ncoef - length(x$unpen.index)) lines(c(0,nrow(plotmat[[i]])-1),c(0,0),col=line.col)

            for (coef.index in 1:ncol(plotmat[[i]])) {
                lines(0:(nrow(plotmat[[i]])-1),plotmat[[i]][,coef.index],col=line.col)
                text(actual.xlim[2],plotmat[[i]][nrow(plotmat[[i]]),coef.index],plot.names[[i]][coef.index],pos=2,cex=label.cex)
            }
        }    
    }
}

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
        if (!require(snowfall)) {
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
        if (!require(parallel)) {
            warning("package 'parallel' not found, i.e., parallelization cannot be performed using this package")
        } else {
            if (multicore > 1) {
                criterion <- mclapply(1:length(folds),eval.fold,mc.preschedule=FALSE,mc.cores=multicore,...)
            } else {
                criterion <- mclapply(1:length(folds),eval.fold,mc.preschedule=FALSE,...)                
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


optimCoxBoostPenalty <- function(time,status,x,minstepno=50,maxstepno=200,start.penalty=9*sum(status==1),
                                 iter.max=10,upper.margin=0.05,parallel=FALSE,trace=FALSE,...)
{
    if (parallel) {
        if (!require(snowfall)) {
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

