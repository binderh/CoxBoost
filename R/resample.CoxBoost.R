#' Performs resampling for a stratified weighted Cox model by componentwise
#' likelihood based boosting for developing subgroup covariate signatures
#'
#' \code{resample.CoxBoost} is used to perform resampling for stratified Cox
#' proportional hazards models by componentwise likelihood based boosting with
#' weights for developing subgroup covariate signatures. Specifically, this
#' weighted Cox regression approach is for performing a subgroup analysis, not
#' for the standard subgroup analysis but for a weighted form of it where the
#' observations of the not analyzed subgroup are down-weighted.
#' \code{resample.CoxBoost} calls \code{cv.CoxBoost} and \code{CoxBoost}.
#'
#' The \code{CoxBoost} function can be used for performing variable selection
#' for time-to-event data based on componentwise likelihood based boosting
#' (Binder and Schumacher, 2009; Binder et al., 2013) for obtaining signatures
#' for high dimensional data, such as gene expression or high throughput
#' methylation measurements.  If there is any heterogeneity due to known
#' patient subgroups, it would be interesting for developing subgroup
#' signatures to account for this heterogeneity. To account for heterogeneity
#' in the data a stratified Cox model via \code{stratum} can be performed,
#' which allows for different baseline hazards in each subgroup. For performing
#' a stratified \code{CoxBoost} for such heterogeneous subgroups the
#' \code{stratum} can be set to a subgroup variable, containing also the
#' subgroup for which a subgroup signature should be developed.  A standard
#' subgroup analysis may lead to a decreased statistical power.  With a
#' weighted approach based on componentwise likelihood-based boosting one can
#' develop subgroup signatures without losing statistical power as it is the
#' case in standard subgroup analysis. This approach focuses on building a risk
#' prediction signature for a specific stratum by down-weighting the
#' observations (Simon, 2002) from the other strata using a range of weights.
#' The not analyzed subgroup has to be indicated in \code{stratnotinfocus}, the
#' observations of this subgroup, which is not in focus but should be retained
#' in the analysis, are weighted down by different relative weights. Binder et
#' al. (2012) propose such a weighting approach for high dimensional data but
#' not for time-to-event endpoints. \code{resample.CoxBoost} performs a
#' weighted regression approach for time-to-event data while using resampling
#' over \code{rep} resampling data sets, specifically for evaluating the
#' results. The different weights can be entered by \code{mix.list} which
#' should be from a range between zero and one.  For examining the variable
#' selection stability for specific covariates as a function of the weights the
#' beta coefficients are an output value of \code{resample.CoxBoost}, so that
#' resampling inclusion frequencies (Sauerbrei et. al, 2011) can be calculated
#' and visualized via the functions \code{stabtrajec} and \code{weightfreqmap}.
#' For mandatory variables (Binder and Schumacher, 2008), which are only
#' included for adjusting the model no penalty should be added, these can be
#' indicated by the indices in \code{unpen.index}.
#'
#' @param time vector of length \code{n} with the observed times of each
#' individual.
#' @param status vector of length \code{n} with entries \code{0} for censored
#' observations and \code{1} for observations having an event.
#' @param x \code{n * p} matrix of mandatory and optional covariates.
#' @param rep number of resampling data sets used for the evaluation of model
#' stability. Should be equal to \code{length(iter)}.
#' @param maxstepno maximum number of boosting steps needed for
#' cross-validation in \code{cv.CoxBoost}. So the optimal step number for
#' fitting the compontentwise boosting has a range between zero and
#' \code{maxstepno}.
#' @param multicore logical value: If \code{TRUE} then the cross-validation for
#' selecting the number of optimal boosting steps within each resampling data
#' set is been parallelized, if \code{FALSE} there will be no parallelization.
#' @param mix.list vector of different relative weights between zero and one,
#' for specifying weights for the individual observations for the weighted
#' regression approach. For each weight, a stratified Cox regression by
#' componentwise boosting is fitted.
#' @param stratum vector specifying different groups of individuals for a
#' stratified Cox regression. In \code{CoxBoost} fit each group gets its own
#' baseline hazard.
#' @param stratnotinfocus value for the group which is not analyzed but should
#' be weighted down in the analysis with the different weights element of
#' \code{mix.list}.
#' @param penalty penalty value for the covariates for the update of an
#' individual element in each boosting step. The mandatory covariates should
#' get a penalty value of zero, which can be realized with \code{unpen.index}.
#' @param criterion indicates the criterion to be used for selection in each
#' boosting step in \code{CoxBoost}. The default is \code{"hscore"} for using a
#' heuristic when evaluating a subset of covariates in each boosting step. For
#' the different other options see \code{CoxBoost}.
#' @param unpen.index index for mandatory covariates, which should get no
#' penalty value in the log-partial likelihood function.
#' @param trace logical value indicating whether progress in estimation should
#' be indicated by printing the number of the cross-validation fold and the
#' index of the covariate updated in the \code{cv.CoxBoost} and \code{CoxBoost}
#' methods.
#' @return \code{resample.CoxBoost} returns a list of length \code{rep} list
#' elements with each list element being further a list of the following two
#' objects:
#'
#' \item{CV.opt}{number of optimal boosting steps obtained from
#' \code{cv.CoxBoost} for each weight.} \item{beta}{beta coefficients obtained
#' from \code{CoxBoost} as vector for each weight.}
#' @author Written by Veronika Weyer \email{weyer@@uni-mainz.de} and Harald
#' Binder \email{binderh@@uni-mainz.de}.
#' @seealso \code{\link{stabtrajec}}, \code{\link{weightfreqmap}}.
#' @references Binder, H., Benner, A., Bullinger, L., Schumacher, M. (2013).
#' Tailoring sparse multivariable regression techniques for prognostic
#' single-nucleotide polymorphism signatures. Statistics in medicine 32(10),
#' 1778-1791
#'
#' Binder, H., Müller, T., Schwender, H., Golka, K., Steffens, M., G., H.J.,
#' Ickstadt, K., and Schumacher, M (2012). Cluster-Localized Sparse Logistic
#' Regression for SNP Data. Statistical applications in genetics and molecular
#' biology, 11(4), 1-31
#'
#' Binder, H. and Schumacher, M. (2009). Incorporating pathway information into
#' boosting estimation of high-dimensional risk prediction models. BMC
#' Bioinformatics. 10:18.
#'
#' Binder, H. and Schumacher, M. (2008). Allowing for mandatory covariates in
#' boosting estimation of sparse high-dimensional survival models. BMC
#' Bioinformatics. 9:14.
#'
#' Sauerbrei, W., Boulesteix, A.-L., Binder, H. (2011). Stability
#' investigations of multivariable regression models derived from low- and
#' high-dimensional data. Journal of biopharmaceutical statistics 21(6),
#' 1206-31
#'
#' Simon, R. (2002). Bayesian subset analysis: application to studying
#' treatment-by-gender interactions. Statistics in medicine 21(19), 2909-16
#' @keywords subgroup signature time-to-event endpoint weighted regression
#' stratified Cox model
#' @examples
#'
#' \donttest{
#' # Generate survival data with five informative covariates for one subgroup
#' n <- 400; p <- 1000
#' set.seed(129)
#' group<-rbinom(n,1,0.5)
#' x <- matrix(rnorm(n*p,0,1),n,p)
#' beta.vec1  <- c(c(1,1,1,1,1),rep(0,p-5))
#' beta.vec0  <- c(c(0,0,0,0,0),rep(0,p-5))
#' linpred<-ifelse(group==1,x %*% beta.vec1,x %*% beta.vec0)
#' set.seed(1234)
#' real.time<- (-(log(runif(n)))/(1/20*exp(linpred)))
#' cens.time <- rexp(n,rate=1/20)
#' obs.status <- ifelse(real.time <= cens.time,1,0)
#' obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)
#'
#'
#' # Fit a stratified weighted Cox proportional hazards model
#' # by \code{resample.CoxBoost}
#'
#' mix.list=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99)
#' RIF <- resample.CoxBoost(
#'   time=obs.time,status=obs.status,x=x,
#'   # use more repetitions (eg `rep = 100`) for more stable results
#'   rep=5,
#'   maxstepno=200,multicore=FALSE,
#'   mix.list=mix.list,
#'   stratum=group,stratnotinfocus=0,penalty=sum(obs.status)*(1/0.02-1),
#'   criterion="hscore",unpen.index=NULL)
#'
#' #   RIF is a list with number of resampling data sets list objects, each list
#' # object has further two list objects, the beta coefficients for all
#' # weights and the selected number of boosting steps.
#' #   For each list object i \code{RIF[[i]]$beta} contains the beta
#' # coefficients with length mix.list*p for p covariates and with
#' # \code{RIF[[i]]$CV.opt}
#' #   For getting an insight in the RIF Distribution of different weights
#' # use the following code:
#'
#' RIF1<-c()
#' for (i in 1: length(RIF)){RIF1<-c(RIF1,RIF[[i]][[1]])}
#' freqmat <-matrix(apply(matrix(unlist(RIF1), ncol=length(RIF))!=0,1,mean),
#' ncol=length(mix.list))
#'
#' #  freqmat is a matrix with p rows and \code{length(mix.list)} columns
#' #  which contains resampling inclusion frequencies for the different
#' #  covariates and different weights.
#'
#' # two plotting functions are available for the resulting object:
#' stabtrajec(RIF)
#' weightfreqmap(RIF)
#' }
#' @export resample.CoxBoost
resample.CoxBoost <- function(time,status,x,rep=100,maxstepno=200,multicore=TRUE,
                       mix.list=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99),
                       stratum,stratnotinfocus=0,
                       penalty=sum(status)*(1/0.02-1),criterion="hscore",unpen.index=NULL,
                       trace=FALSE)
{
  rep <- rep
  trainind <- list()
  for (i in 1:rep){
    trainind[[length(trainind)+1]]  <- sample(1:nrow(x),round(nrow(x)*0.632),replace = F)
  }

  out <- list()
  for (iter in 1:rep) {
    message('iter=', iter)
    outbeta<-c()
    outCV.opt<-c()
    for (mix.prop in mix.list) {
      message('weight=', mix.prop)
      obs.weights <- rep(1,length(status))
      case.weights <- ifelse(stratum == stratnotinfocus,mix.prop,1)
      obs.weights <- case.weights/sum(case.weights)*length(case.weights)

      CV <- cv.CoxBoost(time=time[trainind[[iter]]],status=status[trainind[[iter]]],x=x[trainind[[iter]],],
                        stratum=stratum[trainind[[iter]]],unpen.index=unpen.index,
                        coupled.strata = FALSE,weights=obs.weights[trainind[[iter]]],
                        maxstepno=maxstepno,K=10,penalty=penalty,
                        standardize=TRUE,trace=trace, multicore=multicore,criterion=criterion)

      CB <- CoxBoost(time=time[trainind[[iter]]],status=status[trainind[[iter]]],x=x[trainind[[iter]],],
                    stratum=stratum[trainind[[iter]]],unpen.index=unpen.index,
                    coupled.strata = FALSE,weights=obs.weights[trainind[[iter]]],
                    stepsize.factor=1,stepno=CV$optimal.step,penalty=penalty,
                    standardize=TRUE,trace=trace,criterion=criterion)
      outbeta <- c(outbeta,CB$model[[1]][[5]][nrow(CB$model[[1]][[5]]),])
      outCV.opt <- c(outCV.opt,CV$optimal.step)
    }
    out[[iter]] <- list(beta=outbeta,CV.opt=outCV.opt)
  }

  out
}




#' Plots stability trajectories from resCoxBoost fits
#'
#' Plots the stability trajectories, i.e. the variable selection stability is
#' plotted in form of resampling inclusion frequencies for different weights
#' and a selection of different stably selected genes as obtained from a
#' resCoxBoost object fitted by \code{\link{resample.CoxBoost}}.
#'
#' The \code{stabtrajec} function produces a visualization tool (stability
#' trajectory plot), which shows the RIFs as a function of weights for stably
#' selected covariates selected via \code{resample.CoxBoost}. With this tool it
#' can be shown that with weights between zero and one and so the additional
#' information from the observations of the not analyzed subgroup in comparison
#' to the standard subgroup analysis (weight near zero) can improve the
#' variable selection stability. The number of plotted covariates in this plot
#' depends on \code{lowerRIFlimit}. How many weights are plotted which are
#' element of \code{mix.list} can be changed with \code{plotmix}.
#'
#' @param RIF list obtained from a \code{\link{resample.CoxBoost}} call.
#' @param mix.list vector of weights which were also entered in the
#' \code{\link{resample.CoxBoost}} call.
#' @param plotmix vector of weights which should be plotted in the stability
#' trajectory plot.
#' @param my.colors vector with length(plotmix) of different colors for
#' plotting the different weights, default are gray shades.
#' @param yupperlim value for the upper y coordinate.
#' @param huge size of the labels plottet on the x-axis.
#' @param lowerRIFlimit cutoff RIF value: All covariates which have an
#' resampling inclusion frequency (RIF) greater or equal this lowerRIFlimit
#' value at any weight are presented in the stability trajectory plot.
#' @param legendval space between the last plotted covariate and the legend on
#' the right side of the plot.
#' @return No value is returned, but the stability trajectory plot is
#' generated.
#' @author Veronika Weyer \email{weyer@@uni-mainz.de} and Harald Binder
#' \email{binderh@@uni-mainz.de}
#' @export
stabtrajec<-function(RIF,mix.list=c(0.001,0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99)
                     ,plotmix=c(0.001,0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99)
                     ,my.colors=grDevices::gray(seq(.99,0,len=10)),
                     yupperlim=1,huge=0.6,lowerRIFlimit=0.6,legendval=4.5)
{
  RIF1<-c()
  for (i in 1: length(RIF)){RIF1<-c(RIF1,RIF[[i]][[1]])}
  freqmat <-matrix(apply(matrix(unlist(RIF1), ncol=length(RIF))!=0,1,mean), ncol=length(mix.list))
  sel.mask <- apply(freqmat,1,function(arg) any(arg  >= lowerRIFlimit & arg < 1.1))
  w5<-c(1:length(which(sel.mask==T)))
  jitmat<-cbind(w5-0.28,w5-0.21,w5-0.14,w5-0.07,w5,w5+0.07,w5+0.14,w5+0.21,w5+0.28,w5+0.35)
  colnames(jitmat)<-mix.list
  colnames(freqmat)<-mix.list
  freqmat<-freqmat[,which(colnames(freqmat)%in%paste("",plotmix,sep=""))]
  jitmat<-jitmat[,which(colnames(jitmat)%in%paste("",plotmix,sep=""))]
  plot(0,xlim=c(0.5,length(which(sel.mask==T))+legendval),ylim=c(0,yupperlim),type="n",main=" ",
       xlab=" ",ylab="resampling inclusion frequency", las=T,xaxt = "n")
  graphics::axis(1, at = c(1:length(which(sel.mask==T))),labels =rownames(freqmat[sel.mask,]),cex.axis=huge)
  for (i in 1:length(plotmix)){ points(jitmat[,i],freqmat[,i][sel.mask],col=my.colors[i],type = 'p',pch=16) }
  for (i in 1:length(plotmix)){ points(jitmat[,i],freqmat[,i][sel.mask],col=1,type = 'p') }
  for (i in c(1:length(which(sel.mask==T)))){ graphics::lines(jitmat[i,],freqmat[sel.mask,][i,],col=1) }
  for (i in 1:length(plotmix)) { graphics::legend("topright",paste("w=",plotmix, sep=""),pch=16,col=my.colors,bty="n") }
  for (i in 1:length(plotmix)) { graphics::legend("topright",paste("w=",plotmix, sep=""),pch=1,col=1,bty="n") }
  graphics::abline(v=c(1:length(which(sel.mask==T))), col=grDevices::gray(0.7))
}


#' @export
#' @rdname stabtrajec
#' @param method passed to \code{\link{hclust}} for the the agglomeration method to be used
weightfreqmap<-function(RIF,mix.list=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99)
                        ,plotmix=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99)
                        ,lowerRIFlimit=0.5,method="complete")
{
  RIF1<-c()
  for (i in 1: length(RIF)){RIF1<-c(RIF1,RIF[[i]][[1]])}
  freqmat <-matrix(apply(matrix(unlist(RIF1), ncol=length(RIF))!=0,1,mean), ncol=length(mix.list))
  colnames(freqmat)<-mix.list
  sel.indz<-apply(freqmat,1,function(arg) any(arg  >= lowerRIFlimit & arg < 1.1))
  heatcol<-grDevices::gray(seq(0,.9,len=100))
  heatmap(freqmat[sel.indz,which(colnames(freqmat)%in%paste("",plotmix,sep=""))],
          col=heatcol,hclustfun = function(x) hclust(x,method=method),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          xlab ="relative weights", scale="row",Colv=NA)
}
