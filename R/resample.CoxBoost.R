resample.CoxBoost<- function(time,status,x,rep=100,maxstepno=200,multicore=TRUE,
                       mix.list=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99),
                       stratum,stratnotinfocus=0,
                       penalty=sum(status)*(1/0.02-1),criterion="hscore",unpen.index=NULL)
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
    print(mix.prop)
    obs.weights <- rep(1,length(status))
    case.weights <- ifelse(stratum == stratnotinfocus,mix.prop,1)
    obs.weights <- case.weights/sum(case.weights)*length(case.weights)
    set.seed(x[1,5]*100+time[19]*10)
    CV <- cv.CoxBoost(time=time[trainind[[iter]]],status=status[trainind[[iter]]],x=x[trainind[[iter]],],
                      stratum=stratum[trainind[[iter]]],unpen.index=unpen.index,
                      coupled.strata = FALSE,weights=obs.weights[trainind[[iter]]],
                      maxstepno=maxstepno,K=10,penalty=penalty,
                      standardize=TRUE,trace=TRUE, multicore=multicore,criterion=criterion)
    set.seed(x[1,5]*100+time[19]*10)
    CB <- CoxBoost(time=time[trainind[[iter]]],status=status[trainind[[iter]]],x=x[trainind[[iter]],],
                   stratum=stratum[trainind[[iter]]],unpen.index=unpen.index,
                   coupled.strata = FALSE,weights=obs.weights[trainind[[iter]]],
                   stepsize.factor=1,stepno=CV$optimal.step,penalty=penalty,
                   standardize=TRUE,trace=TRUE,criterion=criterion)
    outbeta<-c(outbeta,CB$model[[1]][[5]][nrow(CB$model[[1]][[5]]),] )
    outCV.opt <- c(outCV.opt,CV$optimal.step)
  }
  out[[iter]] <- list(beta=outbeta,CV.opt=outCV.opt)
  }

  out
}


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
  axis(1, at = c(1:length(which(sel.mask==T))),labels =rownames(freqmat[sel.mask,]),cex.axis=huge)
  for (i in 1:length(plotmix)){ points(jitmat[,i],freqmat[,i][sel.mask],col=my.colors[i],type = 'p',pch=16)}
  for (i in 1:length(plotmix)){ points(jitmat[,i],freqmat[,i][sel.mask],col=1,type = 'p')}
  for (i in c(1:length(which(sel.mask==T)))){lines(jitmat[i,],freqmat[sel.mask,][i,],col=1)}
  for (i in 1:length(plotmix)) {legend("topright",paste("w=",plotmix, sep=""),pch=16,col=my.colors,bty="n")}
  for (i in 1:length(plotmix)) {legend("topright",paste("w=",plotmix, sep=""),pch=1,col=1,bty="n")}
  abline(v=c(1:length(which(sel.mask==T))), col=grDevices::gray(0.7))
}



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
