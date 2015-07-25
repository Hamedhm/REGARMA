search.arima <- function(x, d=NA, D=NA, max.p=5, max.q=5,
                         max.P=2, max.Q=2, max.order=5, stationary=FALSE, ic=c("aic","aicc","bic"),
                         trace=FALSE,approximation=FALSE,xreg=NULL,offset=offset,allowdrift=TRUE,
                         parallel=FALSE, num.cores=NULL)
{
  require(tseries)
  #dataname <- substitute(x)
  ic <- match.arg(ic)
  m <- frequency(x)
  
  oldwarn <- options()$warn
  options(warn=-1)
  on.exit(options(warn=oldwarn))
  
  if(allowdrift)
    maxK <- (d+D <= 1)
  else
    maxK <- ((d+D) == 0)
  
  # Choose model orders
  #Serial - technically could be combined with the code below
  if (parallel==FALSE)
  {
    best.ic <- 1e20
    for(i in 0:max.p)
    {
      for(j in 0:max.q)
      {
        for(I in 0:max.P)
        {
          for(J in 0:max.Q)
          {
            if(i+j+I+J <= max.order)
            {
              for(K in 0:maxK)
              {
                fit <- myarima(x,order=c(i,d,j),seasonal=c(I,D,J),constant=(K==1),trace=trace,ic=ic,approximation=approximation,offset=offset,xreg=xreg)
                if(fit$ic < best.ic)
                {
                  best.ic <- fit$ic
                  bestfit <- fit
                  constant <- (K==1)
                }
              }
            }
          }
        }
      }
    }
  } else
    ################################################################################
  # Parallel
  if (parallel==TRUE){
    
    to.check <- WhichModels(max.p, max.q, max.P, max.Q, maxK)
    
    par.all.arima <- function(l){
      .tmp <- UndoWhichModels(l)
      i <- .tmp[1]; j <- .tmp[2]; I <- .tmp[3]; J <- .tmp[4]; K <- .tmp[5]==1
      
      if (i+j+I+J <= max.order){
        fit <- myarima(x,order=c(i,d,j),seasonal=c(I,D,J),constant=(K==1),trace=trace,ic=ic,approximation=approximation,offset=offset,xreg=xreg)
      }
      if (exists("fit")){
        return(cbind(fit, K))
      } else return(NULL)
    }
    
    if(is.null(num.cores)) {
      num.cores <- detectCores()
    }
    
    # clusterApplyLB() for Windows, mclapply() for POSIX
    if (Sys.info()[1] == "Windows"){
      cl <- makeCluster(num.cores)
      all.models <- parLapply(cl=cl, X=to.check, fun=par.all.arima)
      stopCluster(cl=cl)
    } else all.models <- mclapply(X=to.check, FUN=par.all.arima, mc.cores=num.cores)
    
    # Removing null elements
    all.models <- all.models[!sapply(all.models, is.null)]
    
    # Choosing best model
    best.ic <- 1e20
    for (i in 1:length(all.models)){
      if(!is.null(all.models[[i]][, 1]$ic) && all.models[[i]][, 1]$ic < best.ic){
        bestfit <- all.models[[i]][, 1]
        best.ic <- bestfit$ic
        constant <- unlist(all.models[[i]][1, 2])
      }
    }
    class(bestfit) <- "Arima"
  }
  ################################################################################
  if(exists("bestfit"))
  {
    # Refit using ML if approximation used for IC
    if(approximation)
    {
      #constant <- length(bestfit$coef) - ncol(xreg) > sum(bestfit$arma[1:4])
      newbestfit <- myarima(x,order=bestfit$arma[c(1,6,2)],
                            seasonal=bestfit$arma[c(3,7,4)],constant=constant,ic,trace=FALSE,approximation=FALSE,xreg=xreg)
      if(newbestfit$ic > 1e19)
      {
        options(warn=oldwarn)
        warning("Unable to fit final model using maximum likelihood. AIC value approximated")
      }
      else
        bestfit <- newbestfit
    }
  }
  else
    stop("No ARIMA model able to be estimated")
  
  bestfit$x <- x
  bestfit$series <- deparse(substitute(x))
  bestfit$ic <- NULL
  bestfit$call <- match.call()
  #bestfit$call$data <- dataname
  #    bestfit$xreg <- xreg
  
  if(trace)
    cat("\n\n")
  
  return(bestfit)
}


is.constant <- function(x)
{
  y <- rep(x[1],length(x))
  isequal <- all.equal(c(x),y)
  return(isequal==TRUE)
}


ndiffs <- function(x,alpha=0.05,test=c("kpss","adf","pp"))
{
  test <- match.arg(test)
  #require(tseries)
  x <- c(na.omit(c(x)))
  d <- 0
  
  if(is.constant(x))
    return(d)
  
  oldwarn <- options(warn=-1)
  on.exit(options(warn=oldwarn$warn))
  if(test=="kpss")
    dodiff <- tseries::kpss.test(x)$p.value < alpha
  else if(test=="adf")
    dodiff <- tseries::adf.test(x)$p.value > alpha
  else if(test=="pp")
    dodiff <- tseries::pp.test(x)$p.value > alpha
  else
    stop("This shouldn't happen")
  if(is.na(dodiff))
  {
    return(d)
  }
  while(dodiff & d<2)
  {
    d <- d+1
    x <- diff(x)
    if(is.constant(x))
      return(d)
    if(test=="kpss")
      dodiff <- tseries::kpss.test(x)$p.value < alpha
    else if(test=="adf")
      dodiff <- tseries::adf.test(x)$p.value > alpha
    else if(test=="pp")
      dodiff <- tseries::pp.test(x)$p.value > alpha
    else
      stop("This shouldn't happen")
    if(is.na(dodiff))
      return(d-1)
  }
  return(d)
}


auto.arima <- function(x, d=NA, D=NA, max.p=5, max.q=5,
                       max.P=2, max.Q=2, max.order=5,
                       start.p=2, start.q=2, start.P=1, start.Q=1,
                       stationary=FALSE, seasonal=TRUE, ic=c("aicc","aic","bic"),
                       stepwise=TRUE, trace=FALSE,
                       approximation=(length(x)>100 | frequency(x)>12), xreg=NULL,
                       test=c("kpss","adf","pp"), seasonal.test=c("ocsb","ch"),
                       allowdrift=TRUE,lambda=NULL,
                       parallel=FALSE, num.cores=NULL)
{
  # Only non-stepwise parallel implemented so far.
  if (stepwise==TRUE & parallel==TRUE)
  {
    warning("Parallel computer is only implemented when stepwise=FALSE, the model will be fit in serial.")
    parallel <- FALSE
  }
  
  series <- deparse(substitute(x))
  
  # Check for constant data
  if(is.constant(x))
  {
    fit <- Arima(x,order=c(0,0,0),fixed=mean(x,na.rm=TRUE))
    fit$x <- x
    fit$series <- series
    fit$call <- match.call()
    fit$call$x <- data.frame(x=x)
    return(fit)
  }
  ic <- match.arg(ic)
  test <- match.arg(test)
  seasonal.test <- match.arg(seasonal.test)
  
  # Only consider non-seasonal models
  if(seasonal)
    m <- frequency(x)
  else
    m <- 1
  if(m < 1)
  {
    #warning("I can't handle data with frequency less than 1. Seasonality will be ignored.")
    m <- 1
  }
  max.p<-ifelse(max.p <= floor(length(x)/3), max.p, floor(length(x)/3))
  max.q<-ifelse(max.q <= floor(length(x)/3), max.q, floor(length(x)/3))
  max.P<-ifelse(max.P <= floor((length(x)/3)/m), max.P, floor((length(x)/3)/m))
  max.Q<-ifelse(max.Q <= floor((length(x)/3)/m), max.Q, floor((length(x)/3)/m))
  
  orig.x <- x
  if(!is.null(lambda))
    x <- BoxCox(x,lambda)
  
  # Choose order of differencing
  if(!is.null(xreg))
  {
    nmxreg <- deparse(substitute(xreg))
    xreg <- as.matrix(xreg)
    if(ncol(xreg)==1 & length(nmxreg) > 1)
      nmxreg <- "xreg"
    if (is.null(colnames(xreg)))
      colnames(xreg) <- if (ncol(xreg) == 1) nmxreg
    else paste(nmxreg, 1:ncol(xreg), sep = "")
    j <- !is.na(x) & !is.na(rowSums(xreg))
    xx <- x
    xx[j] <- residuals(lm(x ~ xreg))
  }
  else
    xx <- x
  if(stationary)
    d <- D <- 0
  if(m == 1)
    D <- max.P <- max.Q <- 0
  else if(is.na(D))
    D <- nsdiffs(xx, m=m, test=seasonal.test)
  if(D > 0)
    dx <- diff(xx,differences=D,lag=m)
  else
    dx <- xx
  if(is.na(d))
    d <- ndiffs(dx,test=test)
  if(d>0)
    dx <- diff(dx,differences=d,lag=1)
  
  if(is.constant(dx))
  {
    if(D>0)
      fit <- Arima(x,order=c(0,d,0),seasonal=list(order=c(0,D,0),period=m), fixed=mean(dx,na.rm=TRUE), include.constant=TRUE)
    else if(d < 2)
      fit <- Arima(x,order=c(0,d,0),fixed=mean(dx,na.rm=TRUE),include.constant=TRUE)
    else
      stop("Data follow a simple polynomial and are not suitable for ARIMA modelling.")
    fit$x <- x
    fit$series <- series
    fit$call <- match.call()
    fit$call$x <- data.frame(x=x)
    return(fit)
  }
  
  
  if(m > 1)
  {
    if(max.P > 0)
      max.p <- min(max.p, m-1)
    if(max.Q > 0)
      max.q <- min(max.q, m-1)
  }
  
  # Find constant offset for AIC calculation using simple AR(1) model
  if(approximation)
  {
    if(D==0)
      fit <- try(arima(x,order=c(1,d,0),xreg=xreg))
    else
      fit <- try(arima(x,order=c(1,d,0),seasonal=list(order=c(0,D,0),period=m,xreg=xreg)))
    if(class(fit) != "try-error")
      offset <- -2*fit$loglik - length(x)*log(fit$sigma2)
    else
    {
      warning("Unable to calculate AIC offset")
      offset <- 0
    }
  }
  else
    offset <- 0
  
  if(!stepwise)
  {
    bestfit <- search.arima(x,d,D,max.p,max.q,max.P,max.Q,max.order,stationary,ic,trace,approximation,xreg=xreg,offset=offset,allowdrift=allowdrift,parallel=parallel, num.cores=num.cores)
    bestfit$call <- match.call()
    bestfit$call$x <- data.frame(x=x)
    bestfit$lamba <- lambda
    bestfit$x <- orig.x
    bestfit$series <- series
    bestfit$fitted <- fitted(bestfit)
    if(!is.null(lambda))
    {
      bestfit$fitted <- InvBoxCox(bestfit$fitted,lambda)
      bestfit$lambda <- lambda
    }
    return(bestfit)
  }
  
  # Starting model
  p <- start.p <- min(start.p,max.p)
  q <- start.q <- min(start.q,max.q)
  P <- start.P <- min(start.P,max.P)
  Q <- start.Q <- min(start.Q,max.Q)
  if(allowdrift)
    constant <- (d+D <= 1)
  else
    constant <- ((d+D) == 0)
  
  results <- matrix(NA,nrow=100,ncol=8)
  
  oldwarn <- options()$warn
  options(warn=-1)
  on.exit(options(warn=oldwarn))
  
  bestfit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
  results[1,] <- c(p,d,q,P,D,Q,constant,bestfit$ic)
  # Null model
  fit <- myarima(x,order=c(0,d,0),seasonal=c(0,D,0),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
  results[2,] <- c(0,d,0,0,D,0,constant,fit$ic)
  if(fit$ic < bestfit$ic)
  {
    bestfit <- fit
    p <- q <- P <- Q <- 0
  }
  # Basic AR model
  if(max.p > 0 | max.P > 0)
  {
    fit <- myarima(x,order=c(max.p>0,d,0),seasonal=c((m>1)&(max.P>0),D,0),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
    results[3,] <- c(1,d,0,m>1,D,0,constant,fit$ic)
    if(fit$ic < bestfit$ic)
    {
      bestfit <- fit
      p <- (max.p>0)
      P <- (m>1) & (max.P>0)
      q <- Q <- 0
    }
  }
  # Basic MA model
  if(max.q > 0 | max.Q > 0)
  {
    fit <- myarima(x,order=c(0,d,max.q>0),seasonal=c(0,D,(m>1)&(max.Q>0)),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
    results[4,] <- c(0,d,1,0,D,m>1,constant,fit$ic)
    if(fit$ic < bestfit$ic)
    {
      bestfit <- fit
      p <- P <- 0
      Q <- (m>1) & (max.Q>0)
      q <- (max.q>0)
    }
  }
  k <- 4
  
  startk <- 0
  while(startk < k & k < 94)
  {
    startk <- k
    if(P > 0 & newmodel(p,d,q,P-1,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q),seasonal=c(P-1,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q,P-1,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        P <- (P-1)
      }
    }
    if(P < max.P & newmodel(p,d,q,P+1,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q),seasonal=c(P+1,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q,P+1,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        P <- (P+1)
      }
    }
    if(Q > 0 & newmodel(p,d,q,P,D,Q-1,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q-1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q,P,D,Q-1,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        Q <- (Q-1)
      }
    }
    if(Q < max.Q & newmodel(p,d,q,P,D,Q+1,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q+1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q,P,D,Q+1,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        Q <- (Q+1)
      }
    }
    if(Q > 0 & P > 0 & newmodel(p,d,q,P-1,D,Q-1,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q),seasonal=c(P-1,D,Q-1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q,P-1,D,Q-1,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        Q <- (Q-1)
        P <- (P-1)
      }
    }
    if(Q < max.Q & P < max.P & newmodel(p,d,q,P+1,D,Q+1,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q),seasonal=c(P+1,D,Q+1),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q,P+1,D,Q+1,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        Q <- (Q+1)
        P <- (P+1)
      }
    }
    
    if(p > 0 & newmodel(p-1,d,q,P,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p-1,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p-1,d,q,P,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        p <- (p-1)
      }
    }
    if(p < max.p & newmodel(p+1,d,q,P,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p+1,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p+1,d,q,P,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        p <- (p+1)
      }
    }
    if(q > 0 & newmodel(p,d,q-1,P,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q-1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q-1,P,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        q <- (q-1)
      }
    }
    if(q < max.q & newmodel(p,d,q+1,P,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p,d,q+1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p,d,q+1,P,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        q <- (q+1)
      }
    }
    if(q > 0 & p > 0 & newmodel(p-1,d,q-1,P,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p-1,d,q-1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p-1,d,q-1,P,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        q <- (q-1)
        p <- (p-1)
      }
    }
    if(q < max.q & p < max.p & newmodel(p+1,d,q+1,P,D,Q,constant,results[1:k,]))
    {
      k <- k + 1
      fit <- myarima(x,order=c(p+1,d,q+1),seasonal=c(P,D,Q),constant=constant,ic,trace,approximation,offset=offset,xreg=xreg)
      results[k,] <- c(p+1,d,q+1,P,D,Q,constant,fit$ic)
      if(fit$ic < bestfit$ic)
      {
        bestfit <- fit
        q <- (q+1)
        p <- (p+1)
      }
    }
    if(allowdrift | (d+D)==0)
    {
      if(newmodel(p,d,q,P,D,Q,!constant,results[1:k,]))
      {
        k <- k + 1
        fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q),constant=!constant,ic,trace,approximation,offset=offset,xreg=xreg)
        results[k,] <- c(p,d,q,P,D,Q,!constant,fit$ic)
        if(fit$ic < bestfit$ic)
        {
          bestfit <- fit
          constant <- !constant
        }
      }
    }
  }
  
  # Refit using ML if approximation used for IC
  if(approximation)
  {
    #constant <- length(bestfit$coef) > sum(bestfit$arma[1:4])
    newbestfit <- myarima(x,order=bestfit$arma[c(1,6,2)],
                          seasonal=bestfit$arma[c(3,7,4)],constant=constant,ic,trace=FALSE,approximation=FALSE,xreg=xreg)
    if(newbestfit$ic > 1e19)
    {
      options(warn=oldwarn)
      warning("Unable to fit final model using maximum likelihood. AIC value approximated")
    }
    else
      bestfit <- newbestfit
  }
  
  if(bestfit$ic > 1e19)
  {
    cat("\n")
    stop("No suitable ARIMA model found")
  }
  
  # Return best fit
  
  bestfit$x <- orig.x
  bestfit$series <- series
  bestfit$ic <- NULL
  bestfit$call <- match.call()
  bestfit$call$x <- data.frame(x=x)
  bestfit$lambda <- lambda
  #bestfit$fitted <- fitted(bestfit)
  
  if(trace)
    cat("\n\n Best model:",arima.string(bestfit),"\n\n")
  
  return(bestfit)
}


# Calls arima from stats package and adds data to the returned object
# Also allows refitting to new data
# and drift terms to be included.
myarima <- function(x, order = c(0, 0, 0), seasonal = c(0, 0, 0), constant=TRUE, ic="aic", trace=FALSE,approximation=FALSE,offset=0,xreg=NULL)
{
  n <- length(x)
  m <- frequency(x)
  use.season <- (sum(seasonal)>0) & m>0
  diffs <- order[2]+seasonal[2]
  if(approximation)
    method <- "CSS"
  else
    method <- "CSS-ML"
  if(diffs==1 & constant)
  {
    xreg <- cbind(drift=1:length(x),xreg)
    if(use.season)
      fit <- try(stats::arima(x=x,order=order,seasonal=list(order=seasonal,period=m),xreg=xreg,method=method),silent=TRUE)
    else
      fit <- try(stats::arima(x=x,order=order,xreg=xreg,method=method),silent=TRUE)
  }
  else
  {
    if(use.season)
      fit <- try(stats::arima(x=x,order=order,seasonal=list(order=seasonal,period=m),include.mean=constant,method=method,xreg=xreg),silent=TRUE)
    else
      fit <- try(stats::arima(x=x,order=order,include.mean=constant,method=method,xreg=xreg),silent=TRUE)
  }
  if(is.null(xreg))
    nxreg <- 0
  else
    nxreg <- ncol(as.matrix(xreg))
  if(class(fit) != "try-error")
  {
    nstar <- n - order[2] - seasonal[2]*m
    if(diffs==1 & constant)
    {
      #fitnames <- names(fit$coef)
      #fitnames[length(fitnames)-nxreg] <- "drift"
      #names(fit$coef) <- fitnames
      fit$xreg <- xreg
    }
    npar <- length(fit$coef) + 1
    if(approximation)
      fit$aic <- offset + nstar * log(fit$sigma2) + 2 * npar
    if(!is.na(fit$aic))
    {
      fit$bic <- fit$aic + npar*(log(nstar) - 2)
      fit$aicc <- fit$aic + 2*npar*(nstar/(nstar-npar-1) - 1)
      fit$ic <- switch(ic,bic=fit$bic,aic=fit$aic,aicc=fit$aicc)
    }
    else
      fit$aic <- fit$bic <- fit$aicc <- fit$ic <- 1e20
    # Check for unit roots
    minroot <- 2
    if(order[1] + seasonal[1] > 0)
    {
      testvec <- fit$model$phi
      last.nonzero <- max(which(abs(testvec)>1e-8))
      if(last.nonzero > 0)
      {
        testvec <- testvec[1:last.nonzero]
        if(last.nonzero > 48)
          warning("Unable to check for unit roots")
        else
          minroot <- min(minroot,abs(polyroot(c(1,-testvec))))
      }
    }
    if(order[3] + seasonal[3] > 0)
    {
      testvec <- fit$model$theta
      last.nonzero <- max(which(abs(testvec)>1e-8))
      if(last.nonzero > 0)
      {
        testvec <- testvec[1:last.nonzero]
        if(last.nonzero > 48)
          warning("Unable to check for unit roots")
        else
          minroot <- min(minroot,abs(polyroot(c(1,testvec))))
      }
    }
    if(minroot < 1 + 1e-3)
      fit$ic <- 1e20 # Don't like this model
    if(trace)
      cat("\n",arima.string(fit),":",fit$ic)
    fit$xreg <- xreg
    return(fit)
  }
  else
  {
    if(trace)
    {
      cat("\n ARIMA(",order[1],",",order[2],",",order[3],")",sep="")
      if(use.season)
        cat("(",seasonal[1],",",seasonal[2],",",seasonal[3],")[",m,"]",sep="")
      if(constant & (order[2]+seasonal[2] == 0))
        cat(" with non-zero mean")
      else if(constant & (order[2]+seasonal[2] == 1))
        cat(" with drift        ")
      else if(!constant & (order[2]+seasonal[2] == 0))
        cat(" with zero mean    ")
      else
        cat("         ")
      cat(" :",1e20,"*")
    }
    return(list(ic=1e20))
  }
}

newmodel <- function(p,d,q,P,D,Q,constant,results)
{
  n <- nrow(results)
  for(i in 1:n)
  {
    if(identical(c(p,d,q,P,D,Q,constant),results[i,1:7]))
      return(FALSE)
  }
  return(TRUE)
}

arima.string <- function(object)
{
  order <- object$arma[c(1,6,2,3,7,4,5)]
  result <- paste("ARIMA(",order[1],",",order[2],",",order[3],")",sep="")
  if(order[7]>1 & sum(order[4:6]) > 0)
    result <- paste(result,"(",order[4],",",order[5],",",order[6],")[",order[7],"]",sep="")
  if(is.element("constant",names(object$coef)) | is.element("intercept",names(object$coef)))
    result <- paste(result,"with non-zero mean")
  else if(is.element("drift",names(object$coef)))
    result <- paste(result,"with drift        ")
  else if(order[2]==0 & order[5]==0)
    result <- paste(result,"with zero mean    ")
  else
    result <- paste(result,"                  ")
  return(result)
}

# summary.Arima <- function(object,...)
# {
#   print(object)
#   cat("\nTraining set error measures:\n")
#   print(accuracy(object))
# }



# Number of seasonal differences
nsdiffs <- function(x, m=frequency(x), test=c("ocsb", "ch"))
{
  
  if(is.constant(x))
    return(0)
  
  test <- match.arg(test)
  if(m==1)
    stop("Non seasonal data")
  else if(m < 1)
  {
    warning("I can't handle data with frequency less than 1. Seasonality will be ignored.")
    return(0)
  }
  
  if(test=="ch")
    return(CHtest(x,m))
  else
    return(OCSBtest(x,m))
  
}

CHtest <- function(x,m)
{
  chstat <- SD.test(x, m)
  crit.values <- c(0.4617146,0.7479655,1.0007818,1.2375350,1.4625240,1.6920200,1.9043096,2.1169602,
                   2.3268562,2.5406922,2.7391007)
  if(m <= 12)
    D <- as.numeric(chstat > crit.values[m-1])
  else if (m == 24)
    D <- as.numeric(chstat > 5.098624)
  else if (m ==52)
    D <- as.numeric(chstat > 10.341416)
  else if (m ==365)
    D <- as.numeric(chstat > 65.44445)
  else
    D <- as.numeric(chstat > 0.269 * m^(0.928))
  return(D)
}

# Return critical values for OCSB test at 5% level
# Approximation based on extensive simulations.
calcOCSBCritVal <- function(seasonal.period)
{
  log.m <- log(seasonal.period)
  return(-0.2937411*exp(-0.2850853*(log.m-0.7656451)+(-0.05983644)*((log.m-0.7656451)^2))-1.652202)
}


OCSBtest <- function(time.series, period)
{
  if(length(time.series) < (2*period+5))
  {
    warning("Time series too short for seasonal differencing")
    return(0)
  }
  
  seas.diff.series <- diff(time.series, lag = period, differences=1)
  diff.series <- diff(seas.diff.series, lag = 1, differences=1)
  
  y.one <- time.series[2:length(time.series)]
  y.one <- diff(y.one, lag=period, differences=1)
  
  y.two <- time.series[(1+period):length(time.series)]
  y.two <- diff(y.two, lag=1, differences=1)
  
  y.one <- y.one[(1+period):(length(y.one)-1)]
  y.two <- y.two[2:(length(y.two)-period)]
  diff.series <- diff.series[(1+period+1):(length(diff.series))]
  contingent.series <- diff.series
  
  x.reg <- cbind(y.one, y.two)
  diff.series <- ts(data = diff.series, frequency=period)
  ##Turn off warnings
  old.warning.level <- options()$warn
  options(warn=-1)
  regression <- try(Arima(diff.series, order=c(3,0,0), seasonal=list(order=c(1,0,0),period=period), xreg=x.reg), silent=TRUE)
  
  if(class(regression) == "try-error" | tryCatch(any(is.nan(sqrt(diag(regression$var.coef)))), error=function(e) TRUE))
  {
    regression <- try(Arima(diff.series, order=c(3,0,0), seasonal=list(order=c(0,0,0),period=period), xreg=x.reg), silent=TRUE)
    
    if(class(regression) == "try-error" | tryCatch(any(is.nan(sqrt(diag(regression$var.coef)))), error=function(e) TRUE))
    {
      regression <- try(Arima(diff.series, order=c(2,0,0), seasonal=list(order=c(0,0,0),period=period), xreg=x.reg), silent=TRUE)
      
      if(class(regression) == "try-error" | tryCatch(any(is.nan(sqrt(diag(regression$var.coef)))), error=function(e) TRUE))
      {
        regression <- try(Arima(diff.series, order=c(1,0,0), seasonal=list(order=c(0,0,0),period=period), xreg=x.reg), silent=TRUE)
        
        if(class(regression) == "try-error" | tryCatch(any(is.nan(sqrt(diag(regression$var.coef)))), error=function(e) TRUE))
        {
          regression <- try(lm(contingent.series ~ y.one + y.two - 1, na.action=NULL), silent=TRUE)
          reg.summary <- summary(regression)
          reg.coefs <- reg.summary$coefficients
          t.two.pos <- grep("t.two", rownames(reg.coefs), fixed = TRUE)
          if(length(t.two.pos) != 0)
            t.two <- reg.coefs[t.two.pos,3]
          else
            t.two <- NA
          
          ###Re-enable warnings
          options(warn=old.warning.level)
          
          if((is.nan(t.two)) | (is.infinite(t.two)) | (is.na(t.two)) | (class(regression) == "try-error"))
            return(1)
          else
            return(as.numeric(t.two >= calcOCSBCritVal(period)))
        }
      }
    }
  }
  
  se <- sqrt(diag(regression$var.coef))
  t.two <- regression$coef[names(regression$coef)=="y.two"]/se[names(se)=="y.two"]
  ###Re-enable warnings
  options(warn=old.warning.level)
  return(as.numeric(t.two >= calcOCSBCritVal(period)))
}



# Set up seasonal dummies using Fourier series
SeasDummy <- function(x)
{
  n <- length(x)
  m <- frequency(x)
  if(m==1)
    stop("Non-seasonal data")
  tt <- 1:n
  fmat <- matrix(NA,nrow=n,ncol=2*m)
  for(i in 1:m)
  {
    fmat[,2*i] <- sin(2*pi*i*tt/m)
    fmat[,2*(i-1)+1] <- cos(2*pi*i*tt/m)
  }
  return(fmat[,1:(m-1)])
}

# CANOVA-HANSEN TEST
# Largely based on uroot package code for CH.test()
SD.test <- function (wts, s=frequency(wts))
{
  if(any(is.na(wts)))
    stop("Series contains missing values. Please choose order of seasonal differencing manually.")
  if(s==1)
    stop("Not seasonal data")
  t0 <- start(wts)
  N <- length(wts)
  if(N <= s)
    stop("Insufficient data")
  frec <- rep(1, as.integer((s+1)/2))
  ltrunc <- round(s * (N/100)^0.25)
  R1 <- as.matrix(SeasDummy(wts))
  lmch <- lm(wts ~ R1, na.action=na.exclude)   # run the regression : y(i)=mu+f(i)'gamma(i)+e(i)
  Fhat <- Fhataux <- matrix(nrow=N, ncol=s-1)
  for (i in 1:(s-1))
    Fhataux[, i] <- R1[,i] * residuals(lmch)
  for (i in 1:N)
  {
    for (n in 1:(s - 1))
      Fhat[i, n] <- sum(Fhataux[1:i, n])
  }
  wnw <- 1 - seq(1, ltrunc, 1)/(ltrunc + 1)
  Ne <- nrow(Fhataux)
  Omnw <- 0
  for (k in 1:ltrunc)
    Omnw <- Omnw + (t(Fhataux)[, (k + 1):Ne] %*% Fhataux[1:(Ne - k), ]) * wnw[k]
  Omfhat <- (crossprod(Fhataux) + Omnw + t(Omnw))/Ne
  sq <- seq(1, s-1, 2)
  frecob <- rep(0,s - 1)
  for (i in 1:length(frec))
  {
    if (frec[i] == 1 && i == as.integer(s/2))
      frecob[sq[i]] <- 1
    if (frec[i] == 1 && i < as.integer(s/2))
      frecob[sq[i]] <- frecob[sq[i] + 1] <- 1
  }
  a <- length(which(frecob == 1))
  A <- matrix(0, nrow=s - 1, ncol=a)
  j <- 1
  for (i in 1:(s - 1)) if (frecob[i] == 1)
  {
    A[i, j] <- 1
    ifelse(frecob[i] == 1, j <- j + 1, j <- j)
  }
  tmp <- t(A) %*% Omfhat %*% A
  problems <- (min(svd(tmp)$d) < .Machine$double.eps)
  if(problems)
    stL <- 0
  else
    stL <- (1/N^2) * sum(diag(solve(tmp, tol=1e-25) %*% t(A) %*% t(Fhat) %*% Fhat %*% A))
  return(stL)
}


forecast.Arima <- function (object, h=ifelse(object$arma[5] > 1, 2 * object$arma[5], 10),
                            level=c(80, 95), fan=FALSE, xreg=NULL, lambda=object$lambda,  bootstrap=FALSE, npaths=5000,...)
{
  #    use.constant <- is.element("constant",names(object$coef))
  use.drift <- is.element("drift", names(object$coef))
  x <- object$x <- getResponse(object)
  usexreg <- (!is.null(xreg) | use.drift | is.element("xreg",names(object)))# | use.constant)
  #    if(use.constant)
  #        xreg <- as.matrix(rep(1,h))
  if(!is.null(xreg))
  {
    xreg <- as.matrix(xreg)
    h <- nrow(xreg)
  }
  if (use.drift)
  {
    n <- length(x)
    if(!is.null(xreg))
      xreg <- cbind((1:h)+n,xreg)
    else
      xreg <- as.matrix((1:h)+n)
  }
  if(usexreg)
  {
    if(is.null(xreg))
      stop("No regressors provided")
    if(!is.null(object$xreg))
      object$call$xreg <- object$xreg
    else # object from arima() rather than Arima()
    {
      xr <- object$call$xreg
      object$call$xreg <- if (!is.null(xr))
        eval.parent(xr)
      else NULL
    }
    if(ncol(xreg) != ncol(object$call$xreg))
      stop("Number of regressors does not match fitted model")
    pred <- predict(object, n.ahead=h, newxreg=xreg)
  }
  else
    pred <- predict(object, n.ahead=h)
  
  if(bootstrap) # Recompute se using simulations
  {
    sim <- matrix(NA,nrow=npaths,ncol=h)
    for(i in 1:npaths)
      sim[i,] <- simulate(object, nsim=h, bootstrap=TRUE, xreg=xreg, lambda=NULL)
    pred$se <- apply(sim,2,sd)
  }
  
  # Fix time series characteristics if there are missing values at end of series.
  if(!is.null(x))
  {
    tspx <- tsp(x)
    nx <- max(which(!is.na(x)))
    if(nx != length(x))
    {
      tspx[2] <- time(x)[nx]
      start.f <- tspx[2]+1/tspx[3]
      pred$pred <- ts(pred$pred,frequency=tspx[3],start=start.f)
      pred$se <- ts(pred$se,frequency=tspx[3],start=start.f)
    }
  }
  
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  
  nint <- length(level)
  lower <- matrix(NA, ncol=nint, nrow=length(pred$pred))
  upper <- lower
  for (i in 1:nint)
  {
    qq <- qnorm(0.5 * (1 + level[i]/100))
    lower[, i] <- pred$pred - qq * pred$se
    upper[, i] <- pred$pred + qq * pred$se
  }
  colnames(lower)=colnames(upper)=paste(level, "%", sep="")
  method <- arima.string(object)
  fits <- fitted(object)
  if(!is.null(lambda))
  {
    pred$pred <- InvBoxCox(pred$pred,lambda)
    lower <- InvBoxCox(lower,lambda)
    upper <- InvBoxCox(upper,lambda)
  }
  return(structure(list(method=method, model=object, level=level,
                        mean=pred$pred, lower=lower, upper=upper, x=x,
                        xname=deparse(substitute(x)), fitted=fits, residuals=residuals(object)),
                   class="forecast"))
}


forecast.ar <- function(object,h=10,level=c(80,95),fan=FALSE, lambda=NULL,  bootstrap=FALSE, npaths=5000,...)
{
  pred <- predict(object,n.ahead=h)
  if(bootstrap) # Recompute se using simulations
  {
    sim <- matrix(NA,nrow=npaths,ncol=h)
    for(i in 1:npaths)
      sim[i,] <- simulate(object, nsim=h, bootstrap=TRUE)
    pred$se <- apply(sim,2,sd)
  }
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nint <- length(level)
  lower <- matrix(NA,ncol=nint,nrow=length(pred$pred))
  upper <- lower
  for(i in 1:nint)
  {
    qq <- qnorm(0.5*(1+level[i]/100))
    lower[,i] <- pred$pred - qq*pred$se
    upper[,i] <- pred$pred + qq*pred$se
  }
  colnames(lower)=colnames(upper)=paste(level,"%",sep="")
  method <- paste("AR(",object$order,")",sep="")
  x <- getResponse(object)
  f <- frequency(x)
  res <- ts(object$resid[-(1:object$order)],start=tsp(x)[1]+object$order/f,frequency=f)
  fits <- x-res
  
  if(!is.null(lambda))
  {
    pred$pred <- InvBoxCox(pred$pred,lambda)
    lower <- InvBoxCox(lower,lambda)
    upper <- InvBoxCox(upper,lambda)
    fits <- InvBoxCox(fits,lambda)
    x <- InvBoxCox(x,lambda)
  }
  
  return(structure(list(method=method,model=object,level=level,mean=pred$pred,lower=lower,upper=upper,
                        x=x, xname=deparse(substitute(x)), fitted=fits,residuals=res)
                   ,class="forecast"))
}

# Extract errors from ARIMA model (as distinct from residuals)
arima.errors <- function(z)
{
  if(!is.list(z))
    stop("z must be a list")
  x <- getResponse(z)
  if(!is.element("xreg",names(z)))
  {
    if(!is.element("xreg",names(z$coef)))
      return(x)
    else
      xreg <- eval.parent(z$coef$xreg)
  }
  else
    xreg <- z$xreg
  norder <- sum(z$arma[1:4])
  if(is.element("intercept",names(z$coef)))
    xreg <- cbind(rep(1,length(x)),xreg)
  return(ts(x - xreg %*% as.matrix(z$coef[(norder+1):length(z$coef)]),frequency=frequency(x),start=start(x)))
}

# Return one-step fits
fitted.Arima <- function(object,...)
{
  x <- getResponse(object)
  if(is.null(x))
  {
    #warning("Fitted values are unavailable due to missing historical data")
    return(NULL)
  }
  if(is.null(object$lambda))
    return(x - object$residuals)
  else
    return(InvBoxCox(BoxCox(x,object$lambda) - object$residuals, object$lambda))
}

# Calls arima from stats package and adds data to the returned object
# Also allows refitting to new data
# and drift terms to be included.
Arima <- function(x, order=c(0, 0, 0),
                  seasonal=list(order=c(0, 0, 0), period=NA),
                  xreg=NULL, include.mean=TRUE, include.drift=FALSE, include.constant, lambda=model$lambda,
                  transform.pars=TRUE,
                  fixed=NULL, init=NULL, method=c("CSS-ML", "ML", "CSS"),
                  n.cond, optim.control=list(), kappa=1e6, model=NULL)
{
  # Remove outliers near ends
  #j <- time(x)
  #x <- na.contiguous(x)
  #if(length(j) != length(x))
  #    warning("Missing values encountered. Using longest contiguous portion of time series")
  
  series <- deparse(substitute(x))
  
  origx <- x
  if(!is.null(lambda))
    x <- BoxCox(x,lambda)
  
  if (!is.null(xreg))
  {
    nmxreg <- deparse(substitute(xreg))
    xreg <- as.matrix(xreg)
    if (is.null(colnames(xreg)))
      colnames(xreg) <- if (ncol(xreg) == 1) nmxreg else paste(nmxreg, 1:ncol(xreg), sep="")
  }
  
  if(!missing(include.constant))
  {
    if(include.constant)
    {
      include.mean <- TRUE
      if((order[2] + seasonal$order[2]) == 1)
        include.drift <- TRUE
    }
    else
    {
      include.mean <- include.drift <- FALSE
    }
  }
  
  if(!is.null(model))
  {
    tmp <- arima2(x,model,xreg=xreg)
    xreg <- tmp$xreg
  }
  else
  {
    if(include.drift)
    {
      drift <- 1:length(x)
      xreg <- cbind(drift=drift,xreg)
    }
    if(is.null(xreg))
      tmp <- stats::arima(x=x,order=order,seasonal=seasonal,include.mean=include.mean,
                          transform.pars=transform.pars,fixed=fixed,init=init,method=method,n.cond=n.cond,optim.control=optim.control,kappa=kappa)
    else
      tmp <- stats::arima(x=x,order=order,seasonal=seasonal,xreg=xreg,include.mean=include.mean,
                          transform.pars=transform.pars,fixed=fixed,init=init,method=method,n.cond=n.cond,optim.control=optim.control,kappa=kappa)
  }
  tmp$series <- series
  tmp$xreg <- xreg
  tmp$call <- match.call()
  tmp$lambda <- lambda
  tmp$x <- origx
  
  return(tmp)
}

# Refits the model to new data x
arima2 <- function (x, model, xreg)
{
  use.drift <- is.element("drift",names(model$coef))
  use.intercept <- is.element("intercept",names(model$coef))
  use.xreg <- is.element("xreg",names(model$call))
  if(use.drift)
  {
    driftmod <- lm(model$xreg[,"drift"] ~ I(time(model$x)))
    newxreg <- driftmod$coeff[1] + driftmod$coeff[2]*time(x)
    if(!is.null(xreg))
      xreg[,"drift"] <- newxreg
    else
      xreg <- as.matrix(data.frame(drift=newxreg))
    use.xreg <- TRUE
  }
  
  if(model$arma[5]>1 & sum(abs(model$arma[c(3,4,7)]))>0) # Seasonal model
  {
    if(use.xreg)
      refit <- Arima(x,order=model$arma[c(1,6,2)],seasonal=list(order=model$arma[c(3,7,4)],period=model$arma[5]),
                     fixed=model$coef,include.mean=use.intercept,xreg=xreg)
    else
      refit <- Arima(x,order=model$arma[c(1,6,2)],seasonal=list(order=model$arma[c(3,7,4)],period=model$arma[5]),
                     fixed=model$coef,include.mean=use.intercept)
  }
  else if(length(model$coef)>0) # Nonseasonal model with some parameters
  {
    if(use.xreg)
      refit <- Arima(x,order=model$arma[c(1,6,2)],fixed=model$coef,xreg=xreg,include.mean=use.intercept)
    else
      refit <- Arima(x,order=model$arma[c(1,6,2)],fixed=model$coef,include.mean=use.intercept)
  }
  else # No parameters
    refit <- Arima(x,order=model$arma[c(1,6,2)],include.mean=FALSE)
  
  refit$var.coef <- matrix(0,length(refit$coef),length(refit$coef))
  if(use.xreg) # Why is this needed?
    refit$xreg <- xreg
  return(refit)
}

# Modified version of function in stats package

print.Arima <- function (x, digits=max(3, getOption("digits") - 3), se=TRUE,
                         ...)
{
  cat("Series:",x$series,"\n")
  cat(arima.string(x),"\n")
  if(!is.null(x$lambda))
    cat("Box Cox transformation: lambda=",x$lambda,"\n")
  #cat("\nCall:", deparse(x$call, width.cutoff=75), "\n", sep=" ")
  #    if(!is.null(x$xreg))
  #    {
  #        cat("\nRegression variables fitted:\n")
  #        xreg <- as.matrix(x$xreg)
  #        for(i in 1:3)
  #            cat("  ",xreg[i,],"\n")
  #        cat("   . . .\n")
  #        for(i in 1:3)
  #            cat("  ",xreg[nrow(xreg)-3+i,],"\n")
  #    }
  if (length(x$coef) > 0) {
    cat("\nCoefficients:\n")
    coef <- round(x$coef, digits=digits)
    if (se && NROW(x$var.coef)) {
      ses <- rep(0, length(coef))
      ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits=digits)
      coef <- matrix(coef, 1L, dimnames=list(NULL, names(coef)))
      coef <- rbind(coef, s.e.=ses)
    }
    print.default(coef, print.gap=2)
  }
  cm <- x$call$method
  if (is.null(cm) || cm != "CSS")
  {
    cat("\nsigma^2 estimated as ", format(x$sigma2, digits=digits),
        ":  log likelihood=", format(round(x$loglik, 2L)),"\n",sep="")
    npar <- length(x$coef) + 1
    nstar <- length(x$residuals) - x$arma[6] - x$arma[7]*x$arma[5]
    bic <- x$aic + npar*(log(nstar) - 2)
    aicc <- x$aic + 2*npar*(nstar/(nstar-npar-1) - 1)
    cat("AIC=", format(round(x$aic, 2L)), sep="")
    cat("   AICc=", format(round(aicc, 2L)), sep="")
    cat("   BIC=", format(round(bic, 2L)), "\n",sep="")
  }
  else cat("\nsigma^2 estimated as ", format(x$sigma2, digits=digits),
           ":  part log likelihood=", format(round(x$loglik, 2)),
           "\n", sep="")
  invisible(x)
}

# Modified version of function in stats package

predict.Arima <- function(object, n.ahead=1, newxreg=NULL, se.fit=TRUE, ...)
{
  myNCOL <- function(x) if (is.null(x))
    0
  else NCOL(x)
  rsd <- object$residuals
  ## LINES ADDED
  if(!is.null(object$xreg))
    object$call$xreg <- object$xreg
  ## END ADDITION
  xr <- object$call$xreg
  xreg <- if (!is.null(xr))
    eval.parent(xr)
  else NULL
  ncxreg <- myNCOL(xreg)
  if (myNCOL(newxreg) != ncxreg)
    stop("'xreg' and 'newxreg' have different numbers of columns: ", ncxreg, " != ", myNCOL(newxreg))
  class(xreg) <- NULL
  xtsp <- tsp(rsd)
  n <- length(rsd)
  arma <- object$arma
  coefs <- object$coef
  narma <- sum(arma[1:4])
  if (length(coefs) > narma) {
    if (names(coefs)[narma + 1] == "intercept") {
      xreg <- cbind(intercept=rep(1, n), xreg)
      newxreg <- cbind(intercept=rep(1, n.ahead), newxreg)
      ncxreg <- ncxreg + 1
    }
    xm <- if (narma == 0)
      drop(as.matrix(newxreg) %*% coefs)
    else drop(as.matrix(newxreg) %*% coefs[-(1:narma)])
  }
  else xm <- 0
  if (arma[2] > 0) {
    ma <- coefs[arma[1] + 1:arma[2]]
    if (any(Mod(polyroot(c(1, ma))) < 1))
      warning("MA part of model is not invertible")
  }
  if (arma[4] > 0) {
    ma <- coefs[sum(arma[1:3]) + 1:arma[4]]
    if (any(Mod(polyroot(c(1, ma))) < 1))
      warning("seasonal MA part of model is not invertible")
  }
  z <- KalmanForecast(n.ahead, object$model)
  pred <- ts(z[[1]] + xm, start=xtsp[2] + deltat(rsd), frequency=xtsp[3])
  if (se.fit) {
    se <- ts(sqrt(z[[2]] * object$sigma2), start=xtsp[2] +
               deltat(rsd), frequency=xtsp[3])
    return(list(pred=pred, se=se))
  }
  else return(pred)
}
## Generic forecast functions
## Part of forecast and demography packages

forecast <- function(object,...) UseMethod("forecast")

forecast.default <- function(object,...) forecast.ts(object,...)

## A function determining the appropriate period, if the data is of unknown period
## Written by Rob Hyndman
find.freq <- function(x)
{
  n <- length(x)
  spec <- spec.ar(c(na.contiguous(x)),plot=FALSE)
  if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) # Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        if(nextmax <= length(spec$freq))
          period <- round(1/spec$freq[nextmax])
        else
          period <- 1
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  
  return(period)
}

forecast.ts <- function(object, h=ifelse(frequency(object)>1, 2*frequency(object), 10), 
                        level=c(80,95), fan=FALSE, robust=FALSE, lambda = NULL, find.frequency = FALSE, ...)
{
  n <- length(object)
  if (find.frequency) {
    object <- ts(object, frequency = find.freq(object))
    obj.freq <- frequency(object)
  } else {
    obj.freq <- frequency(object)
  }
  if(robust)
    object <- tsclean(object, replace.missing=TRUE, lambda = lambda)
  
  if(n > 3)
  {
    if(obj.freq < 13)
      forecast(ets(object,lambda = lambda, ...),h=h,level=level,fan=fan)
    else
      stlf(object,h=h,level=level,fan=fan,lambda = lambda, ...)
  }
  else
    meanf(object,h=h,level=level,fan=fan,lambda = lambda, ...)
}

as.data.frame.forecast <- function(x,...)
{
  nconf <- length(x$level)
  out <- matrix(x$mean, ncol=1)
  ists <- is.ts(x$mean)
  if(ists)
  {
    out <- ts(out)
    attributes(out)$tsp <- attributes(x$mean)$tsp
  }
  names <- c("Point Forecast")
  if (!is.null(x$lower) & !is.null(x$upper) & !is.null(x$level))
  {
    x$upper <- as.matrix(x$upper)
    x$lower <- as.matrix(x$lower)
    for (i in 1:nconf)
    {
      out <- cbind(out, x$lower[, i], x$upper[, i])
      names <- c(names, paste("Lo", x$level[i]), paste("Hi", x$level[i]))
    }
  }
  colnames(out) <- names
  rownames(out) <- time(x$mean)
  # Rest of function borrowed from print.ts(), but with header() omitted
  if(!ists)
    return(as.data.frame(out))
  
  x <- as.ts(out)
  fr.x <- frequency(x)
  calendar <- any(fr.x == c(4, 12)) && length(start(x)) ==  2L
  Tsp <- tsp(x)
  if (is.null(Tsp))
  {
    warning("series is corrupt, with no 'tsp' attribute")
    print(unclass(x))
    return(invisible(x))
  }
  nn <- 1 + round((Tsp[2L] - Tsp[1L]) * Tsp[3L])
  if (NROW(x) != nn)
  {
    warning(gettextf("series is corrupt: length %d with 'tsp' implying %d", NROW(x), nn), domain=NA, call.=FALSE)
    calendar <- FALSE
  }
  if (NCOL(x) == 1)
  {
    if (calendar)
    {
      if (fr.x > 1)
      {
        dn2 <- if (fr.x == 12)
          month.abb
        else if (fr.x == 4)
          c("Qtr1", "Qtr2", "Qtr3", "Qtr4")
        else paste("p", 1L:fr.x, sep="")
        if (NROW(x) <= fr.x && start(x)[1L] == end(x)[1L])
        {
          dn1 <- start(x)[1L]
          dn2 <- dn2[1 + (start(x)[2L] - 2 + seq_along(x))%%fr.x]
          x <- matrix(format(x, ...), nrow=1L, byrow=TRUE,
                      dimnames=list(dn1, dn2))
        }
        else
        {
          start.pad <- start(x)[2L] - 1
          end.pad <- fr.x - end(x)[2L]
          dn1 <- start(x)[1L]:end(x)[1L]
          x <- matrix(c(rep.int("", start.pad), format(x, ...), rep.int("", end.pad)), ncol=fr.x,
                      byrow=TRUE, dimnames=list(dn1, dn2))
        }
      }
      else
      {
        tx <- time(x)
        attributes(x) <- NULL
        names(x) <- tx
      }
    }
    else
      attr(x, "class") <- attr(x, "tsp") <- attr(x, "na.action") <- NULL
  }
  else
  {
    if (calendar && fr.x > 1)
    {
      tm <- time(x)
      t2 <- 1 + round(fr.x * ((tm + 0.001)%%1))
      p1 <- format(floor(zapsmall(tm)))
      rownames(x) <- if (fr.x == 12)
        paste(month.abb[t2], p1, sep=" ")
      else paste(p1, if (fr.x == 4)
        c("Q1", "Q2", "Q3", "Q4")[t2]
        else format(t2), sep=" ")
    }
    else
      rownames(x) <- format(time(x))
    attr(x, "class") <- attr(x, "tsp") <- attr(x, "na.action") <- NULL
  }
  return(as.data.frame(x))
}

print.forecast <- function(x ,...)
{
  print(as.data.frame(x))
}


summary.forecast <- function(object,...)
{
  cat(paste("\nForecast method:",object$method))
  #    cat(paste("\n\nCall:\n",deparse(object$call)))
  cat(paste("\n\nModel Information:\n"))
  print(object$model)
  cat("\nError measures:\n")
  print(accuracy(object))
  if(is.null(object$mean))
    cat("\n No forecasts\n")
  else
  {
    cat("\nForecasts:\n")
    print(object)
  }
}

plotlmforecast <- function(object, plot.conf, shaded, shadecols, col, fcol, pi.col, pi.lty,
                           xlim=NULL, ylim, main, ylab, xlab, ...)
{
  xvar <- attributes(terms(object$model))$term.labels
  if(length(xvar) > 1)
    stop("Forecast plot for regression models only available for a single predictor")
  else if(ncol(object$newdata)==1) # Make sure column has correct name
    colnames(object$newdata) <- xvar
  if(is.null(xlim))
    xlim <- range(object$newdata[,xvar],model.frame(object$model)[,xvar])
  if(is.null(ylim))
    ylim <- range(object$upper,object$lower,fitted(object$model)+residuals(object$model))
  plot(formula(object$model),data=model.frame(object$model), 
       xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,col=col,...)
  abline(object$model)
  nf <- length(object$mean)
  if(plot.conf)
  {
    nint <- length(object$level)
    idx <- rev(order(object$level))
    if(is.null(shadecols))
    {
      #require(colorspace)
      if(min(object$level) < 50) # Using very small confidence levels.
        shadecols <- rev(colorspace::sequential_hcl(100)[object$level])
      else # This should happen almost all the time. Colors mapped to levels.
        shadecols <- rev(colorspace::sequential_hcl(52)[object$level-49])
    }
    if(length(shadecols)==1)
    {
      if(shadecols=="oldstyle") # Default behaviour up to v3.25.
        shadecols <- heat.colors(nint+2)[switch(1+(nint>1),2,nint:1)+1]
    }
    
    for(i in 1:nf)
    {
      for(j in 1:nint)
      {
        if(shaded)
          lines(rep(object$newdata[i,xvar],2), c(object$lower[i,idx[j]],object$upper[i,idx[j]]), col=shadecols[j],lwd=6)
        else
          lines(rep(object$newdata[i,xvar],2), c(object$lower[i,idx[j]],object$upper[i,idx[j]]), col=pi.col, lty=pi.lty)
      }
    }
  }
  points(object$newdata[,xvar],object$mean,pch=19,col=fcol)
}

plot.forecast <- function(x, include, plot.conf=TRUE, shaded=TRUE, shadebars=(length(x$mean)<5),
                          shadecols=NULL, col=1, fcol=4, pi.col=1, pi.lty=2, ylim=NULL, main=NULL, ylab="",
                          xlab="", type="l",  flty = 1, flwd = 2, ...)
{
  if(is.element("x",names(x))) # Assume stored as x
    xx <- x$x
  else
    xx=NULL
  if(is.null(x$lower) | is.null(x$upper) | is.null(x$level))
    plot.conf=FALSE
  if(!shaded)
    shadebars <- FALSE
  if(is.null(main))
    main <- paste("Forecasts from ",x$method,sep="")
  if(plot.conf)
  {
    x$upper <- as.matrix(x$upper)
    x$lower <- as.matrix(x$lower)
  }
  
  if(is.element("lm",class(x$model)) & !is.element("ts",class(x$mean))) # Non time series linear model
  {
    plotlmforecast(x, plot.conf=plot.conf, shaded=shaded, shadecols=shadecols, col=col, fcol=fcol, pi.col=pi.col, pi.lty=pi.lty,
                   ylim=ylim, main=main, xlab=xlab, ylab=ylab, ...)
    if(plot.conf)
      return(invisible(list(mean=x$mean,lower=as.matrix(x$lower),upper=as.matrix(x$upper))))
    else
      return(invisible(list(mean=x$mean)))
  }
  
  # Otherwise assume x is from a time series forecast
  n <- length(xx)
  if(n==0)
    include <- 0
  else if(missing(include))
    include <- length(xx)
  
  # Check if all historical values are missing
  if(n > 0)
  {
    if(sum(is.na(xx))==length(xx))
      n <- 0
  }
  if(n > 0)
  {
    xx <- as.ts(xx)
    freq <- frequency(xx)
    strt <- start(xx)
    nx <- max(which(!is.na(xx)))
    xxx <- xx[1:nx]
    include <- min(include,nx)
  }
  else
  {
    freq <- frequency(x$mean)
    strt <- start(x$mean)
    nx <- include <- 1
    xx <- xxx <- ts(NA,frequency=freq,end=tsp(x$mean)[1]-1/freq)
  }
  pred.mean <- x$mean
  
  if(is.null(ylim))
  {
    ylim <- range(c(xx[(n-include+1):n],pred.mean),na.rm=TRUE)
    if(plot.conf)
      ylim <- range(ylim,x$lower,x$upper,na.rm=TRUE)
  }
  npred <- length(pred.mean)
  tsx <- is.ts(pred.mean)
  if(!tsx)
  {
    pred.mean <- ts(pred.mean,start=nx+1,frequency=1)
    type <- "p"
  }
  plot(ts(c(xxx[(nx-include+1):nx], rep(NA, npred)), end=tsp(xx)[2] + (nx-n)/freq + npred/freq, frequency=freq),
       xlab=xlab,ylim=ylim,ylab=ylab,main=main,col=col,type=type, ...)
  if(plot.conf)
  {
    xxx <- tsp(pred.mean)[1] - 1/freq + (1:npred)/freq
    idx <- rev(order(x$level))
    nint <- length(x$level)
    if(is.null(shadecols))
    {
      #require(colorspace)
      if(min(x$level) < 50) # Using very small confidence levels.
        shadecols <- rev(colorspace::sequential_hcl(100)[x$level])
      else # This should happen almost all the time. Colors mapped to levels.
        shadecols <- rev(colorspace::sequential_hcl(52)[x$level-49])
    }
    if(length(shadecols)==1)
    {
      if(shadecols=="oldstyle") # Default behaviour up to v3.25.
        shadecols <- heat.colors(nint+2)[switch(1+(nint>1),2,nint:1)+1]
    }
    for(i in 1:nint)
    {
      if(shadebars)
      {
        for(j in 1:npred)
        {
          polygon(xxx[j] + c(-0.5,0.5,0.5,-0.5)/freq, c(rep(x$lower[j,idx[i]],2),rep(x$upper[j,idx[i]],2)),
                  col=shadecols[i], border=FALSE)
        }
      }
      else if(shaded)
      {
        polygon(c(xxx,rev(xxx)), c(x$lower[,idx[i]],rev(x$upper[,idx[i]])),
                col=shadecols[i], border=FALSE)
      }
      else if(npred == 1)
      {
        lines(xxx+c(-0.5,0.5)/freq,rep(x$lower[,idx[i]],2),col=pi.col,lty=pi.lty)
        lines(xxx+c(-0.5,0.5)/freq,rep(x$upper[,idx[i]],2),col=pi.col,lty=pi.lty)
      }
      else
      {
        lines(xxx,x$lower[,idx[i]],col=pi.col,lty=pi.lty)
        lines(xxx,x$upper[,idx[i]],col=pi.col,lty=pi.lty)
      }
    }
  }
  if(npred > 1 & !shadebars & tsx)
    lines(pred.mean, lty = flty, lwd=flwd, col = fcol)
  else
    points(pred.mean, col=fcol, pch=19)
  if(plot.conf)
    invisible(list(mean=pred.mean,lower=x$lower,upper=x$upper))
  else
    invisible(list(mean=pred.mean))
}

predict.default <- function(object, ...)
{
  forecast(object, ...)
}

# The following function is for when users don't realise they already have the forecasts. 
# e.g., with the dshw(), meanf() or rwf() functions.

forecast.forecast <- function(object, ...)
{
  return(object)
}
ses <- function (x, h = 10, level = c(80, 95), fan = FALSE, initial=c("optimal","simple"),
                 alpha=NULL, ...)
{
  initial <- match.arg(initial)
  
  if(initial=="optimal")
    fcast <- forecast(ets(x, "ANN", alpha=alpha, opt.crit="mse"), h, level = level, fan = fan, ...)
  else
    fcast <- forecast(HoltWintersZZ(x, alpha=alpha, beta=FALSE, gamma=FALSE), h, level = level, fan = fan, ...)
  
  fcast$method <- "Simple exponential smoothing"
  fcast$model$call <- match.call()
  return(fcast)
}
accuracy <- function(f,x,test=NULL,d=NULL,D=NULL)
{
  if(is.element("mforecast", class(f)))
    return(accuracy.mforecast(f,x,test,d,D))
  
  trainset <- (is.list(f))
  testset <- (!missing(x))
  if(testset & !is.null(test))
    trainset <- FALSE
  if(!trainset & !testset)
    stop("Unable to compute forecast accuracy measures")
  
  # Find d and D
  if(testset)
  {
    d <- as.numeric(frequency(x) == 1)
    D <- as.numeric(frequency(x) > 1)
  }
  else if(trainset)
  {
    d <- as.numeric(frequency(f$mean) == 1)
    D <- as.numeric(frequency(f$mean) > 1)
  }
  else
  {
    d <- as.numeric(frequency(f)==1)
    D <- as.numeric(frequency(f) > 1)
  }
  
  
  if(trainset)
  {
    trainout <- trainingaccuracy(f,test,d,D)
    trainnames <- names(trainout)
  }
  else
    trainnames <- NULL
  if(testset)
  {
    testout <- testaccuracy(f,x,test,d,D)
    testnames <- names(testout)
  }
  else
    testnames <- NULL
  outnames <- unique(c(trainnames,testnames))
  
  out <- matrix(NA,nrow=2,ncol=length(outnames))
  colnames(out) <- outnames
  rownames(out) <- c("Training set","Test set")
  if(trainset)
    out[1,names(trainout)] <- trainout
  if(testset)
    out[2,names(testout)] <- testout
  
  if(!testset)
    out <- out[1,,drop=FALSE]
  if(!trainset)
    out <- out[2,,drop=FALSE]
  return(out)
}
#### ------------------ These functions belongs to Forecast package ----------------####
BoxCox <- function(x,lambda)
{
  if(lambda < 0)
    x[x < 0] <- NA
  if(lambda==0)
    out <- log(x)
  else
    out <- (sign(x)*abs(x)^lambda - 1)/lambda
  if(!is.null(colnames(x)))
    colnames(out) <- colnames(x)
  return(out)  
}

# Mean forecast
meanf <- function(x,h=10,level=c(80,95),fan=FALSE, lambda=NULL)
{
  xname <- deparse(substitute(x))
  n <- length(x)
  #if(!is.ts(x))
  if(!is.null(lambda))
  {
    origx <- x
    x <- BoxCox(x,lambda)
  }
  meanx <- mean(x, na.rm=TRUE)
  fits <- rep(meanx,length(x))
  res <- x-fits
  f <- rep(meanx,h)
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nconf <- length(level)
  lower <- upper <- matrix(NA,nrow=h,ncol=nconf)
  s <- sd(x,na.rm=TRUE)
  for(i in 1:nconf)
  {
    if(n > 1)
      tfrac <- qt( 0.5 - level[i]/200, n-1)
    else
      tfrac <- -Inf
    w <- -tfrac * s*sqrt(1+1/n)
    lower[,i] <- f-w
    upper[,i] <- f+w
  }
  colnames(lower) <- colnames(upper) <- paste(level,"%",sep="")
  if(is.ts(x))
  {
    fits <- ts(fits)
    res <- ts(res)
    tsp(fits) <- tsp(res) <- tsp(x)
    freq <- frequency(x)
    f <- ts(f,start=tsp(x)[2]+1/freq,frequency=freq)
    lower <- ts(lower,start=tsp(x)[2]+1/freq,frequency=freq)
    upper <- ts(upper,start=tsp(x)[2]+1/freq,frequency=freq)
    #fits <- ts(rep(NA,n))
    #if(n > 1)
    #{
    #  for(i in 2:n)
    #    fits[i] <- mean(x[1:(i-1)],na.rm=TRUE)
    #}
    #res <- x - fits	
  }
  
  if(!is.null(lambda))
  {
    fits <- InvBoxCox(fits,lambda)
    x <- origx
    f <- InvBoxCox(f,lambda)
    lower <- InvBoxCox(lower,lambda)
    upper <- InvBoxCox(upper,lambda)
  }	
  
  junk <- list(method="Mean",level=level,x=x,xname=xname,mean=f,lower=lower,upper=upper,
               model=list(mu=f[1],mu.se=s/sqrt(length(x)),sd=s), lambda=lambda, fitted=fits, residuals=res)
  junk$model$call <- match.call()
  
  return(structure(junk,class="forecast"))
}

thetaf <- function(x,h=10,level=c(80,95),fan=FALSE)
{
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  fcast <- ses(x,h=h)
  tmp2 <- lsfit(0:(length(x)-1),x)$coef[2]/2
  alpha <- fcast$model$par["alpha"]
  n <- length(x)
  fcast$mean <- fcast$mean + tmp2*(0:(h-1) + (1-(1-alpha)^n)/alpha)
  fcast.se <- sqrt(fcast$model$sigma) * sqrt((0:(h-1))*alpha^2+1)
  nconf <- length(level)
  fcast$lower <- fcast$upper <- matrix(NA,nrow=h,ncol=nconf)
  for(i in 1:nconf)
  {
    zt <- -qnorm( 0.5 - level[i]/200)
    fcast$lower[,i] <- fcast$mean - zt*fcast.se
    fcast$upper[,i] <- fcast$mean + zt*fcast.se
  }
  fcast$level <- level
  fcast$method <- "Theta"
  fcast$model <- list(alpha=alpha,drift=tmp2,sigma=fcast$model$sigma)
  fcast$model$call <- match.call()
  return(fcast)
}

# Random walk
rwf <- function(x,h=10,drift=FALSE,level=c(80,95),fan=FALSE,lambda=NULL)
{
  xname <- deparse(substitute(x))
  n <- length(x)
  freq <- frequency(x)
  nn <- 1:h
  if(!is.ts(x))
    x <- ts(x)
  if(!is.null(lambda))
  {
    origx <- x
    x <- BoxCox(x,lambda)
  }
  if(drift)
  {
    fit <- summary(lm(diff(x) ~ 1,na.action=na.exclude))
    b <- fit$coefficients[1,1]
    b.se <- fit$coefficients[1,2]
    s <- fit$sigma
    #    res <- ts(c(NA,residuals(fit)))
    method <- "Random walk with drift"
  }
  else
  {
    b <- b.se <- 0
    s <- sd(diff(x),na.rm=TRUE)
    #    fits <- ts(x[-n],start=tsp(x)[1]+1/freq,frequency=freq)
    #    res <- ts(c(NA,diff(x)))
    method <- "Random walk"
  }
  fits <- ts(c(NA,x[-n]) + b, start=tsp(x)[1], frequency=freq)
  res <- x - fits
  #tsp(res) <- tsp(fits) <- tsp(x)
  f <- x[n] + nn*b
  se <- sqrt((nn*s^2) + (nn*b.se)^2)
  
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nconf <- length(level)
  z <- qnorm(.5 + level/200)
  lower <- upper <- matrix(NA,nrow=h,ncol=nconf)
  for(i in 1:nconf)
  {
    lower[,i] <- f - z[i]*se
    upper[,i] <- f + z[i]*se
  }
  lower <- ts(lower,start=tsp(x)[2]+1/freq,frequency=freq)
  upper <- ts(upper,start=tsp(x)[2]+1/freq,frequency=freq)
  colnames(lower) <- colnames(upper) <- paste(level,"%",sep="")
  fcast <- ts(f,start=tsp(x)[2]+1/freq,frequency=freq)
  #fits <- x - res
  if(!is.null(lambda))
  {
    x <- origx
    fcast <- InvBoxCox(fcast,lambda)
    fits <- InvBoxCox(fits,lambda)
    upper <- InvBoxCox(upper,lambda)
    lower <- InvBoxCox(lower,lambda)
  }
  
  junk <- list(method=method,level=level,x=x,xname=xname,mean=fcast,lower=lower,upper=upper,
               model=list(drift=b,drift.se=b.se,sd=s), fitted = fits, residuals = res, lambda=lambda)
  junk$model$call <- match.call()
  
  return(structure(junk,class="forecast"))
}

BoxCox <- function(x,lambda)
{
  if(lambda < 0)
    x[x < 0] <- NA
  if(lambda==0)
    out <- log(x)
  else
    out <- (sign(x)*abs(x)^lambda - 1)/lambda
  if(!is.null(colnames(x)))
    colnames(out) <- colnames(x)
  return(out)	
}

InvBoxCox <- function(x,lambda)
{
  if(lambda < 0)
    x[x > -1/lambda] <- NA
  if(lambda==0)
    out <- exp(x)
  else
  {
    xx <- x*lambda + 1
    out <- sign(xx)*abs(xx)^(1/lambda)
  }
  if(!is.null(colnames(x)))
    colnames(out) <- colnames(x)
  return(out)	
}

forecast.StructTS <- function(object,h=ifelse(object$coef["epsilon"]>1e-10, 2*object$xtsp[3], 10),level=c(80,95),fan=FALSE,lambda=NULL,...)
{
  xname <- deparse(substitute(x))
  x <- object$data
  pred <- predict(object,n.ahead=h)
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nint <- length(level)
  lower <- matrix(NA,ncol=nint,nrow=length(pred$pred))
  upper <- lower
  for(i in 1:nint)
  {
    qq <- qnorm(0.5*(1+level[i]/100))
    lower[,i] <- pred$pred - qq*pred$se
    upper[,i] <- pred$pred + qq*pred$se
  }
  colnames(lower) = colnames(upper) = paste(level,"%",sep="")
  if(is.element("seas",names(object$coef)))
    method <- "Basic structural model"
  else if(is.element("slope",names(object$coef)))
    method <- "Local linear structural model"
  else 
    method <- "Local level structural model"
  
  fits <- x - residuals(object)
  if(!is.null(lambda))
  {
    fits <- InvBoxCox(fits,lambda)
    x <- InvBoxCox(x,lambda)
    pred$pred <- InvBoxCox(pred$pred,lambda)
    lower <- InvBoxCox(lower,lambda)
    upper <- InvBoxCox(upper,lambda)
  }
  
  
  return(structure(list(method=method,model=object,level=level,mean=pred$pred,lower=lower,upper=upper,
                        x=x,xname=xname,fitted=fits,residuals=residuals(object)),
                   class="forecast"))
}

forecast.HoltWinters <- function(object, h=ifelse(frequency(object$x)>1,2*frequency(object$x),10),
                                 level=c(80,95), fan=FALSE, lambda=NULL,...)
{
  xname <- deparse(substitute(x))
  x <- object$x
  if(!is.null(object$exponential))
    if(object$exponential)
      stop("Forecasting for exponential trend not yet implemented.")
  
  pred <- predict(object,n.ahead=h,prediction.interval=TRUE,level=level[1]/100)
  pmean <- pred[,1]
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nint <- length(level)
  upper <- lower <- matrix(NA,ncol=nint,nrow=length(pred[,1]))
  se <- (pred[,2]-pred[,3])/(2*qnorm(0.5*(1+level[1]/100)))
  for(i in 1:nint)
  {
    qq <- qnorm(0.5*(1+level[i]/100))
    lower[,i] <- pmean - qq*se
    upper[,i] <- pmean + qq*se
  }
  colnames(lower) = colnames(upper) = paste(level,"%",sep="")
  
  
  if(!is.null(lambda))
  {
    object$fitted[,1] <- InvBoxCox(object$fitted[,1],lambda)
    x <- InvBoxCox(x,lambda)
    pmean <- InvBoxCox(pmean,lambda)
    lower <- InvBoxCox(lower,lambda)
    upper <- InvBoxCox(upper,lambda)
  }	
  
  return(structure(list(method="HoltWinters", model=object, level=level,
                        mean=pmean, lower=lower, upper=upper,
                        x=x, xname=xname, fitted=object$fitted[,1], residuals=residuals(object)),
                   class="forecast"))
}


## CROSTON

croston <- function(x,h=10,alpha=0.1)
{
  if(sum(x<0) > 0)
    stop("Series should not contain negative values")
  out <- croston2(x,h,alpha)
  out$x <- x
  out$xname <- deparse(substitute(x))
  if(!is.null(out$fitted))
    out$residuals <- x-out$fitted
  out$method <- "Croston's method"
  return(structure(out,class="forecast"))
}

croston2 <- function(x,h=10,alpha=0.1,nofits=FALSE)
{
  x <- as.ts(x)
  y <- x[x>0]
  tsp.x <- tsp(x)
  freq.x <- tsp.x[3]
  start.f <- tsp.x[2] + 1/freq.x
  if(length(y)==0) # All historical values are equal to zero
  {
    fc <- ts(rep(0,h), start=start.f, frequency=freq.x)
    if(nofits)
      return(fc)
    else
      return(list(mean=fc, fitted=ts(x*0, start=tsp.x[1], frequency=freq.x)))
  }
  tt <- diff(c(0,(1:length(x))[x>0])) # Times between non-zero observations
  if(length(y)==1 & length(tt)==1) # Only one non-zero observation
  {
    y.f <- list(mean=ts(rep(y,h), start=start.f, frequency=freq.x))
    p.f <- list(mean=ts(rep(tt,h), start=start.f, frequency=freq.x))
  }
  else if(length(y)<=1 | length(tt)<=1) # length(tt)==0 but length(y)>0. How does that happen?
    return(list(mean=ts(rep(NA,h), start=start.f, frequency=freq.x)))
  else
  {
    y.f <- ses(y,alpha=alpha,initial="simple",h=h)
    p.f <- ses(tt,alpha=alpha,initial="simple",h=h)
  }
  ratio <- ts(y.f$mean/p.f$mean,start=start.f, frequency = freq.x)
  if(nofits)
    return(ratio)
  else
  {
    n <- length(x)
    junk <- x*NA
    if(n > 1)
    {
      for(i in 1:(n-1))
        junk[i+1] <- croston2(x[1:i],h=1,alpha=alpha,nofits=TRUE)
    }
    junk <- ts(junk)
    tsp(junk) <- tsp.x
    return(list(mean = ratio, fitted = junk, model=list(demand=y.f,period=p.f)))
  }
}


snaive <- function(x,h=2*frequency(x),level=c(80,95),fan=FALSE, lambda=NULL)
{
  fc <- forecast(Arima(x,seasonal=list(order=c(0,1,0),period=frequency(x)), lambda=lambda),h=h,level=level,fan=fan)
  # Remove initial fitted values and error
  m <- frequency(x)
  fc$fitted[1:m] <- NA
  fc$residuals[1:m] <- NA
  fc$method <- "Seasonal naive method"
  return(fc)
}

naive <- function(x,h=10,level=c(80,95),fan=FALSE, lambda=NULL)
{
  fc <- forecast(Arima(x,order=c(0,1,0),lambda=lambda),h,level=level,fan=fan)
  # Remove initial fitted values and error
  fc$fitted[1] <- NA
  fc$residuals[1] <- NA
  fc$method <- "Naive method"
  return(fc)
}
getResponse <- function(object,...) UseMethod("getResponse")

getResponse.default <- function(object,...){
  if(is.list(object))
    return(object$x)
  else
    return(NULL)
}

getResponse.lm <- function(object,...) {
  responsevar <- as.character(formula(object))[2]
  ans <- model.frame(object$model)[,responsevar]
  return(ans)
}

getResponse.Arima <- function(object,...) {
  if (is.element("x", names(object))) 
    x <- object$x
  else
  {
    series.name <- object$series
    if(is.null(series.name))
      stop("missing component series in Arima model")
    else
    {
      x <- try(eval.parent(parse(text = series.name)),silent=TRUE)
      if(class(x)=="try-error") # Try one level further up the chain
        x <- try(eval.parent(parse(text = series.name),2),silent=TRUE)
      if(class(x)=="try-error") # give up
        return(NULL)
    }
  }
  return(as.ts(x))
}

getResponse.fracdiff <- function(object, ...) {
  if (is.element("x", names(object))) 
    x <- object$x
  else
    x <- eval.parent(parse(text=as.character(object$call)[2]))
  if(is.null(tsp(x)))
    x <- ts(x,frequency=1,start=1)
  return(x)
}

getResponse.ar <- function(object, ...) {
  if (is.element("x", names(object))) 
    x <- object$x
  else
    x <- eval.parent(parse(text=object$series))
  if(is.null(tsp(x)))
    x <- ts(x,frequency=1,start=1)
  return(x)
}

croston2 <- function(x,h=10,alpha=0.1,nofits=FALSE)
{
  x <- as.ts(x)
  y <- x[x>0]
  tsp.x <- tsp(x)
  freq.x <- tsp.x[3]
  start.f <- tsp.x[2] + 1/freq.x
  if(length(y)==0) # All historical values are equal to zero
  {
    fc <- ts(rep(0,h), start=start.f, frequency=freq.x)
    if(nofits)
      return(fc)
    else
      return(list(mean=fc, fitted=ts(x*0, start=tsp.x[1], frequency=freq.x)))
  }
  tt <- diff(c(0,(1:length(x))[x>0])) # Times between non-zero observations
  if(length(y)==1 & length(tt)==1) # Only one non-zero observation
  {
    y.f <- list(mean=ts(rep(y,h), start=start.f, frequency=freq.x))
    p.f <- list(mean=ts(rep(tt,h), start=start.f, frequency=freq.x))
  }
  else if(length(y)<=1 | length(tt)<=1) # length(tt)==0 but length(y)>0. How does that happen?
    return(list(mean=ts(rep(NA,h), start=start.f, frequency=freq.x)))
  else
  {
    y.f <- ses(y,alpha=alpha,initial="simple",h=h)
    p.f <- ses(tt,alpha=alpha,initial="simple",h=h)
  }
  ratio <- ts(y.f$mean/p.f$mean,start=start.f, frequency = freq.x)
  if(nofits)
    return(ratio)
  else
  {
    n <- length(x)
    junk <- x*NA
    if(n > 1)
    {
      for(i in 1:(n-1))
        junk[i+1] <- croston2(x[1:i],h=1,alpha=alpha,nofits=TRUE)
    }
    junk <- ts(junk)
    tsp(junk) <- tsp.x
    return(list(mean = ratio, fitted = junk, model=list(demand=y.f,period=p.f)))
  }
}
WhichModels <- function(max.p, max.q, max.P, max.Q, maxK)
{
  total.models <- (max.p+1)*(max.q+1)*(max.P+1)*(max.Q+1)*length(0:maxK)
  x <- numeric(total.models)
  i <- 1
  
  for (x1 in 0:max.p) for (x2 in 0:max.q){
    for (x3 in 0:max.P) for (x4 in 0:max.Q){
      for (K in 0:maxK)
      {
        x[i] <- paste(x1, "f", x2, "f", x3, "f", x4, "f", K, sep="")
        i <- i+1
      }
    }
  }
  return(x)
}


UndoWhichModels <- function(n)
{
  as.numeric(unlist(strsplit(n, split="f")))
}
accuracy.mforecast <- function(object, x, test=NULL, d, D)
{
  fc <- object
  class(fc) <- "forecast"
  vnames <- names(object$mean)
  nox <- missing(x)
  for(i in 1:length(object$mean))
  {
    fc$mean <- object$mean[[i]]
    fc$x <- object$x[,i]
    fc$fitted <- object$fitted[,i]
    if(nox)
      out1 <- accuracy(fc, test=test, d, D)
    else
      out1 <- accuracy(fc, x[,i], test, d, D)
    rownames(out1) <- paste(vnames[i],rownames(out1))
    if(i==1)
      out <- out1
    else
      out <- rbind(out, out1)
  }
  return(out)
}
HoltWintersZZ  <- function (x,
                            # smoothing parameters
                            alpha    = NULL, # level
                            beta     = NULL, # trend
                            gamma    = NULL, # seasonal component
                            seasonal = c("additive", "multiplicative"),
                            exponential = FALSE, # exponential
                            phi = NULL # damp
)
{
  x <- as.ts(x)
  seasonal <- match.arg(seasonal)
  m <- frequency(x)
  lenx <- length(x)
  
  if(is.null(phi) || !is.numeric(phi))
    phi <- 1
  if(!is.null(alpha) && !is.numeric(alpha))
    stop ("cannot fit models without level ('alpha' must not be 0 or FALSE).")
  if(!all(is.null(c(alpha, beta, gamma))) &&
       any(c(alpha, beta, gamma) < 0 || c(alpha, beta, gamma) > 1))
    stop ("'alpha', 'beta' and 'gamma' must be within the unit interval.")
  if((is.null(gamma) || gamma > 0)) {
    if (seasonal == "multiplicative" && any(x <= 0))
      stop ("data must be positive for multiplicative Holt-Winters.")
  }
  
  if(m<=1)
    gamma <- FALSE
  
  ## initialise l0, b0, s0
  if(!is.null(gamma) && is.logical(gamma) && !gamma) {
    seasonal <- "none"
    l.start <- x[1L]
    s.start <- 0
    if(is.null(beta) || !is.logical(beta) || beta){
      if(!exponential)
        b.start <- x[2L] - x[1L]
      else
        b.start <- x[2L]/x[1L]
    }
  } else {
    ## seasonal Holt-Winters
    l.start <- mean(x[1:m])
    b.start <- (mean(x[m+(1:m)]) - l.start)/m
    if(seasonal=="additive")
      s.start <- x[1:m]-l.start
    else
      s.start <- x[1:m]/l.start
  }
  
  #initialise smoothing parameters
  #lower=c(rep(0.0001,3), 0.8)
  #upper=c(rep(0.9999,3),0.98)
  lower <- c(0,0,0,0)
  upper <- c(1,1,1,1)
  
  if(!is.null(beta) && is.logical(beta) && !beta)
    trendtype <- "N"
  else if(exponential)
    trendtype <- "M"
  else
    trendtype <- "A"
  
  if(seasonal=="none")
    seasontype <- "N"
  else if(seasonal=="multiplicative")
    seasontype <- "M"
  else
    seasontype <- "A"
  
  ## initialise smoothing parameter
  optim.start <- initparam(alpha = alpha, beta = beta, gamma=gamma, phi=1,
                           trendtype=trendtype, seasontype=seasontype, damped=FALSE, lower=lower, upper=upper, m=m)
  
  # if(!is.na(optim.start["alpha"]))
  # 	alpha2 <- optim.start["alpha"]
  # else
  # 	alpha2 <- alpha
  # if(!is.na(optim.start["beta"]))
  # 	beta2 <- optim.start["beta"]
  # else
  # 	beta2 <- beta
  # if(!is.na(optim.start["gamma"]))
  # 	gamma2 <- optim.start["gamma"]
  # else
  # 	gamma2 <- gamma
  
  #	if(!check.param(alpha = alpha2,beta = beta2, gamma = gamma2,phi=1,lower,upper,bounds="haha",m=m))
  #	{
  #		print(paste("alpha=", alpha2, "beta=",beta2, "gamma=",gamma2))
  #		stop("Parameters out of range")
  #	}
  
  ###################################################################################
  #optimisation: alpha, beta, gamma, if any of them is null, then optimise them
  error <- function (p, select)
  {
    if(select[1]>0)
      alpha <- p[1L]
    if(select[2]>0)
      beta <- p[1L+select[1]]
    if(select[3]>0)
      gamma <- p[1L+select[1]+select[2]]
    
    zzhw(x,lenx=lenx, alpha = alpha, beta=beta, gamma=gamma, seasonal=seasonal, m=m,
         dotrend=(!is.logical(beta) || beta), doseasonal=(!is.logical(gamma) || gamma),
         exponential=exponential, phi=phi, l.start=l.start, b.start=b.start, s.start=s.start)$SSE
  }
  select <- as.numeric(c(is.null(alpha),is.null(beta),is.null(gamma)))
  
  if(sum(select)>0) # There are parameters to optimize
  {
    sol <- optim(optim.start, error, method = "L-BFGS-B", lower = lower[select], upper = upper[select], select=select)
    if(sol$convergence || any(sol$par < 0 | sol$par > 1)) {
      if (sol$convergence > 50) {
        warning(gettextf("optimization difficulties: %s", sol$message), domain = NA)
      } else stop("optimization failure")
    }
    if(select[1]>0)
      alpha <- sol$p[1L]
    if(select[2]>0)
      beta <- sol$p[1L+select[1]]
    if(select[3]>0)
      gamma <- sol$p[1L+select[1]+select[2]]
  }
  
  final.fit <- zzhw(x, lenx=lenx, alpha = alpha, beta=beta, gamma=gamma, seasonal=seasonal, m=m,
                    dotrend=(!is.logical(beta) || beta), doseasonal=(!is.logical(gamma) || gamma),
                    exponential=exponential, phi=phi, l.start=l.start, b.start=b.start, s.start=s.start)
  
  tspx <- tsp(x)
  fitted <- ts(final.fit$fitted,frequency=m,start=tspx[1])
  states <- matrix(final.fit$level,ncol=1)
  colnames(states) <- "l"
  if(trendtype!="N")
    states <- cbind(states,b=final.fit$trend)
  if(seasontype!="N")
  {
    nr <- nrow(states)
    nc <- ncol(states)
    for(i in 1:m)
      states <- cbind(states,final.fit$season[(m-i)+(1:nr)])
    colnames(states)[nc+(1:m)] <- paste("s",1:m,sep="")
  }
  states <- ts(states,frequency=m,start=tspx[1]-1/m)
  
  # Package output as HoltWinters class
  # structure(list(fitted    = fitted,
  # 				x        = x,
  # 				alpha    = alpha,
  # 				beta     = beta,
  # 				gamma    = gamma,
  # 				coefficients = c(a = final.fit$level[lenx],
  # 						b = if (!is.logical(beta) || beta) final.fit$trend[lenx],
  # 						s = if (!is.logical(gamma) || gamma) final.fit$season[lenx - m + 1L:m]),
  # 				seasonal  = seasonal,
  # 				exponential = exponential,
  # 				SSE       = final.fit$SSE,
  # 				call      = match.call(),
  # 				level = final.fit$level,
  # 				trend = final.fit$trend,
  # 				season = final.fit$season,
  # 				phi = phi
  # 		),
  # 		class = "HoltWinters"
  # )
  # Package output as ets class
  damped <- (phi<1.0)
  if(seasonal=="additive") # This should not happen
    components <- c("A",trendtype,seasontype,damped)
  else if(seasonal=="multiplicative")
    components <- c("M",trendtype,seasontype, damped)
  else if(seasonal=="none" & exponential)
    components <- c("M",trendtype,seasontype,damped)
  else# if(seasonal=="none" & !exponential)
    components <- c("A",trendtype,seasontype, damped)
  
  initstate <- states[1,]
  param <- alpha
  names(param) <- "alpha"
  if(trendtype!="N")
  {
    param <- c(param,beta=beta)
    names(param)[length(param)] <- "beta"
  }
  if(seasontype!="N")
  {
    param <- c(param,gamma=gamma)
    names(param)[length(param)] <- "gamma"
  }
  if(damped)
  {
    param <- c(param,phi=phi)
    names(param)[length(param)] <- "phi"
  }
  
  if(components[1]=="A")
    sigma2 <- mean(final.fit$residuals^2)
  else
    sigma2 <- mean((final.fit$residuals/fitted)^2)
  structure(list(fitted    = fitted,
                 residuals=final.fit$residuals,
                 components=components,
                 x        = x,
                 par=c(param,initstate),
                 initstate=initstate,
                 states=states,
                 SSE       = final.fit$SSE,
                 sigma2 = sigma2,
                 call      = match.call(),
                 m = m
  ),
  class = "ets"
  )
}

ets <- function(y, model="ZZZ", damped=NULL,
                alpha=NULL, beta=NULL, gamma=NULL, phi=NULL, additive.only=FALSE, lambda=NULL,
                lower=c(rep(0.0001,3), 0.8), upper=c(rep(0.9999,3),0.98),
                opt.crit=c("lik","amse","mse","sigma","mae"), nmse=3, bounds=c("both","usual","admissible"),
                ic=c("aicc","aic","bic"),restrict=TRUE, use.initial.values=FALSE, ...)
{
  #dataname <- substitute(y)
  opt.crit <- match.arg(opt.crit)
  bounds <- match.arg(bounds)
  ic <- match.arg(ic)
  
  #if(max(y,na.rm=TRUE) > 1e6)
  #    warning("Very large numbers which may cause numerical problems. Try scaling the data first")
  
  if(any(class(y) %in% c("data.frame","list","matrix","mts")))
    stop("y should be a univariate time series")
  y <- as.ts(y)
  
  # Check if data is constant
  if (is.constant(y)) 
    return(ses(y, alpha=0.99999, initial='simple')$model)
  
  # Remove missing values near ends
  ny <- length(y)
  y <- na.contiguous(y)
  if(ny != length(y))
    warning("Missing values encountered. Using longest contiguous portion of time series")
  
  orig.y <- y
  if(class(model)=="ets" & is.null(lambda))
    lambda <- model$lambda
  if(!is.null(lambda))
  {
    y <- BoxCox(y,lambda)
    additive.only=TRUE
  }
  
  if(nmse < 1 | nmse > 10)
    stop("nmse out of range")
  m <- frequency(y)
  
  if(sum((upper-lower)>0)<4)
    stop("Lower limits must be less than upper limits")
  
  # If model is an ets object, re-fit model to new data
  if(class(model)=="ets")
  {
    alpha <- model$par["alpha"]
    beta <- model$par["beta"]
    if(is.na(beta))
      beta <- NULL
    gamma <- model$par["gamma"]
    if(is.na(gamma))
      gamma <- NULL
    phi <- model$par["phi"]
    if(is.na(phi))
      phi <- NULL
    modelcomponents <- paste(model$components[1],model$components[2],model$components[3],sep="")
    damped <- (model$components[4]=="TRUE")
    if(use.initial.values)
    {
      errortype  <- substr(modelcomponents,1,1)
      trendtype  <- substr(modelcomponents,2,2)
      seasontype <- substr(modelcomponents,3,3)
      
      # Recompute errors from pegelsresid.C
      e <- pegelsresid.C(y, m, model$initstate, errortype, trendtype, seasontype, damped, alpha, beta, gamma, phi)
      
      # Compute error measures
      np <- length(model$par)
      model$loglik <- -0.5*e$lik
      model$aic <- e$lik + 2*np
      model$bic <- e$lik + log(ny)*np
      model$aicc <- e$lik + 2*ny*np/(ny-np-1)
      model$mse <- e$amse[1]
      model$amse <- mean(e$amse[1:nmse])
      
      # Compute states, fitted values and residuals
      tsp.y <- tsp(y)
      model$states=ts(e$states,frequency=tsp.y[3],start=tsp.y[1]-1/tsp.y[3])
      colnames(model$states)[1] <- "l"
      if(trendtype!="N")
        colnames(model$states)[2] <- "b"
      if(seasontype!="N")
        colnames(model$states)[(2+(trendtype!="N")):ncol(model$states)] <- paste("s",1:m,sep="")
      if(errortype=="A")
        model$fitted <- ts(y-e$e,frequency=tsp.y[3],start=tsp.y[1])
      else
        model$fitted <- ts(y/(1+e$e),frequency=tsp.y[3],start=tsp.y[1])
      model$residuals <- ts(e$e,frequency=tsp.y[3],start=tsp.y[1])
      model$sigma2 <- mean(model$residuals^2,na.rm=TRUE)
      model$x <- orig.y
      if(!is.null(lambda))
      {
        model$fitted <- InvBoxCox(model$fitted,lambda)
      }
      model$lambda <- lambda
      
      # Return model object
      return(model)
    }
    else
    {
      model <- modelcomponents
    }
  }
  
  errortype  <- substr(model,1,1)
  trendtype  <- substr(model,2,2)
  seasontype <- substr(model,3,3)
  
  if(!is.element(errortype,c("M","A","Z")))
    stop("Invalid error type")
  if(!is.element(trendtype,c("N","A","M","Z")))
    stop("Invalid trend type")
  if(!is.element(seasontype,c("N","A","M","Z")))
    stop("Invalid season type")
  
  if(m < 1 | length(y) <= m)
  {
    #warning("I can't handle data with frequency less than 1. Seasonality will be ignored.")
    seasontype <- "N"
  }
  if(m == 1)
  {
    if(seasontype=="A" | seasontype=="M")
      stop("Nonseasonal data")
    else
      substr(model,3,3) <- seasontype <- "N"
  }
  if(m > 24)
  {
    if(is.element(seasontype,c("A","M")))
      stop("Frequency too high")
    else if(seasontype=="Z")
    {
      warning("I can't handle data with frequency greater than 24. Seasonality will be ignored. Try stlf() if you need seasonal forecasts.")
      substr(model,3,3) <- seasontype <- "N"
      #m <- 1
    }
  }
  
  # Check inputs
  if(restrict)
  {
    if((errortype=="A" & (trendtype=="M" | seasontype=="M")) |
         (errortype=="M" & trendtype=="M" & seasontype=="A") |
         (additive.only & (errortype=="M" | trendtype=="M" | seasontype=="M")))
      stop("Forbidden model combination")
  }
  
  data.positive <- (min(y) > 0)
  
  if(!data.positive & errortype=="M")
    stop("Inappropriate model for data with negative or zero values")
  
  if(!is.null(damped))
  {
    if(damped & trendtype=="N")
      stop("Forbidden model combination")
  }
  
  # Check we have enough data to fit a model
  n <- length(y)
  npars <- 2L # alpha + l0
  if(trendtype=="A" | trendtype=="M") 
    npars <- npars + 2L # beta + b0
  if(seasontype=="A" | seasontype=="M")
    npars <- npars + m # gamma + s
  if(!is.null(damped))
    npars <- npars + as.numeric(damped)
  if(n <= npars + 1)
    stop("You've got to be joking. I need more data!")
  
  # Fit model (assuming only one nonseasonal model)
  if(errortype=="Z")
    errortype <- c("A","M")
  if(trendtype=="Z")
    trendtype <- c("N","A","M")
  if(seasontype=="Z")
    seasontype <- c("N","A","M")
  if(is.null(damped))
    damped <- c(TRUE,FALSE)
  best.ic <- Inf
  for(i in 1:length(errortype))
  {
    for(j in 1:length(trendtype))
    {
      for(k in 1:length(seasontype))
      {
        for(l in 1:length(damped))
        {
          if(trendtype[j]=="N" & damped[l])
            next
          if(restrict)
          {
            if(errortype[i]=="A" & (trendtype[j]=="M" | seasontype[k]=="M"))
              next
            if(errortype[i]=="M" & trendtype[j]=="M" & seasontype[k]=="A")
              next
            if(additive.only & (errortype[i]=="M" | trendtype[j]=="M" | seasontype[k]=="M"))
              next
          }
          if(!data.positive & errortype[i]=="M")
            next
          fit <- etsmodel(y,errortype[i],trendtype[j],seasontype[k],damped[l],alpha,beta,gamma,phi,
                          lower=lower,upper=upper,opt.crit=opt.crit,nmse=nmse,bounds=bounds, ...)
          fit.ic <- switch(ic,aic=fit$aic,bic=fit$bic,aicc=fit$aicc)
          if(!is.na(fit.ic))
          {
            if(fit.ic < best.ic)
            {
              model <- fit
              best.ic <- fit.ic
              best.e <- errortype[i]
              best.t <- trendtype[j]
              best.s <- seasontype[k]
              best.d <- damped[l]
            }
          }
        }
      }
    }
  }
  if(best.ic == Inf)
    stop("No model able to be fitted")
  
  model$m <- m
  model$method <- paste("ETS(",best.e,",",best.t,ifelse(best.d,"d",""),",",best.s,")",sep="")
  model$components <- c(best.e,best.t,best.s,best.d)
  model$call <- match.call()
  model$initstate <- model$states[1,]
  model$sigma2 <- mean(model$residuals^2,na.rm=TRUE)
  model$x <- orig.y
  model$lambda <- lambda
  if(!is.null(lambda))
  {
    model$fitted <- InvBoxCox(model$fitted,lambda)
  }
  
  #model$call$data <- dataname
  
  return(structure(model,class="ets"))
}
