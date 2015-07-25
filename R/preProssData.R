`preProssData`<-function(y,x,dif=FALSE,remove.outliers=TRUE,remove.na=TRUE,dfn=1,dfs=0,ses.l=1,plots=TRUE,outlierCoef=1.5,kpss.Test=TRUE) {
  org.y=y
  org.x=x
  if(!require('tseries')) stop('It needs tseries package to do KPSS.test',call. = 0)
  #----- REMOVE NA ------------
  if (remove.na==TRUE){
    tmp=na.omit(cbind(y,x))
    y=tmp[,1]
    x=tmp[,-1]
  }
  #--- DIFFERENTIATING -----
  if(dif==TRUE){
    if(dfn>0){
      if (dfs>0){
        y=diff(x = y,lag = ses.l,differences = dfs)
        y=diff(x = y,lag = 1    ,differences = dfn)
        
        x=apply(x,2,difff1<-function(x,l=ses.l,s=dfs)(diff(x,lag = l,differences = s)))
        x=apply(x,2,difff2<-function(x,l=1,s=dfn)    (diff(x,lag = l,differences = s)))
        
        
      }else{
        y=diff(y,lag = 1,differences = dfn)
        x=apply(x,2,difff3<-function(x,d=dfn)(diff(x,lag = 1,differences = d)))
      }
    }else{
      if (dfs>0){
        y=diff(y,lag = ses.l,differences = dfs)
        x=apply(x,2,difff4<-function(x,s=dfs,l=ses.l){(diff(x,lag = l,differences = s))})
      }
    }
  }
  #---- REMOVEING OUTLIERS ---------
  totalm=cbind(y,x) 
  if(remove.outliers==TRUE){
    o=c()
    qnt <- quantile(y, probs=c(.25, .75), na.rm = TRUE)  
    H <- outlierCoef * IQR(y, na.rm = TRUE)
    t <- y
    otl=which(y < (qnt[1] - H)) 
    otg=which(y > (qnt[2] + H)) 
    
    if(length(otl)>0){
      if(length(otg)>0){
        o=c(otl,otg)
        
      }else{
        o=c(otl)
      }    
    }else{
      if(length(otg)>0){
        o=c(otg)      
      }
    }  
    if (length(o)>0){
      totalm=totalm[-c(o),]
      cat(length(o),' outlier(s) removed [index:',o,']\n',sep=' ')
    }
  }
  if (plots==TRUE){
    par(mfrow=c(2,1))
    plot(org.y,type='b',ylab='y before transformations')
    plot(totalm[,1],type='b',ylab='y after transformations')
    par(mfrow=c(1,1))
  }
  
  #----- KPSS test ----#
  if(kpss.Test == TRUE){
    cat('\n\n KPSS.TEST in progress ... \n')
    kprL=apply(as.matrix(totalm),2,f3<-function(x){kpss.test(x,null='Level',lshort = 1)$p.value})
    kprT=apply(as.matrix(totalm),2,f4<-function(x){kpss.test(x,null='Trend',lshort = 1)$p.value})
    print(kprL)
    print(kprT)
    #mat=rbind(kprL,kprT)
    #rownames(mat)=c('Level','Trend')
    #cat('\n',round(mat,3),'\n\n')
  }
  return(totalm)
}

vectorOutliers<-function(x,percent=FALSE){
  if(length(x)==0 || is.null(x)==TRUE){
    stop ('X must be a vector')
  }
  if(percent==TRUE){
    n=length(x)
  }
  o=c()
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)  
  H <- 1.5 * IQR(x, na.rm = TRUE)
  otl=which(x < (qnt[1] - H)) 
  otg=which(x > (qnt[2] + H)) 
  o=c(otl,otg)
  if (percent==FALSE){n=1}
  lo=length(o)/n
  if(length(o)>=0){
    return(lo)
  }else{
    return(0)
  }
}

plotOfoutliers<-function(x){
  if (is.matrix(x)==FALSE){
    stop ('x must be a matrix');
  } 
  outp=apply(x,2,vectorOutliers)
  barplot(outp[which(outp>0)],ylab='number of outliers',las = 2, cex.names=.7)
  return(outp)
}

HuberScale<-function(x){
  f1<-function(x){MASS::huber(x)$mu}
  f2<-function(x){MASS::huber(x)$s}
  meanx=apply(x,2,f1)
  normx=apply(x,2,f2)
  xs=scale(x,meanx,normx) 
  return(xs)
}

OOacf<-function(y,main='ACF',R=1, ...){
  a=acf(y,plot=FALSE,...)$acf[-1]
  ci=qnorm(.975)/sqrt(length(y/R))
  aMax=max(abs(a))
  #print(R)
  if(aMax<ci){
    plot(1:length(a),a,main=main,xlab='Lag',type='h',ylab='ACF',ylim=c(-ci*1.5,ci*1.5))
  }else{
    plot(1:length(a),a,main=main,xlab='Lag',type='h',ylab='ACF',ylim=c(-aMax*1.2,aMax*1.2))
  }
  abline(h=c(-ci,0,ci),lty=2)
}

OOpacf<-function(y,main='ACF',R=1,...){
  a=pacf(y,plot=FALSE,...)$acf
  ci=qnorm(.975)/sqrt(length(y/R))
  #print(R)
  aMax=max(abs(a))
  if(aMax<ci){
    plot(1:length(a),a,main=main,xlab='Lag',type='h',ylab='PACF',ylim=c(-ci*1.5,ci*1.5))
  }else{
    plot(1:length(a),a,main=main,xlab='Lag',type='h',ylab='PACF',ylim=c(-aMax*1.2,aMax*1.2))
  }
  abline(h=c(-ci,0,ci),lty=2)
}

#REGRESSION WITH MOORE-PENROSE inverse
mplm<-function(x=1,y=1){
  x=as.matrix(x)
  y=as.matrix(y)
  beta= pinv(t(x) %*% x) %*% t(x) %*% y 
  error=y - x %*% beta
  return(list
         (
         #--------- Data ---------#
         y=y, x=x,
         coef=beta,
         error=error
         )
  )
}
mplm2<-function(x=1,y=1){
  x=as.matrix(x)
  y=as.matrix(y)
  beta= pinv(t(x) %*% x) %*% t(x) %*% y 
  error=y - x %*% beta
  return(list
         (
         #--------- Data ---------#
         y=y, x=x,
         coefficients=beta,
         error=error
         )
  )
}
MinBox.test<-function(y,lag=c(3,30),df=1){
  if(min(lag)>max(lag)) stop('Lag are not in proper orders')
  if(min(lag)<1) stop('Lag zero in not allowded')
  r=c()
  for (i in lag[1]:lag[2]){
    r[i]=Box.test(y,lag=i,fitdf=df,type='Ljung-Box')$p.value
  }
  return(min(r,na.rm=TRUE))
}

RRC<-function(x,nr=0,nc=0,first=TRUE){
  if (is.matrix(x)==FALSE){
    stop('X is not matrix!')
  }
  xnc=ncol(x)
  xnr=nrow(x)
  if (first==TRUE){
    if(nr>0){
      x=x[-(1:nr),]
    }
    if(nc>0){
      x=x[,-(1:nc)]
    }
  }else{
    if(nr>0){
      x=x[-((xnr-nr):xnr),]
    }
    if(nc>0){
      x=x[,-((xnc-nc):xnc)]
    }
  }
  return(x)  
}

Nprint<-function(x,display=FALSE){
  if (display==TRUE){
    print(x)
  }
}
Ncat<-function(x,display=FALSE){
  if (display==TRUE){
    cat(x)
  }
}


CBIND<-function(x1,x2=NULL,x3=NULL,n=0){
  if (is.matrix(x1)==FALSE){stop('X1 is not a matrix!')}
  cx1=nrow(x1)
  if (is.null(x2)==FALSE){cx2=nrow(x2)}else{cx2=Inf}
  if (is.null(x3)==FALSE){cx3=nrow(x3)}else{cx3=Inf}
  c=c(cx1,cx2,cx3);
  rmin=min(c,na.rm=TRUE)
  minIndex=which(c==rmin)
  if(minIndex!=1){
    x=x1[(cx1-rmin+1):cx1,]
  }else{
    x=x1
  }
  if(is.finite(cx2)==TRUE){
    if(minIndex!=2){
      Newx2=x2[(cx2-rmin+1):cx2,]
      x=cbind(x,Newx2)
    }else{
      x=cbind(x,x2)
    }
  }
  if(is.finite(cx3)==TRUE){
    if(minIndex!=3){
      Newx3=x3[(cx3-rmin+1):cx3,]
      x=cbind(x,Newx3)
    }else{
      x=cbind(x,x3)
    }
  }
  if(n>0){
    x=x[-(1:n),]
  }
  return(x)
}

MAE<-function(x1=NULL,x2=NULL){
  x1=as.vector(x1)
  x2=as.vector(x2)
  if(length(x1)!=length(x2) || is.null(x1)==TRUE || is.null(x2)==TRUE){stop('two vector must be of the same length!')}
  return(sum(abs(x1-x2))/length(x1))
}

MSE<-function(x1=NULL,x2=NULL,v.root=FALSE){
  x1=as.vector(x1)
  x2=as.vector(x2)
  if(length(x1)!=length(x2) || is.null(x1)==TRUE || is.null(x2)==TRUE){stop('two vector must be of the same length!')}
  if(v.root==TRUE){
    s=sqrt(sum(abs(x1-x2)^2)/length(x1))
  }else{
    s=(sum(abs(x1-x2)^2)/length(x1))
  }
  return(s)
}


Make.LagMatrix<-function(nc=10){
  m=matrix(0,ncol=nc,nrow=nc)
  nzl=1
  for (i in 2:nc){
    m[i,nzl]=1
    nzl=nzl+1   
  }
  return(m)
}

#--- Discrete uniform random number
rdunif<-function(n=1,a=0,b=1){
  x=runif(n,a,b)
  return(round(x))
}

# Computes power of a matrix   #
"%^%" <- function(S, power) {
  with(eigen(S), vectors %*% (values^power * t(vectors))) 
}


# ERF function
erf<-function(x,s=1){
  r=(2*pnorm(x/s,0,1/sqrt(2))-1)
  return(r)
}
# --- ERF with high accuracy ----
erf2<-function(x){
  if(x<4) stop('This method does not work for this value!')
  j=0:50
  if(require('Rmpfr')==TRUE){
    e=mpfr(exp(1),100)
    pi2=mpfr(pi,100)
    s=mpfr((-1)^j * factorial(2*j)/ factorial(j) * (2*x)^(-2*j),100)
  }else{
    e=exp(1)
    pi2=pi
    s=(-1)^j * factorial(2*j)/ factorial(j) * (2*x)^(-2*j)
  }
  
  k=which(abs(s)==min(abs(s)))
  y=e^(-x^2)/(x*sqrt(pi2)) * sum(s)
  return(1-y)
}

# ERF approximation using Kummer
erfAp<-function(z=0,a=1,b=1.5,maxt=5){
  z=z^2
  s=1:maxt
  a=a+s
  b=b+s
  ss=1
  for(i in s){
    mt=prod(a[1:i]/b[1:i])
    ss = ss + mt * (z^i) / factorial(i)
  }
  ss=2*sqrt(z)/sqrt(pi)*exp(-z)*ss
  return(ss)
}
# Another ERF approximation
erfa3<-function(x=1,s=1){
  x=x/s
  e=exp(1)
  y=tanh(39*x/sqrt(4*pi)-111/2*atan(35*x/111/sqrt(pi)))
  return(y)
}

# ---- PRODUCE Hp and Hq ----#
`makematrix`<-function (y,p=3,beforestart.mean=mean(y),rname='H'){
  if (p>0 && is.null(p)==FALSE && length(p)>0){
    ny=length(y)
    
    H=matrix(ncol=p,nrow=ny)
    progress=1
    
    for (i in 1:ny){
      if(progress <= p){
        for (j in progress:1){
          H[i,progress-j]=y[j]  #+1
        }
      }else {
        col=1
        
        for (j in progress:(progress - p +1)){ 
          H[i,col]=y[j-1] #-1
          col=col+1
        }
      }
      progress=progress+1
    }
    H[is.na(H)]=beforestart.mean
    colnames(H)=paste(rname,1:p,sep='.')
    return(H)
  }else{
    return (0)
  }
}


#---- MOVING AVERAGE FUNCTION -----#
movingAverage<-function(x,order=1,Draw.plot=FALSE){
  n=length(x)
  mav=c()
  if(order>0 && n>0){
    for(i in order:(n)){
      mav[i-order+1]=mean(x[(i-order+1):(i)])
      
    }
  }
  if(length(mav)>0){
    if(Draw.plot==TRUE){
      plot(mav,type='b')
    }
    return (mav)
  }else{
    return('ERROR IN CALCULATIONS!')
  }
}

MmovingAverage<-function(x,ord=1){
  if (is.matrix(x)==FALSE){
    stop('x must be matrix')    
  }
  mov.avr<-function(x,order=ord){
    mov.avr =movingAverage(x,order,Draw.plot=FALSE) 
    return(mov.avr)
  }
  return(apply(x,2,mov.avr))
}

# simulate regression model
sim.regression<-function(n.obs=10,coefficients=runif(10,-5,5),s.deviation=.1,normalize=TRUE,n.zero=0,Coliniarity=0,shuffle=FALSE){
  
  n.var=length(coefficients) 
  
  M=matrix(rnorm(n.var*n.obs),ncol=n.var,nrow=n.obs) 
  
  beta=as.matrix(coefficients)
  
  if(round(n.zero)>0){
    
    cp=round(n.zero)
    
    M2=matrix(rnorm(cp*n.obs,0,.1),ncol=cp,nrow=n.obs)
    
    M=cbind(M,M2)
    
    beta=as.matrix(c(coefficients,rep(0,cp)))
    
  }
  
  nco=round(Coliniarity*n.var)
  
  if(Coliniarity>0 & nco>0 & nco<=n.var){
    
    mM=ncol(M);nM=nrow(M)
    
    for(k in 1:nco){
      
      SaC=sample(1:mM, 2, replace = FALSE)
      
      SaCC=rowSums(M[,SaC])+rnorm(length(SaC),0,.5)
      
      M=cbind(M,SaCC)
      
      beta=rbind(beta,1)
      
    }
    
  }
  
  
  if(normalize==TRUE){
    y=scale(M)  %*% beta +rnorm(n.obs,0,runif(1,0,s.deviation))
  }else{
    y=M  %*% beta +rnorm(n.obs,0,runif(1,0,s.deviation))
  }
  if(shuffle==TRUE){
    col.sa=sample(1:ncol(M),replace = FALSE)
    M=M[,col.sa]
    beta=as.matrix(beta[col.sa])
  }
  return (list(x=M,y=y,coeff=beta))
  
}
#set.seed(1);s=sim.regression(n.obs = 5,coefficients = 1,n.zero = 3);dim(s$x)
#--- Moore penros pseudoinverse
pinv<-function(H){
  x=t(H) %*% H
  s=svd(x)
  xp=s$d
  for (i in 1:length(xp)){
    if (xp[i] != 0){
      xp[i]=1/xp[i]
    }
    else{
      xp[i]=0
    }
  }
  return(s$u %*% diag(xp) %*% t(s$v) %*% t(H))
}
#------ SDD ------------
sdd<-function (A, kmax = 100, alphamin = 0.01, lmax = 100, rhomin = 1e-19) 
{
  betabar = 1
  yinit = 1
  idx = 1
  m = nrow(A)
  n = ncol(A)
  Xsav = matrix(0, nrow = m, ncol = 0)
  Dsav = matrix(0, nrow = 0, ncol = 1)
  Ysav = matrix(0, nrow = n, ncol = 0)
  rhosav = numeric(0)
  rho = sum(A^2)
  iitssav = numeric(0)
  for (k in 1:kmax) {
    s = matrix(0, nrow = m, ncol = 1)
    iits = 0
    while ((svd(s)$d[1])^2 < (rho/n)) {
      y = matrix(0, nrow = n, ncol = 1)
      y[idx] = 1
      s = A %*% y
      if (k > 1) {
        s = s - (Xsav %*% (Dsav * (t(Ysav) %*% y)))
      }
      idx = idx%%n + 1
      iits = iits + 1
    }
    iitssav[k] = iits
    for (l in 1:lmax) {
      s = A %*% y
      if (k > 1) {
        s = s - (Xsav %*% (Dsav * (t(Ysav) %*% y)))
      }
      tmp = sddsolve(s, m)
      x = tmp[[1]]
      xcnt = tmp[[2]]
      s = t(A) %*% x
      if (k > 1) {
        s = s - (Ysav %*% (Dsav * (t(Xsav) %*% x)))
      }
      tmp2 = sddsolve(s, n)
      y = tmp2$x
      ycnt = tmp2$imax
      fmax = tmp2$fmax
      d = sqrt(fmax %*% ycnt)/(ycnt %*% xcnt)
      beta = d^2 %*% ycnt %*% xcnt
      if (l > 1) {
        alpha = (beta - betabar)/betabar
        if (alpha <= alphamin || abs(beta - betabar) < 
            1e-19) {
          break
        }
      }
      betabar = beta
    }
    Xsav = cbind(Xsav, x)
    Ysav = cbind(Ysav, y)
    Dsav = rbind(Dsav, d)
    rho = pmax(rho - beta, 0)
    rhosav[k] = rho
    if (rho <= rhomin) {
      break
    }
    if (k >= 2 && (rhosav[k] - rhosav[k - 1] == 0)) {
      break
    }
  }
  ret = list()
  colnames(Xsav) = rep("", ncol(Xsav))
  ret$x = Xsav
  ret$d = as.numeric(Dsav)
  colnames(Ysav) = rep("", ncol(Ysav))
  ret$y = Ysav
  return(ret)
}

### overlay histogram ######
plotOverlappingHist <- function(a, b, colors=c("white","gray20","gray50"),
                                breaks=NULL, xlim=NULL, ylim=NULL){
  
  ahist=NULL
  bhist=NULL
  
  if(!(is.null(breaks))){
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  } else {
    ahist=hist(a,plot=F)
    bhist=hist(b,plot=F)
    
    dist = ahist$breaks[2]-ahist$breaks[1]
    breaks = seq(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks),dist)
    
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  }
  
  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }
  
  if(is.null(ylim)){
    ylim = c(0,max(ahist$counts,bhist$counts))
  }
  
  overlap = ahist
  for(i in 1:length(overlap$counts)){
    if(ahist$counts[i] > 0 & bhist$counts[i] > 0){
      overlap$counts[i] = min(ahist$counts[i],bhist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }
  
  plot(ahist, xlim=xlim, ylim=ylim, col=colors[1])
  plot(bhist, xlim=xlim, ylim=ylim, col=colors[2], add=T)
  plot(overlap, xlim=xlim, ylim=ylim, col=colors[3], add=T)
}


#----- Permutation -------
perm = function(n, x) {
  return(factorial(n) / factorial(n-x))
}
#----- Combination -------
comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}

#----- Make a matrix symmetrix -----
makeSymmetric<-function(x){
  x=as.matrix(x)
  o1=diag(rep(0,nrow(x)))
  o2=diag(rep(0,ncol(x)))
  x=rbind(cbind(o1,x),cbind(t(x),o2))
  return(x)
}

#------ Update all packages in R directory ------
UpdateR<-function(){
  lib <- .libPaths()[1]
  install.packages( 
    lib  = lib,
    pkgs = as.data.frame(installed.packages(lib), stringsAsFactors=FALSE)$Package,
    type = 'source'
  )
}

erfAprSin<-function(x){
  if(abs(x)<= 1.5){
    erf=1/sqrt(4*pi) * (sin(pi*x/10)+sin(x))+.5
  }
  if(x > 1.5){
    erf=1-exp(-1.775*x) + x/exp(x+4)
  }
  if(x< (-1.5)){
    x=abs(x)
    erf=exp(-1.775*x) - x/exp(x+4)
  }
  return(erf)
}


#---- THIS IS DIFFERENT FROM MAKEMATRIX! (JUST A TRIANGULAR MATRIX WITH DATA IN FIRST COLUMN) ------#
MakeTriangularMatrix<-function(v,p=1,justOne=FALSE){
  x=as.matrix(v)
  outp=x
  nv=length(v)
  if(p>nv){p=nv}  
  for(i in 1:p){
    if(i<nv){
      s=c(rep(0,i),x[(i+1):nv])
    }else if(i==nv){      
      s=c(rep(0,i))
    }
    outp=cbind(outp,s)  
  }
  colnames(outp)=paste('X.',1:(p+1))
  if(justOne){outp[outp!=0]=1}
  return(outp)
}
#----------- MAKE LAG MATRIX -------------
makeLagmatrix<-function(y){
  ny=length(y)
  x=matrix(0,ncol=ny,nrow=ny)
  for(i in 2:ny){
    x[i,i-1]=1
  }
  return(x)
}

# find mode of a distribution
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
