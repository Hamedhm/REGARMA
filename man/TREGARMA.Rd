
\name{regARMA}
\alias{regARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
    A function for penalised estimation of a REGARMA model from time-series regression data.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    This function estimates the parameters of a REGARMA model using lasso and adaptive lasso penalties for a given dataset with or without replications. It also selects the optimal model using AIC, BIC and GCV.
}
\usage{
  regARMA(data, ar = 0, ma = 0, method = c("alasso"), mselection = "BIC", 
  alpha = 0, Ndf = 0, Sdf = 0, normalize = FALSE, auto.order = FALSE, 
  with.rep = NA, non.penalized = FALSE, iteration=10, auto.prune = TRUE, 
  debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    %%     ~~Describe \code{data} here~~
      Raw data of the form of a matrix. The first column should contain ids of replications. Note that the points belonging to each  replication id must be sorted with respect to time. The second column should contain the response variable (y), the remaining columns should contain the predictors (x). In the case of no replication, the first column should be simply 1:T, with T the number of time points.
    
    
  }
  \item{ar}{
    %%     ~~Describe \code{ar} here~~
      AR order of model
  }
  \item{ma}{
    %%     ~~Describe \code{ma} here~~
      MA order of model 
  }
  \item{method}{
    %%     ~~Describe \code{method} here~~
      Type of penalties. Choose between elastic net (enet) and adaptive lasso (alasso). Note that lasso can be computed by choosing enet for method and setting alpha equals to zero.
  }
  \item{mselection}{
    %%     ~~Describe \code{mselection} here~~
      Model selection criteria. Choose between AIC, BIC and GCV for the selection of the tuning parameters. %%(See also msgps function in msgps package for more details) 
  }
  \item{alpha}{
    %%     ~~Describe \code{alpha} here~~
      Weight in elastic net.
  }
  \item{Ndf}{
    %%     ~~Describe \code{Ndf} here~~
      Degree of differencing if necessary
  }
  \item{Sdf}{
    %%     ~~Describe \code{Sdf} here~~
      Degree of seasonal differencing if necessary
  }
  \item{normalize}{
    %%     ~~Describe \code{normalize} here~~
      Set to TRUE in order to normalize the data to mean zero and unit variance 
  }
  \item{auto.order}{
    %%     ~~Describe \code{auto.order} here~~
      If TRUE, the optimal orders of the model are chosen according to AIC (if mselection="AIC") or BIC (if mselection="BIC" or "GCV"). %( See also auto.arima function in forecast package for more detail)
  }
  \item{with.rep}{
    %%     ~~Describe \code{auto.order} here~~
      If set to 1 then model with replications is applied. If set to 0 then model without replication is applies. If set to NA then function decides about the model. %( See also auto.arima function in forecast package for more detail)
  }
  \item{non.penalized}{
    %%     ~~Describe \code{debug} here~~
      If set TRUE then non-penalized model is applied. \cr
    Notice: in this case method, mselection, alpha and with.rep are not working.
    
  }
  \item{iteration}{
    The number of iterations in estimating parameters\cr
    Note: setting to a high number might cause a very long running time.
  }
  \item{auto.prune}{
   If set to TRUE, then function prune the orders to get the minimum orders. Recommended to set an initial order of 5 and 5 for AR and MA respectively.
  }
  \item{debug}{
    %%     ~~Describe \code{debug} here~~
      To see a detailed report of process set this parameter to TRUE.
  }
}
\details{ More details about the parameter estimation approach are given in Haselimashhadi and Vinciotti (2014).
          %%  ~~ If necessary, more details than the description above ~~
            %%   * Note that the number of replications in each group of data MUST be more that two if data are replicated.
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
\references{
  %% ~put references to the literature/web site here ~
    Haselimashhadi, H. and Vinciotti, V. (2014) Penalised estimation of autoregressive moving average models with time-dependent predictors (Working Paper). Brunel University London, UK
}
\author{
  %%  ~~who you are~~
    Hamed Haselimashhadi <hamed.haselimashhadi@brunel.ac.uk>
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
      \code{\link{sim.regARMA}}
  }
\examples{
  # --------------- FIRST EXAMPLE: regARMA without replications ----------------#  
  # True regression coefficients
  beta=c(1:6,rep(0,200))
  # Simulating data (100 time points) 
  simdata=sim.regARMA(n=100,beta=beta,phi=generateAR(2,min.distance = .3),theta=generateAR(1),
                      x.independent=TRUE,var.error=.1)
  #Fitting a REGARMA(2,1) model 
  regarma.fit=regARMA(data = simdata$rawdata,ar = 2,ma = 1,method = 'alasso', mselection = 'BIC')
  #Printing true and estimated coefficients
  cat('coefficients:',beta[1:6],',0,..., ',0)
  
  print(round(as.vector(regarma.fit$Betas),2))  
  print('--------------   OUTPUT structure --------------')
  str(regarma.fit)
  
  
  
  cat ("Press [enter] to see the second example");line <- readline()
  
  
  
  # ----------------------- SECOND EXAMPLE: regARMA with replications ------------------------#
  data(naics5809)
  R=10;dif1=2;dif2=0;ar=2;ma=1
  
  # Select the first 10 ids. 52 time points are available for each id
  # Note that naics5809 contains two response variables namely,
  # emp and pay. The analysis below considers pay for the response variable.
  # Data are sorted for each id (based on time). 
  
  data=as.matrix(naics5809)[1:(52*R),-(2:3)] 
  data[is.na(data)==TRUE]=0
  # ----------------------- Fitting model ------------------------# 
  regarma.fit=regARMA(data=data,ar=ar,ma=ma,method='alasso',Ndf=dif1,Sdf=dif2)
  # ----------------------- Printing outputs ------------------------# 
  print('--------------   OUTPUT structure --------------')
  str(regarma.fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
