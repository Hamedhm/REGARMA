\name{sim.regARMA}
\alias{sim.regARMA}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  %%  ~~function to do ... ~~
   A function to simulate data from a REGARMA model
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    This function simulates a mean zero Gaussian lagged regression variable in the presence
    of autocorrelationS among residuals (REGARMA)
}
\usage{
  sim.regARMA(n = 5, beta = c(0.62), x.independent = TRUE, 
    phi = c(0.3), theta = c(0.5), var.error = 1, n.z.coeffs=0,
    shuffle=FALSE, draw.plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    %%     ~~Describe \code{n} here~~
      The number of datapoints to be simulated
  }
  \item{beta}{
    %%     ~~Describe \code{beta} here~~
      Regression coefficients\cr
      * Note that coefficients must be in [-1,1] for time-series data (x.independent=FALSE). 
  }
  \item{x.independent}{
    %%     ~~Describe \code{x.independent} here~~
      Set this option to true if regression variables are assumed to be independent
  }
  \item{phi}{
    %%     ~~Describe \code{phi} here~~
      Autoregressive coefficients\cr
    * Note that, for stationarity, all roots of the polynomial must site outside of a unit circle
  }
  \item{theta}{
    %%     ~~Describe \code{theta} here~~
      Residuals coefficients\cr
    * Note that, for stationarity, all roots of the polynomial include coefficients must site outside of a unit circle
  }
  \item{var.error}{
    %%     ~~Describe \code{var.error} here~~
      Variance of the process
  }
  \item{n.z.coeffs}{
    The number of zero coefficients if needed
  }
  \item{shuffle}{
    Shuffle coefficients if set to TRUE. Otherwise the first coefficients are non-zero and the rest are zero.
  }
  \item{draw.plot}{
    %%     ~~Describe \code{draw.plot} here~~
      If TRUE, a time series plot is drawn
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
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
    Hamed Haselimashhadi <hamed.haselimashhadi@brunel.ac.u>
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
    %%\code{\link{regarma}}, \code{\link{iregarma}}
    \code{\link{regARMA}}
  }
\examples{
  simdata=sim.regARMA(n = 100,beta = .1,x.independent = TRUE,phi = .4,
    theta = -.4,var.error = 1,draw.plot = TRUE)
  str(simdata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
