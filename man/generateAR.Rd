\name{generateAR}
\alias{generateAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
    A function to generate autoregressive coefficients
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    A function to generate autoregressive coefficients in a certain interval and with a mimum distance amongst any two coefficients.
}
\usage{
  generateAR(n = 1, l = -1, u = 1, min.distance = .Machine$double.eps, sort.coeff=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    %%     ~~Describe \code{n} here~~
      The number of coefficients
  }
  \item{l}{
    %%     ~~Describe \code{l} here~~
      Lower bound of coefficients
  }
  \item{u}{
    %%     ~~Describe \code{u} here~~
      Upper bound of coefficients
  }
  \item{min.distance}{
    %%     ~~Describe \code{min.distance} here~~
      Minimum distance among each two coefficients
  }
  \item{sort.coeff}{
      If TRUE then function generates ordered coefficients
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
    Hamed Haselimashhadi <hamed.haselimashhadi@brunel.ac.uk>
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
    \code{\link{regARMA}}
  }
\examples{
  generateAR(5,-1,1,.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
