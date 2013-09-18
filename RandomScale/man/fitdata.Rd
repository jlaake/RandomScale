\name{fitdata}
\alias{fitdata}
\alias{flnl}
\title{Random Scale Detection Function Fitting}
\usage{
  fitdata(x,beta=NULL,beta_eps=-3,bounds=NULL,w=Inf,weps=5,
    wrong=FALSE,debug=FALSE,method="L-BFGS-B",hessian=FALSE)

  flnl(par,x,w,weps=5,wrong=FALSE,debug)
}
\arguments{
  \item{x}{vector of observed distances}

  \item{beta}{mean value for log(sigma)}

  \item{beta_eps}{log(sigma_epsilon); if set to NULL uses
  fixed effect}

  \item{bounds}{2 x 2 matrix; rows: beta,beta_eps; columns:
  lower,upper bounds; only 1x2 if beta_eps=NULL}

  \item{par}{parameter vector for flnl}

  \item{w}{half-width of strip}

  \item{weps}{range (-weps,weps) used for integration of
  std normal distribution}

  \item{wrong}{if TRUE uses incorrect likelihood}

  \item{debug}{if TRUE output parameter and -lnl values
  during iterations}

  \item{method}{optim method; if bounds specified uses
  L-BFGS-B regardless}

  \item{hessian}{if TRUE, returns hessian for v-c matrix}
}
\description{
  Fits a half-normal detection function with fixed or
  random scale for line transect sampling data. flnl
  computes negative log-likelihood used from fitdata
  through optim or optimize for fixed effect only.
}
\details{
  Beta and bounds for beta are set based on the measured
  distance scale and they are adjusted internally by
  scaling distances such that max(x)=1 Likewise the
  likelihood that is shown when debug=T is for scaled
  distances but the reported value for the output is
  adjusted back to the original scale If method="L-BFGS-B"
  and no bounds are set, default scaled values are used:
  bounds=matrix(c(-3,2,-10,2),nrow=2,byrow=TRUE)
}
\examples{
#simulate some data using rejection sampling and n=500
set.seed(123)
x=simdata(n=500,w=Inf,beta_eps=-.5)
#fit and plot model with g and f1 likelihood
par(mfrow=c(1,3))
results_random=fitdata(x,w=Inf,beta_eps=-.5)
plotfit(x,w=max(x),results_random$model$par,nclass=30,
                main="eq 4 likelihood",adjust=FALSE)
#Because data are untruncated the estimates from the f1 likelihood will match
# once the intercept is adjusted
results_random_wrong=fitdata(x,w=Inf,beta_eps=-.5,wrong=TRUE)
plotfit(x,w=max(x),results_random_wrong$model$par,nclass=30,
      main="eq 8 likelihood",adjust=FALSE)
plotfit(x,w=max(x),results_random_wrong$model$par,nclass=30,
      main="eq 8 likelihood\\nadjusted beta")
}
\author{
  Jeff Laake
}

