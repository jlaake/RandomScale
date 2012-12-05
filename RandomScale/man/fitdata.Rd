\name{fitdata}
\alias{fitdata}
\alias{flnl}
\title{Random Scale Detection Function Fitting}
\usage{
  fitdata(x,beta=2,beta_eps=-1,w=Inf,weps=5,
    lower=beta/2,upper=2*beta,wrong=FALSE)

  flnl(par,x,w,weps=5,wrong=FALSE)
}
\arguments{
  \item{x}{vector of observed distances}

  \item{beta}{mean value for log(sigma)}

  \item{beta_eps}{log(sigma_epsilon); if set to NULL uses
  fixed effect}

  \item{par}{parameter vector for flnl}

  \item{w}{half-width of strip}

  \item{weps}{range (-weps,weps) used for integration of
  std normal distribution}

  \item{lower}{lower bound for beta in fixed effect model
  for optimize}

  \item{upper}{upper bound for beta in fixed effect model
  for optimize}

  \item{wrong}{if TRUE uses incorrect likelihood}
}
\description{
  Fits a half-normal detection function with fixed or
  random scale for line transect sampling data. flnl
  computes negative log-likelihoo used from fitdata through
  optim or optimize for fixed effect only.
}
\examples{
#simulate some data using rejection sampling and n=500
set.seed(123)
x=simdata(n=500,w=Inf,beta_eps=-.5)
#fit and plot model with g and f1 likelihood
par(mfrow=c(1,3))
results_random=fitdata(x,w=Inf,beta_eps=-.5)
plotfit(x,w=max(x),results_random$model$par,nclass=30,
                main="eq 4 likelihood")
#Because data are untruncated the estimates from the f1 likelihood
results_random_wrong=fitdata(x,w=Inf,beta_eps=-.5,wrong=TRUE)
param=results_random_wrong$model$par
plotfit(x,w=max(x),param,nclass=30,main="eq 8 likelihood")
plotfit(x,w=max(x),c(param[1]-exp(2*param[2]),param[2]),nclass=30,
      main="eq 8 likelihood\\nadjusted beta")
}
\author{
  Jeff Laake
}

