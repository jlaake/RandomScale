\name{fitdata}
\alias{fitdata}
\alias{flnl}
\title{Random Scale Detection Function Fitting}
\usage{
  fitdata(x,beta=2,beta_eps=-1,w=Inf,weps=5,
    lower=beta/2,upper=2*beta,wrong=FALSE)
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
\author{
  Jeff Laake
}

