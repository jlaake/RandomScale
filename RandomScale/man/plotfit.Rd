\name{plotfit}
\alias{plotfit}
\title{Random Scale Detection Function Plotting}
\usage{
  plotfit(x, w, par, nclass = NULL, weps = 5, dm = NULL,
    main = NULL, adjust = TRUE)
}
\arguments{
  \item{x}{vector of observed distances}

  \item{w}{half-width of strip}

  \item{nclass}{number of histogram bins}

  \item{par}{parameter vector (beta, beta_eps)}

  \item{weps}{range (-weps,weps) used for integration of
  std normal distribution}

  \item{dm}{the design matrix if one was used}

  \item{main}{title for plot}

  \item{adjust}{if TRUE, adjusts intercept}
}
\description{
  Plots histogram of data and fitted half-normal detection
  function with fixed or random scale for line transect
  sampling data. Also computes abundance if w!=Inf
}
\author{
  Jeff Laake
}

