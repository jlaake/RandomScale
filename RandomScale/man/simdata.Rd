\name{simdata}
\alias{simdata}
\title{Simulates data from random scale half-normal detection function}
\usage{
  simdata(n = 500, beta = 2, beta_eps = -1, w = Inf,
    reject = TRUE, b = 10000, fixed = TRUE)
}
\arguments{
  \item{n}{sample size}

  \item{beta}{scale intercept of log(sigma)}

  \item{beta_eps}{log(sigma_epsilon)}

  \item{w}{half-width of line transect}

  \item{reject}{if TRUE uses brute force rejection sampling
  to generate distances;if w=Inf, it uses a large w to
  avoid any appreciable truncation}

  \item{b}{number of deviates generated in a batch for
  rejection sampling}

  \item{fixed}{if TRUE and reject==TRUE, returns n
  observations; otherwise if reject=TRUE uses n as sample
  of bernoulli trials}
}
\description{
  If reject==FALSE, it generates distances directly with
  rnorm but if this approach is used, the estimated beta
  with the correct likelihood will not match the generating
  beta but this is not bias in the usual sense. It will
  match the beta from the incorrect likelihood.
}
\author{
  Jeff Laake
}

