\name{gx}
\alias{avg_gx}
\alias{gx}
\alias{gx.eps}
\title{Half-normal detection and related functions}
\usage{
  gx(x,sigma)

  gx.eps(x,eps,par)

  avg_gx(x,par,weps=5)
}
\arguments{
  \item{x}{distance}

  \item{sigma}{half-normal scale (std deviation)}

  \item{eps}{std normal deviate}

  \item{par}{parameter vector (beta, beta_eps)}

  \item{weps}{range (-weps,weps) for std normal
  integration}
}
\description{
  Half-normal detection and related functions
}
\author{
  Jeff Laake
}

