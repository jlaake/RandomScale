\name{mu}
\alias{avg_mu}
\alias{mu}
\title{Integration of Half-normal detection}
\usage{
  mu(eps,w,par,dm=NULL)

  avg_mu(w,par,weps=5,dm=NULL)
}
\arguments{
  \item{eps}{std normal deviate}

  \item{par}{parameter vector (beta, beta_eps)}

  \item{dm}{design matrix for fixed effect if one is used}

  \item{w}{transect half-width}

  \item{weps}{range (-weps,weps) for std normal
  integration}
}
\description{
  Integration of Half-normal detection
}
\author{
  Jeff Laake
}

