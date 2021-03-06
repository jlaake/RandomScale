% Generated by roxygen2 (4.0.1): do not edit by hand
\name{compute_Nhat.se}
\alias{compute_Nhat.se}
\title{Random Scale Detection Function Abundance Standard Error}
\usage{
compute_Nhat.se(par, vcov, x, w, weps = 5, dm = NULL, adjust = TRUE)
}
\arguments{
\item{par}{parameter vector (beta, beta_eps)}

\item{vcov}{variance-covariance matrix of parameters (beta, beta_eps)}

\item{x}{vector of observed distances}

\item{w}{half-width of strip}

\item{weps}{range (-weps,weps) used for integration of std normal distribution}

\item{dm}{the design matrix if one was used}

\item{adjust}{if TRUE, intercept is adjusted}
}
\description{
Compute estimate of se(N) in the covered area; does not account for variation in n
}
\author{
Jeff Laake
}

