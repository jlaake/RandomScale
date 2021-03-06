% Generated by roxygen2 (4.0.1): do not edit by hand
\name{compute_Nhat}
\alias{compute_Nhat}
\title{Random Scale Detection Function Abundance}
\usage{
compute_Nhat(par, x, w, weps = 5, dm = NULL, both = FALSE,
  adjust = TRUE)
}
\arguments{
\item{par}{parameter vector (beta, beta_eps)}

\item{x}{vector of observed distances}

\item{w}{half-width of strip}

\item{weps}{range (-weps,weps) used for integration of std normal distribution}

\item{dm}{the design matrix if one was used}

\item{both}{if TRUE returns both avg_mu_est and Nhat in a list; otherwise just Nhat}

\item{adjust}{if TRUE, adjusts intercept}
}
\description{
Compute estimate of N in the covered area.
}
\author{
Jeff Laake
}

