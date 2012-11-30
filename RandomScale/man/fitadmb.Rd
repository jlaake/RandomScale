\name{fitadmb}
\alias{fitadmb}
\title{Random Scale Detection Function Fitting with admb}
\usage{
  fitadmb(x, w, formula = ~1, likelihood = "f2",
    extra.args = "-est -gh 10", verbose = FALSE)
}
\arguments{
  \item{x}{vector of distances or dataframe containing
  observed distances (distance) and other covariates}

  \item{w}{half-width of strip}

  \item{formula}{formula for scale function}

  \item{likelihood}{character string "g","f1","f2"}

  \item{extra.args}{for admb run}

  \item{verbose}{for compile and run}
}
\description{
  Fits a half-normal detection function with random scale
  for line transect sampling data. Uses one of 3
  likelihoods. g,f1,f2
}
\author{
  Jeff Laake
}

