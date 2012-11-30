\name{fitadmb}
\alias{fitadmb}
\title{Random Scale Detection Function Fitting with admb}
\usage{
  fitadmb(x, w, likelihood = "f2",
    extra.args = "-est -gh 10", verbose = FALSE)
}
\arguments{
  \item{x}{vector of observed distances}

  \item{w}{half-width of strip}

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

