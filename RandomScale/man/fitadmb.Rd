\name{fitadmb}
\alias{fitadmb}
\title{Random Scale Detection Function Fitting with admb}
\usage{
  fitadmb(x, w = Inf, formula = ~1, beta = NULL, sigma = 0,
    likelihood = "f2", extra.args = "-gh 10",
    verbose = TRUE, nsteps = 8, keep = FALSE,
    debug = FALSE)
}
\arguments{
  \item{x}{vector of distances or dataframe containing
  observed distances (distance) and other covariates}

  \item{w}{half-width of strip; if infinite w routine sets
  w to 2*max(x)}

  \item{formula}{formula for scale function}

  \item{beta}{starting values for beta}

  \item{sigma}{starting value for log sigma}

  \item{likelihood}{character string "g","f1","f2","fixed"}

  \item{extra.args}{for admb run}

  \item{verbose}{for compile and run}

  \item{nsteps}{adromb integration argument; default 8.}

  \item{keep}{if TRUE, uses existing tpl and exe}

  \item{debug}{if TRUE output parameter and -lnl values
  during iterations}
}
\description{
  Fits a half-normal detection function with random scale
  for line transect sampling data. Uses one of 3
  likelihoods. g,f1,f2,fixed
}
\details{
  Beta[1] and bounds for beta[1] are set based on the
  measured distance scale and they are adjusted internally
  by scaling distances such that max(x)=1. Likewise the
  likelihood that is shown when debug=T is for scaled
  distances but the reported value for the output is
  adjusted back to the original scale
}
\examples{
#fit simulated data with random scale
set.seed(123)
x=simdata(n=500,w=Inf,beta_eps=-.5)
par(mfrow=c(1,3))
glike=fitadmb(x,w=Inf,likelihood="g")
plotfit(x,w=Inf,c(glike$coeflist[[1]],glike$coeflist[[2]]),nclass=30,
		main="eq 4 likelihood")
f2like=fitadmb(x,w=Inf,likelihood="f2")
param=c(f2like$coeflist[[1]],f2like$coeflist[[2]])
plotfit(x,w=Inf,c(param[1]-exp(2*param[2]),param[2]),nclass=30,
		main="eq 7 likelihood")
f1like=fitadmb(x,w=Inf,likelihood="f1")
param=c(f1like$coeflist[[1]],f1like$coeflist[[2]])
plotfit(x,w=Inf,c(param[1]-exp(2*param[2]),param[2]),nclass=30,
		main="eq 8 likelihood")
dev.new()
par(mfrow=c(1,2))
#Mixed effect model
x1=simdata(n=2000,w=50,beta_eps=-.5,beta=2,fixed=FALSE,reject=TRUE)
x2=simdata(n=1000,w=50,beta_eps=-.5,beta=1,fixed=FALSE,reject=TRUE)
df=data.frame(covariate=c(rep(0,length(x1)),rep(1,length(x2))),
		distance=c(x1,x2))
fwlike=fitadmb(df,w=50,formula=~covariate,likelihood="f2")
param=c(fwlike$coeflist[[1]],fwlike$coeflist[[2]])
Nhatwcov=plotfit(df$distance,w=50,
		par=c(param[1]-exp(2*param[3]),param[2],param[3]),
		nclass=30,dm=model.matrix(~covariate,df),main="With covariate")
flike=fitadmb(df,w=50,formula=~1,likelihood="f2")
param=c(flike$coeflist[[1]],flike$coeflist[[2]])
Nhatwocov=plotfit(df$distance,w=50,
		par=c(param[1]-exp(2*param[2]),param[2]),nclass=30,
		main="Without covariate")
dev.new()
par(mfrow=c(2,2))
param=c(fwlike$coeflist[[1]],fwlike$coeflist[[2]])
Nhatwcov0=plotfit(df$distance[df$covariate==0],w=50,
		par=c(param[1]-exp(2*param[3]),param[2],param[3]),
		nclass=30,dm=model.matrix(~covariate,df[df$covariate==0,]),
		main="With covariate value=0")
Nhatwcov1=plotfit(df$distance[df$covariate==1],w=50,
		par=c(param[1]-exp(2*param[2]),param[2],param[3]),nclass=30,
		dm=model.matrix(~covariate,df[df$covariate==1,]),
		main="With covariate value=1")
param=c(flike$coeflist[[1]],flike$coeflist[[2]])
Nhatwocov0=plotfit(df$distance[df$covariate==0],w=50,
		par=c(param[1]-exp(2*param[2]),param[2]),nclass=30,
		main="Without covariate value=0")
Nhatwocov1=plotfit(df$distance[df$covariate==1],w=50,
		par=c(param[1]-exp(2*param[2]),param[2]),nclass=30,
		main="Without covariate value=1")
}
\author{
  Jeff Laake
}

