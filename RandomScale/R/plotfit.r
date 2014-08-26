#' Random Scale Detection Function Plotting
#' 
#' Plots histogram of data and fitted half-normal detection function with fixed or random scale for line
#' transect sampling data. Also computes abundance if w!=Inf
#' 
#' @param x vector of observed distances
#' @param w half-width of strip 
#' @param nclass number of histogram bins
#' @param par parameter vector (beta, beta_eps)
#' @param weps range (-weps,weps) used for integration of std normal distribution
#' @param dm the design matrix if one was used
#' @param main title for plot
#' @param adjust if TRUE, adjusts intercept
#' @export plotfit
#' @author Jeff Laake
plotfit=function(x,w,par,nclass=NULL,weps=5,dm=NULL,main=NULL,adjust=TRUE)
{
#  Compute Nhat for plotting
   Nlist=compute_Nhat(par,x,w,weps=5,dm=dm,both=TRUE,adjust=adjust)
   Nhat=Nlist$Nhat
   avg_mu_est=Nlist$avg_mu_est
# Create plot of fit; uses histline from mrds
	if(adjust)
		par[1]=par[1]-exp(2*par[length(par)])
	if(w==Inf)w=max(x)
	max_x=w
	if(is.null(nclass))
		nints=ceiling(sqrt(length(x)))
	else
		nints=nclass
	int_width=max_x/nints
	breaks=int_width*(0:nints)
	ints=(0:100)*max_x/100
	hh=hist(x,plot=F,breaks=breaks)
	expected.counts=int_width*Nhat/max_x
	bars=hh$counts/expected.counts
	if(length(par)<2)
	{
		gx=gx(ints,sigma=exp(par))
	}else
	{
		if(is.null(dm))
		{
			gx=avg_gx(ints,par=par)
		}else
		{
			gx=vector("numeric",length=length(ints))
			gx=0
			for(i in 1:nrow(dm))
				gx=gx+avg_gx(ints,par=par,dm=dm[i,,drop=FALSE])*((w/avg_mu_est[i])/Nhat)
		}
    }  
	histline(bars,breaks,ylim=c(0,max(c(1,bars))),
			xlab="Distance",ylab="Detection probability",main=main)
	lines(ints,gx)
	invisible(Nhat)
}
#' Random Scale Detection Function Abundance
#' 
#' Compute estimate of N in the covered area.
#' 
#' @param par parameter vector (beta, beta_eps)
#' @param x vector of observed distances
#' @param w half-width of strip 
#' @param weps range (-weps,weps) used for integration of std normal distribution
#' @param dm the design matrix if one was used
#' @param both if TRUE returns both avg_mu_est and Nhat in a list; otherwise just Nhat
#' @param adjust if TRUE, adjusts intercept
#' @export compute_Nhat
#' @author Jeff Laake
compute_Nhat=function(par,x,w,weps=5,dm=NULL,both=FALSE,adjust=TRUE)
{
if(length(par)<2)
{
	# compute fixed mu
	avg_mu_est=exp(par)*sqrt(2*pi)*
			(pnorm(w,0,exp(par))-.5)
}else
{
	if(adjust)
	   par[1]=par[1]-exp(2*par[length(par)])
	# Compute average_mu
	if(is.null(dm))
		avg_mu_est=avg_mu(par=par,w=w,weps=weps,dm=NULL)
	else
	{
		avg_mu_est=vector("numeric",length=nrow(dm))
		for(i in 1:nrow(dm))
			avg_mu_est[i]=avg_mu(par=par,w=w,weps=weps,dm=dm[i,,drop=FALSE])
	}
}
# Compute Nhat
if(w==Inf)w=max(x)
if(w!=Inf)
{
	if(is.null(dm))
	{
		Nhat=length(x)/(avg_mu_est/w)
	}else
	{
		Nhat=sum(w/avg_mu_est)
	}
}else
	Nhat=NULL
if(both)
   return(list(Nhat=Nhat,avg_mu_est=avg_mu_est))
else
	return(Nhat)
}

#' Random Scale Detection Function Abundance Standard Error
#' 
#' Compute estimate of se(N) in the covered area; does not account for variation in n
#' 
#' @param par parameter vector (beta, beta_eps)
#' @param vcov variance-covariance matrix of parameters (beta, beta_eps)
#' @param x vector of observed distances
#' @param w half-width of strip 
#' @param weps range (-weps,weps) used for integration of std normal distribution
#' @param dm the design matrix if one was used
#' @param adjust if TRUE, intercept is adjusted
#' @export compute_Nhat.se
#' @author Jeff Laake
compute_Nhat.se=function(par,vcov,x,w,weps=5,dm=NULL,adjust=TRUE)
{
   Nlist=compute_Nhat(par,x,w,weps=5,dm=dm,both=TRUE,adjust=adjust)
   p=Nlist$avg_mu_est/w
   if(length(p)==1)
   return(sqrt(Nlist$Nhat*(1-p)/p+
	     DeltaMethod(par,compute_Nhat,vcov,delta=0.001,x=x,w=w,weps=weps,dm=dm,both=FALSE,adjust=adjust)$variance))
   else
	   return(sqrt(sum((1-p)/p^2)+
	      DeltaMethod(par,compute_Nhat,vcov,delta=0.001,x=x,w=w,weps=weps,dm=dm,both=FALSE,adjust=adjust)$variance))
   
}

#' Numeric Delta Method approximation for the variance-covariance matrix (from mrds package)
#'
#' Computes delta method variance-covariance matrix of results of any generic
#' function \code{fct} that computes a vector of estimates as a function of a
#' set of estimated parameters \code{par}.
#'
#' The delta method (aka propogation of
#' errors is based on Taylor series approximation - see Seber's book on
#' Estimation of Animal Abundance).  It uses the first derivative of \code{fct}
#' with respect to \code{par} which is computed in this function numerically
#' usign the central-difference formula.  It also uses the variance-covariance
#' matrix of the estimated parameters which is derived in estimating the
#' parameters and is an input argument.
#'
#' The first argument of \code{fct} should be \code{par} which is a vector of
#' parameter estimates. It should return a single value (or vector) of
#' estimate(s).  The remaining arguments of \code{fct} if any can be passed to
#' \code{fct} by including them at the end of the call to \code{DeltaMethod} as
#' name=value pairs separated by commas as with all arguments.
#'
#' @param par vector of parameter values at which estimates should be
#'   constructed
#' @param fct function that constructs estimates from parameters \code{par}
#' @param vcov variance-covariance matrix of the parameters
#' @param delta proportional change in parameters used to numerically estimate
#'   first derivative with central-difference formula
#' @param \dots any additional arguments needed by \code{fct}
#' @export
#' @return a list with values \item{variance}{estimated variance-covariance
#'   matrix of estimates derived by \code{fct}} \item{partial}{ matrix (or
#'   vector) of partial derivatives of \code{fct} with respect to the
#'   parameters \code{par}}
#' @note This is a generic function that can be used in any setting beyond the
#'   mrds package.  However this is an internal function for mrds and the user
#'   does not need to call it explicitly.
#' @author Jeff Laake
#' @keywords utility
DeltaMethod <- function(par, fct, vcov, delta, ...){
	
	# Construct theta call to fct
	theta <- function(par) fct(par,...)
	
	# Numerically compute the first derivative of the estimator with
	# respect to each parameter using a central difference formula
	savepar <- par
	value1 <- theta(par)
	partial <- matrix(0,nrow=length(par),ncol=length(value1))
	for(i in 1:length(par)){
		# Store the original parameters into par and then adjust the
		# ith parameter by adding the proportion based on delta
		if(savepar[i]!=0){
			deltap <- delta*savepar[i]
		}else{
			deltap <- delta
		}
		
		par <- savepar
		par[i] <- savepar[i]+deltap
		
		# With this new value call the function to compute the estimate
		value1 <- theta(par)
		
		# Next do the same thing as above except substract off the proportion
		# delta from the original parameter.
		
		par <- savepar
		par[i] <- savepar[i]-deltap
		value2 <- theta(par)
		
		# Compute the central difference formula for the first derivative
		# with respect to the ith parameter
		
		partial[i,] <- (value1-value2)/(2*deltap)
	}
	variance <- t(partial)%*%vcov%*%partial
	
	# return the v-c matrix and the first partial vector(matrix)
	return(list(variance = variance,
					partial  = partial))
}


