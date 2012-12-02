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
#' @export plotfit
#' @author Jeff Laake
plotfit=function(x,w,par,nclass=NULL,weps=5,dm=NULL,main=NULL)
{
# Create plot of fit; uses histline from mrds
	if(length(par)<2)
	{
		# compute fixed mu
		avg_mu_est=exp(par)*sqrt(2*pi)*
				(pnorm(w,0,exp(par))-.5)
	}else
	{
		# Compute average_mu
		if(is.null(dm))
		  avg_mu_est=avg_mu(w=w,par=par,weps=weps,dm=NULL)
	    else
		{
			avg_mu_est=vector("numeric",length=nrow(dm))
			for(i in 1:nrow(dm))
				avg_mu_est[i]=avg_mu(w=w,par=par,weps=weps,dm=dm[i,,drop=FALSE])
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
