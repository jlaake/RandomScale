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
#' @param main title for plot
#' @export plotfit
#' @author Jeff Laake
plotfit=function(x,w,par,nclass=NULL,weps=5,main=NULL)
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
		avg_mu_est=avg_mu(w=w,par=par,weps=weps)
	}
# Compute Nhat if W<Inf
	if(w!=Inf)
	{
		Nhat=length(x)/(avg_mu_est/w)
	}else
		Nhat=NULL
	max_x=max(x)
	if(is.null(nclass))
		nints=ceiling(sqrt(length(x)))
	else
		nints=nclass
	int_width=max_x/nints
	breaks=int_width*(0:nints)
	ints=(0:100)*max_x/100
	hh=hist(x,plot=F,breaks=breaks)
	bars=(hh$counts/(sum(hh$counts)*int_width))*avg_mu_est
	if(length(par)<2)
	{
		gx=gx(ints,sigma=exp(par))
	}else
		gx=avg_gx(ints,par=par)
	histline(bars,breaks,ylim=c(0,max(c(1,bars))),
			xlab="Distance",ylab="Detection probability",main=main)
	lines(ints,gx)
	invisible(Nhat)
}
