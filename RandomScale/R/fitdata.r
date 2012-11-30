#' Random Scale Detection Function Fitting
#' 
#' Fits a half-normal detection function with fixed or random scale for line
#' transect sampling data. flnl computes negative log-likelihoo used from fitdata
#' through optim or optimize for fixed effect only.
#' 
#' @usage fitdata(x,beta=2,beta_eps=-1,w=Inf,weps=5,
#'             lower=beta/2,upper=2*beta,wrong=FALSE)
#' @usage flnl(par,x,w,weps=5,wrong=FALSE)
#' 
#' @aliases fitdata flnl
#' @export fitdata flnl
#' @param x vector of observed distances
#' @param beta mean value for log(sigma)
#' @param beta_eps log(sigma_epsilon); if set to NULL uses fixed effect
#' @param par parameter vector for flnl
#' @param w half-width of strip 
#' @param weps range (-weps,weps) used for integration of std normal distribution
#' @param lower lower bound for beta in fixed effect model for optimize
#' @param upper upper bound for beta in fixed effect model for optimize
#' @param wrong if TRUE uses incorrect likelihood
#' @author Jeff Laake
fitdata=function(x,beta=2,beta_eps=-1,w=Inf,
		weps=5,lower=beta/2,upper=2*beta,
		wrong=FALSE)
{
	n=length(x)
	if(is.null(beta_eps))
	{
		# fit fixed scale model
		model=optimize(flnl,lower=lower,upper=upper,
				x=x,w=w,weps=weps)
		model$par=model$minimum
	}else
	{
#   fit random scale model 
		model=optim(par=c(beta,beta_eps),
				flnl,x=x,w=w,weps=weps,wrong=wrong)
	}
	model$x=x
	model$w=w
	return(list(model=model))
}
# Negative log-likelihood function 
flnl=function(par,x,w,weps=5,wrong=FALSE)
{
	lnl=0
	for(i in 1:length(x))
	{
		if(length(par)>1)
			if(wrong)
				lnl=lnl-log(integrate(fx.eps,-weps,weps,
								x=x[i],w=w,par=par)$value)
			else
				lnl=lnl-log(avg_gx(x[i],par)/
								avg_mu(w,par,weps))
		else
		{
			sigma=exp(par[1])
			mu=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)
			lnl=lnl-log(gx(x[i],sigma=sigma))+log(mu)
		}
	}
	lnl
}
