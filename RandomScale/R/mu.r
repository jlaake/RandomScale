#' Integration of Half-normal detection 
#' 
#' @usage mu(eps,w,par,dm=NULL)
#' 
#'      avg_mu(w,par,weps=5,dm=NULL)
#' @export mu avg_mu
#' @aliases mu avg_mu       
#' @param eps std normal deviate
#' @param par parameter vector (beta, beta_eps)
#' @param dm design matrix for fixed effect if one is used
#' @param w transect half-width
#' @param weps range (-weps,weps) for std normal integration 
#' @author Jeff Laake
mu=function(eps,w,par,dm=NULL)
{
# Average detection probability * density eps
	mu=vector("numeric",length=length(eps))
	for(i in 1:length(eps))
	{
		if(is.null(dm))
		   sigma=exp(par[1]+eps[i]*exp(par[2]))
	    else
			sigma=exp(dm%*%par[1:ncol(dm)]+eps[i]*exp(par[ncol(dm)+1]))
		mu[i]=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)*dnorm(eps[i],0,1)
	}
	mu
}
# Average detection probability - integrated over eps
avg_mu=function(w,par,weps=5,dm=NULL)
	integrate(mu,-weps,weps,w=w,par=par,dm=dm)$value
