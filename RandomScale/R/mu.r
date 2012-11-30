#' Integration of Half-normal detection 
#' 
#' @usage mu(eps,w,par)
#' 
#'      avg_mu(w,par,weps=5)
#' @export mu avg_mu
#' @aliases mu avg_mu       
#' @param x distance
#' @param sigma half-normal scale (std deviation)
#' @param eps std normal deviate
#' @param par parameter vector (beta, beta_eps)
#' @param weps range (-weps,weps) for std normal integration 
#' @author Jeff Laake
mu=function(eps,w,par)
{
# Average detection probability * density eps
	mu=vector("numeric",length=length(eps))
	for(i in 1:length(eps))
	{
		sigma=exp(par[1]+eps[i]*exp(par[2]))
		mu[i]=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)*dnorm(eps[i],0,1)
	}
	mu
}
# Average detection probability - integrated over eps
avg_mu=function(w,par,weps=5)
	integrate(mu,-weps,weps,w=w,par=par)$value
