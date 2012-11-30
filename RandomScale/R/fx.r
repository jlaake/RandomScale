#' Half-normal pdf and related functions
#' 
#' @usage fx(x,sigma,w)
#' 
#'        fx.eps(eps, x, par, w)
#' 
#' @aliases fx fx.eps       
#' @export fx fx.eps
#' @param x distance
#' @param sigma half-normal scale (std deviation)
#' @param eps std normal deviate
#' @param par parameter vector (beta, beta_eps)
#' @param w half-width of transect
#' @param weps range (-weps,weps) for std normal integration 
#' @author Jeff Laake# pdf - for incorrect likelihood 
fx=function(x,sigma,w)
{
	mu=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)
	exp(-.5*(x/sigma)^2)/mu
}
# pdf * density(eps) - for incorrect likelihood 
fx.eps=function(eps,x,par,w)
{
	sigma=exp(par[1]+eps*exp(par[2]))
	fx(x,sigma,w)*dnorm(eps,0,1)
}
