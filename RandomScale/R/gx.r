#' Half-normal detection and related functions
#' 
#' @usage gx(x,sigma)
#' 
#'        gx.eps(x,eps,par)
#' 
#'        avg_gx(x,par,weps=5)
#' @export gx gx.eps avg_gx
#' @aliases gx gx.eps avg_gx       
#' @param x distance
#' @param sigma half-normal scale (std deviation)
#' @param eps std normal deviate
#' @param par parameter vector (beta, beta_eps)
#' @param weps range (-weps,weps) for std normal integration 
#' @author Jeff Laake
# Half-normal detection function
gx=function(x,sigma)
{
	exp(-.5*(x/sigma)^2)
}
# Half-normal * density(eps)
gx.eps=function(x,eps,par)
{
	sigma=exp(par[1]+eps*exp(par[2]))
	gx(x,sigma)*dnorm(eps,0,1)
}
# Average half-normal integrated over distribution of eps
avg_gx=function(x,par,weps=5)
{
	gx=vector("numeric",length=length(x))
	geps.x=function(eps,x,par)gx.eps(x,eps,par)
	for(i in 1:length(x))
		gx[i]=integrate(geps.x,-weps,weps,x=x[i],par=par)$value
	gx
}
