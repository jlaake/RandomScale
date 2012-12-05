#' Half-normal detection and related functions
#' 
#' @usage gx(x,sigma)
#' 
#'        gx.eps(x,eps,par,dm=NULL)
#' 
#'        avg_gx(x,par,weps=5,dm=NULL)
#' @export gx gx.eps avg_gx
#' @aliases gx gx.eps avg_gx       
#' @param x distance
#' @param sigma half-normal scale (std deviation)
#' @param eps std normal deviate
#' @param par parameter vector (beta, beta_eps)
#' @param dm design matrix for fixed effect parameters if one is used
#' @param weps range (-weps,weps) for std normal integration 
#' @author Jeff Laake
# Half-normal detection function
gx=function(x,sigma)
{
	exp(-.5*(x/sigma)^2)
}
# Half-normal * density(eps)
gx.eps=function(x,eps,par,dm=NULL)
{
	if(is.null(dm))
	   sigma=exp(par[1]+eps*exp(par[2]))
   else
	   sigma=exp(dm%*%par[1:ncol(dm)]+eps*exp(par[ncol(dm)+1]))
	gx(x,sigma)*dnorm(eps,0,1)
}
# Average half-normal integrated over distribution of eps
avg_gx=function(x,par,weps=5,dm=NULL)
{
	gx=vector("numeric",length=length(x))
	geps.x=function(eps,x,par,dm)gx.eps(x,eps,par,dm)
	for(i in 1:length(x))
		gx[i]=integrate(geps.x,-weps,weps,x=x[i],par=par,dm=dm)$value
	gx
}
