#' Simulates data from random scale half-normal detection function
#' 
#' If reject==FALSE, it generates distances directly with 
#' rnorm but if this approach is used, the estimated beta with 
#' the correct likelihood will not match the generating beta but 
#' this is not bias in the usual sense. It will match the beta
#' from the incorrect likelihood.
#' @param n sample size
#' @param beta scale intercept of log(sigma)
#' @param beta_eps log(sigma_epsilon)
#' @param w half-width of line transect
#' @param reject if TRUE uses brute force rejection sampling to 
#'          generate distances;if w=Inf, it uses a large w to avoid any appreciable truncation
#' @param b number of deviates generated in a batch for rejection sampling
#' @param fixed if TRUE and reject==TRUE, returns n observations; otherwise if reject=TRUE uses n as sample of bernoulli trials
#' @export         
#' @author Jeff Laake
simdata=function(n=500,beta=2,beta_eps=-1,w=Inf,
		reject=TRUE,b=10000,fixed=TRUE)
{
	if(reject)
	{
		if(w==Inf)
			w=3*exp(beta+3*exp(beta_eps))
		if(fixed)
		{
			x=NULL
			while(length(x)<n)
			{
				u=runif(b,0,w)
				sigma=exp(beta+rnorm(b,0,1)*exp(beta_eps))
				seen=gx(u,sigma)>runif(b,0,1)
				x=c(x,u[seen])
			}
			return(x[1:n])
		} else
		{
			u=runif(n,0,w)
			sigma=exp(beta+rnorm(n,0,1)*exp(beta_eps))
			seen=gx(u,sigma)>runif(n,0,1)
			return(u[seen])
		}
	}else
	{
		x=NULL
		while(length(x)<n)
		{
			sigma=exp(beta+rnorm(b,0,1)*exp(beta_eps))
			z=abs(rnorm(n,0,sigma))
			x=c(x,z[z<w])
		}
		return(x[1:n])
	}
}
