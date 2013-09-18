#' Random Scale Detection Function Fitting
#' 
#' Fits a half-normal detection function with fixed or random scale for line
#' transect sampling data. flnl computes negative log-likelihood used from fitdata
#' through optim or optimize for fixed effect only.
#' 
#' Beta and bounds for beta are set based on the measured distance scale and
#' they are adjusted internally by scaling distances such that max(x)=1
#' Likewise the likelihood that is shown when debug=T is for scaled distances
#' but the reported value for the output is adjusted back to the original scale
#' If method="L-BFGS-B" and no bounds are set, default scaled values are
#' used: bounds=matrix(c(-3,2,-10,2),nrow=2,byrow=TRUE)
#' 
#' @usage fitdata(x,beta=NULL,beta_eps=-3,bounds=NULL,w=Inf,weps=5,
#'             wrong=FALSE,debug=FALSE,method="L-BFGS-B",hessian=FALSE)
#' 
#'        flnl(par,x,w,weps=5,wrong=FALSE,debug)
#' 
#' @aliases fitdata flnl
#' @export fitdata flnl
#' @param x vector of observed distances
#' @param beta mean value for log(sigma)
#' @param beta_eps log(sigma_epsilon); if set to NULL uses fixed effect
#' @param bounds 2 x 2 matrix; rows: beta,beta_eps; columns: lower,upper bounds; only 1x2 if beta_eps=NULL
#' @param par parameter vector for flnl
#' @param w half-width of strip 
#' @param weps range (-weps,weps) used for integration of std normal distribution
#' @param wrong if TRUE uses incorrect likelihood
#' @param debug if TRUE output parameter and -lnl values during iterations
#' @param method optim method; if bounds specified uses L-BFGS-B regardless
#' @param hessian if TRUE, returns hessian for v-c matrix
#' @author Jeff Laake
#' @examples
#' #simulate some data using rejection sampling and n=500
#' set.seed(123) 
#' x=simdata(n=500,w=Inf,beta_eps=-.5)
#' #fit and plot model with g and f1 likelihood
#' par(mfrow=c(1,3)) 
#' results_random=fitdata(x,w=Inf,beta_eps=-.5)
#' plotfit(x,w=max(x),results_random$model$par,nclass=30,
#'                 main="eq 4 likelihood",adjust=FALSE)
#' #Because data are untruncated the estimates from the f1 likelihood will match
#' # once the intercept is adjusted
#' results_random_wrong=fitdata(x,w=Inf,beta_eps=-.5,wrong=TRUE)
#' plotfit(x,w=max(x),results_random_wrong$model$par,nclass=30,
#'       main="eq 8 likelihood",adjust=FALSE)
#' plotfit(x,w=max(x),results_random_wrong$model$par,nclass=30,
#'       main="eq 8 likelihood\nadjusted beta")
fitdata=function(x,beta=NULL,beta_eps=-3,bounds=NULL,w=Inf,
		weps=5,wrong=FALSE,debug=FALSE,method="L-BFGS-B",hessian=FALSE)
{
	if(w<max(x))
		x=x[x<=w]
	n=length(x)
    scale=max(x)
    x=x/scale
	w=w/scale
	if(is.null(beta))
		beta=log(sqrt(mean((x)^2)))
    else
	{
	   if(length(beta)!=1)stop("beta not of length 1\n")
	   beta=beta-log(scale)	
	}
	if(!is.null(beta_eps) && length(beta_eps)!=1)stop("beta_eps not of length 1\n")	
	if(is.null(beta_eps))
	{
		if(is.null(bounds))bounds=matrix(c(beta-abs(beta)/2,beta+abs(beta)*2),nrow=1,ncol=2)
		# fit fixed scale model
		model=optimize(flnl,lower=bounds[1,1],upper=bounds[1,2],
				x=x,w=w,weps=weps,debug=debug)
        if(abs(model$par-bounds[1,1])<1e-10) warning("beta at lower bound =",bounds[1,1],"\n")
		if(abs(model$par-bounds[1,2])<1e-10) warning("beta at upper bound =",bounds[1,2],"\n")
		model$par=model$minimum+log(scale)
		model$objective=model$objective + n*log(scale)

	}else
	{
#   fit random scale model
		if(!is.null(bounds) | method=="L-BFGS-B")
		{
			if(is.null(bounds))
				bounds=matrix(c(-3,2,-10,2),nrow=2,byrow=TRUE)
			else
			    bounds[1,]=bounds[1,]-log(scale)
			model=optim(par=c(beta,beta_eps),flnl,x=x,w=w,
					lower=bounds[,1],upper=bounds[,2],weps=weps,
					wrong=wrong,debug=debug,method="L-BFGS-B",hessian=hessian)
		}
	    else
			model=optim(par=c(beta,beta_eps),flnl,x=x,w=w,weps=weps,
					wrong=wrong,debug=debug,method=method)
		if(abs(model$par[1]-bounds[1,1])<1e-10) warning("beta at lower bound =",bounds[1,1],"\n")
		if(abs(model$par[1]-bounds[1,2])<1e-10) warning("beta at upper bound =",bounds[1,2],"\n")	
		model$par[1]=model$par[1]+log(scale)
		model$value=model$value + n*log(scale)
	}
	model$x=x*scale
	model$w=w*scale
	return(list(model=model))
}
# Negative log-likelihood function 
flnl=function(par,x,w,weps=5,wrong=FALSE,debug)
{
	lnl=0
	if(debug)
	{
		cat("beta = ",par[1],"\n")
		cat("beta_eps = ",par[2],"\n")
	}	
	for(i in 1:length(x))
	{
		if(length(par)>1)
			if(wrong)
				lnl=lnl-log(integrate(fx.eps,-weps,weps,
								x=x[i],w=w,par=par)$value+1e-10)
			else
				lnl=lnl-log(avg_gx(x[i],par)/
								avg_mu(par,w,weps))
		else
		{
			sigma=exp(par[1])
			mu=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)
			lnl=lnl-log(gx(x[i],sigma=sigma)+1e-10)+log(mu+1e-10)
		}
	}
	if(debug)cat("lnl=",lnl,"\n")
	lnl
}


