#' Random Scale Detection Function Fitting with admb
#' 
#' Fits a half-normal detection function with random scale for line
#' transect sampling data. Uses one of 3 likelihoods. g,f1,f2,fixed
#' 
#' Beta[1] and bounds for beta[1] are set based on the measured distance scale and
#' they are adjusted internally by scaling distances such that max(x)=1.
#' Likewise the likelihood that is shown when debug=T is for scaled distances
#' but the reported value for the output is adjusted back to the original scale
#' 
#' @export 
#' @param x vector of distances or dataframe containing observed distances (distance) and other covariates
#' @param w half-width of strip; if infinite w routine sets w to 2*max(x)
#' @param formula formula for scale function 
#' @param beta starting values for beta
#' @param sigma starting value for log sigma
#' @param likelihood character string "g","f1","f2","fixed"; see likelihoods.pdf in the package directory
#' @param extra.args for admb run
#' @param verbose for compile and run
#' @param nsteps adromb integration argument; default 8.
#' @param keep if TRUE, uses existing tpl and exe
#' @param debug if TRUE output parameter and -lnl values during iterations
#' @author Jeff Laake
#' @examples 
#' # random effect example in paper
#'dev.new()
#'par(mfrow=c(1,3))
#'set.seed(123)
#'# simulate data
#'x=simdata(n=500,w=Inf,beta=2,beta_eps=-.5)
#'# fit data with g likelihood eq(6) using R code 
#'results_random=fitdata(x,w=Inf)
#'plotfit(x,w=max(x),results_random$model$par,nclass=30,
#'		main="R code\neq 6 likelihood",adjust=FALSE)
#'# fit data with g likelihood eq (6) using ADMB
#'glike=fitadmb(x,w=Inf,likelihood="g",verbose=FALSE)
#'plotfit(x,w=Inf, glike$coefficients[1:2],nclass=30,
#'		main="ADMB code\neq 6 likelihood",adjust=FALSE)
#'# fit data with f likelihood eq (9) using ADMB
#'f2like=fitadmb(x,w=Inf,likelihood="f2",verbose=FALSE)
#'plotfit(x,w=Inf,f2like$coefficients[1:2],nclass=30,
#'		main="ADMB code\neq 9 likelihood")
#' #results in table 1
#' # R code
#' sprintf("%7.2f",-results_random$model$val)
#' sprintf("%5.3f",results_random$model$par[1])
#' sprintf("%5.3f",results_random$model$par[2])
#' # Admb code g likelihood
#' sprintf("%7.2f",glike$loglik)
#' sprintf("%5.3f",glike$coeflist[[1]])
#' sprintf("%5.3f",glike$coeflist[[2]])
#' # Admb code with f likelihood (f2)
#' sprintf("%7.2f",f2like$loglik)
#' sprintf("%5.3f",f2like$coeflist[[1]])
#' sprintf("%5.3f",f2like$coeflist[[1]]-exp(2*f2like$coeflist[[2]]))
#' sprintf("%5.3f",f2like$coeflist[[2]])
#'# mixed efffect example in paper
#'dev.new()
#'par(mfrow=c(1,2))
#'# simulate data
#'x1=simdata(n=2000,w=50,beta_eps=-.5,beta=2,
#'		fixed=FALSE,reject=TRUE)
#'x2=simdata(n=1000,w=50,beta_eps=-.5,beta=1,
#'		fixed=FALSE,reject=TRUE)
#'df=data.frame(covariate=c(rep(0,length(x1)),
#'				rep(1,length(x2))),distance=c(x1,x2))
#'# fit data with covariate
#'fwlike=fitadmb(df,w=50,formula=~covariate,
#'		likelihood="f2",verbose=FALSE)
#'param=fwlike$coefficients[1:3]
#'# plot and get estimates of abundance in covered area and its std error 
#'Nhatwcov=plotfit(df$distance,w=50,par=param,nclass=30,
#'		dm=model.matrix(~covariate,df),
#'		main="With covariate")
#'Nhatwcov.se=compute_Nhat.se(par=param,fwlike$vcov[1:3,1:3],
#'		x=df,w=50,dm=model.matrix(~covariate,df))
#'# fit data without covariate
#'flike=fitadmb(df,w=50,formula=~1,
#'		likelihood="f2",verbose=FALSE)
#'param=flike$coefficients[1:2]
#'# plot and get estimates of abundance in covered area and its std error 
#'Nhatwocov=plotfit(df$distance,w=50,par=param,nclass=30,
#'		main="Without covariate")
#'Nhatwocov.se=compute_Nhat.se(par=param,flike$vcov[1:2,1:2],
#'		x=df,w=50,dm=model.matrix(~1,df))
#'# The code to show delta AIC, abundance and std errors and sigma estimates is
#'round(-2*flike$loglik+2*2-(-2*fwlike$loglik+2*3),2)
#'round(Nhatwcov,0)
#'round(Nhatwcov.se,1)
#'round(Nhatwocov,0)
#'round(Nhatwocov.se,1)
#'round(exp(fwlike$coefficients[3]),2)
#'round(exp(flike$coefficients[2]),2)
#' # plots in figure 3 and results in paper
#'dev.new()
#'par(mfrow=c(2,2))
#'param=fwlike$coefficients[1:3]
#'Nhatwcov0=plotfit(df$distance[df$covariate==0],w=50,par=param,
#'		nclass=30,dm=model.matrix(~covariate,df[df$covariate==0,]),
#'		main="Model with covariate\ncovariate value=0")
#'Nhatwcov0.se=compute_Nhat.se(par=param,fwlike$vcov[1:3,1:3],
#'		x=df[df$covariate==0,],w=50,
#'		dm=model.matrix(~covariate,df[df$covariate==0,]))
#'Nhatwcov1=plotfit(df$distance[df$covariate==1],w=50,par=param,
#'		nclass=30,dm=model.matrix(~covariate,df[df$covariate==1,]),
#'		main="Model with covariate\ncovariate value=1")
#'Nhatwcov1.se=compute_Nhat.se(par=param,fwlike$vcov[1:3,1:3],
#'		x=df[df$covariate==1,],w=50,
#'		dm=model.matrix(~covariate,df[df$covariate==1,]))
#'param=flike$coefficients[1:2]
#'Nhatwocov0=plotfit(df$distance[df$covariate==0],w=50,par=param,
#'		nclass=30, main="Model without covariate\ncovariate value=0")
#'Nhatwocov0.se=compute_Nhat.se(par=param,flike$vcov[1:2,1:2],
#'		x=df[df$covariate==0,],w=50,
#'		dm=model.matrix(~1,df[df$covariate==0,]))
#'Nhatwocov1=plotfit(df$distance[df$covariate==1],w=50,par=param,
#'		nclass=30,main="Model without covariate\ncovariate value=1")
#'Nhatwocov1.se=compute_Nhat.se(par=param,flike$vcov[1:2,1:2],
#'		x=df[df$covariate==1,],w=50,
#'		dm=model.matrix(~1,df[df$covariate==1,]))
#'round(Nhatwcov0,0)
#'round(Nhatwcov0.se,1)
#'round(Nhatwcov1,1)
#'round(Nhatwcov1.se,1)
#'round(Nhatwocov0,0)
#'round(Nhatwocov0.se,1)
#'round(Nhatwocov1,1)
#'round(Nhatwocov1.se,1)
fitadmb=function(x,w=Inf,formula=~1,beta=NULL,sigma=0,likelihood="f2",
		extra.args="-gh 10",verbose=TRUE,nsteps=8,keep=FALSE,debug=FALSE)
{
	sdir=system.file(package="RandomScale")
	if(!likelihood%in%c("g","f1","f2","fixed"))stop("incorrect likelihood string")
	if(!file.exists("df1b2gh.cpp"))
		file.copy(file.path(sdir,"df1b2gh.cpp"),file.path(getwd(),"df1b2gh.cpp"),overwrite=TRUE)
	if(!file.exists("minfil.cpp"))
		file.copy(file.path(sdir,"minfil.cpp"),file.path(getwd(),"minfil.cpp"),overwrite=TRUE)
	if(!file.exists("xmodelm5.cpp"))
		file.copy(file.path(sdir,"xmodelm5.cpp"),file.path(getwd(),"xmodelm5.cpp"),overwrite=TRUE)
    if(likelihood=="fixed")
		tpl="distcov"
    else
	{
		if(formula!=~1)
			tpl="mixed_hnre_f2"
		else
			tpl=paste("hnre_",likelihood,sep="")
	}
	if(!keep)
	{
		if(file.exists(paste(tpl,".exe",sep=""))) unlink(paste(tpl,".exe",sep=""))
		if(file.exists(paste(tpl,".tpl",sep=""))) unlink(paste(tpl,".tpl",sep=""))
	}
	clean_admb(tpl)
#####################	
#   Create DAT file
#####################
	con=file(paste(tpl,".dat",sep=""),open="wt")
    dm=NULL
	if(!is.vector(x))
	{
		if(is.null(x$distance))
			stop("missing distance field")
		else
		{
			if(w<max(x$distance))
				x=x[x$distance<=w,]
			dm=model.matrix(formula,data=x)
			x=x$distance
		}
	} else
	{
		if(w<max(x))
			x=x[x<=w]
	}		
	n=length(x)
	write(n,con,append=FALSE)
	if(w==Inf)w=2*max(x)
	scale=max(x)
	write(w/scale,con,append=TRUE)
    # these bounds are hard-coded in tpl files
	bounds=matrix(c(-3,2,-10,1),nrow=2,byrow=TRUE)
	write(as.numeric(debug),con,append=TRUE)
	write(x/scale,con,ncolumns=1,append=TRUE)
	if(likelihood=="fixed")
	{
		writeLines("1",con)
		writeLines("1",con)
		writeLines("2",con)
		write(nsteps,con,append=TRUE)
	}
	if((formula!=~1) | (likelihood=="fixed"))
	{
		if(is.null(dm))dm=matrix(1,ncol=1,nrow=length(x))
		write(ncol(dm),con,append=TRUE)
		write(t(dm),con,ncolumns=ncol(dm),append=TRUE)
		npar=ncol(dm)
	}else
		npar=1
	close(con)
#####################	
#   Create PIN file
#####################
	con=file(paste(tpl,".pin",sep=""),open="wt")
	if(is.null(beta))
	{
		beta=log(sqrt(mean((x/scale)^2)))
		if(formula!=~1)beta=c(beta,rep(0,ncol(dm)))
	} else
	{
		if(length(beta)!=(npar-1))stop("beta not of proper length = ",npar,"\n")
		beta[1]=beta[1]-log(scale)
	}
	write(beta,ncolumns=length(beta),con)
	if(length(sigma)!=1)stop("sigma not of length 1\n")
	write(sigma,con,append=TRUE)
	if(formula!=~1)
	   write(rep(0,2*length(x)),ncolumns=2*length(x),con,append=TRUE)
    else
	{
		if(likelihood=="f1")
			write(rep(0,length(x)),ncolumns=length(x),con,append=TRUE)
		else
		{
			if(likelihood!="fixed")
				write(rep(0,length(x)+1),ncolumns=length(x)+1,con,append=TRUE)
		}
	}
	close(con)
##################
# Compile TPL file
##################
    if(!file.exists(paste(tpl,".tpl",sep="")))
	    file.copy(file.path(sdir,paste(tpl,".tpl",sep="")),file.path(getwd(),paste(tpl,".tpl",sep="")),overwrite=TRUE)
    if(!file.exists(paste(tpl,".exe",sep="")))
	{
		if(file.exists(file.path(sdir,paste(tpl,".exe",sep=""))))
		   file.copy(file.path(sdir,paste(tpl,".exe",sep="")),file.path(getwd(),paste(tpl,".exe",sep="")),overwrite=TRUE)
		else
		{
		    if(likelihood=="fixed")
			   compile_admb(tpl,verbose=verbose)	
		    else
		       compile_admb(tpl,re=TRUE,verbose=verbose)
		    if(!file.exists(paste(tpl,".exe",sep="")))stop("problem with tpl compile")
	   }
	}
##############################
# Run program and read results
##############################
    run_admb(tpl,extra.args=extra.args,verbose=verbose)
	results=read_admb(tpl,checkterm=FALSE)
	if(formula==~1)
	{
		if(abs(results$coeflist$beta[1]-bounds[1,1])<1e-10) warning("beta at lower bound =",bounds[1,1],"\n")
		if(abs(results$coeflist$beta[1]-bounds[1,2])<1e-10) warning("beta at upper bound =",bounds[1,2],"\n")	
	}
	if(likelihood!="fixed")
	{
		if(abs(results$coeflist$sigeps[1]-bounds[2,1])<1e-10) warning("sigeps at lower bound =",bounds[2,1],"\n")
	    if(abs(results$coeflist$sigeps[1]-bounds[2,2])<1e-10) warning("sigeps at upper bound =",bounds[2,2],"\n")	
	}
	results$coeflist$beta[1]=results$coeflist$beta[1]+log(scale)
	results$coefficients[1]=results$coefficients[1]+log(scale)
	results$loglik=results$loglik - n*log(scale)
	if(likelihood=="fixed")
	{
		cnames=paste("scale:",colnames(dm),sep="")
		names(results$coefficients)=cnames
		if(!is.null(results$vcov))
		{
			rownames(results$vcov)=cnames
			colnames(results$vcov)=cnames
		}
	}
	return(results)
}
