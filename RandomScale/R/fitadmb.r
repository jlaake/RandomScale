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
#' @param likelihood character string "g","f1","f2","fixed"
#' @param extra.args for admb run
#' @param verbose for compile and run
#' @param nsteps adromb integration argument; default 8.
#' @param keep if TRUE, uses existing tpl and exe
#' @param debug if TRUE output parameter and -lnl values during iterations
#' @author Jeff Laake
#' @examples 
#' #fit simulated data with random scale
#' set.seed(123)
#' x=simdata(n=500,w=Inf,beta_eps=-.5)
#' par(mfrow=c(1,3)) 
#' glike=fitadmb(x,w=Inf,likelihood="g")
#' plotfit(x,w=Inf,c(glike$coeflist[[1]],glike$coeflist[[2]]),nclass=30,
#'		main="eq 4 likelihood")
#' f2like=fitadmb(x,w=Inf,likelihood="f2")
#' param=c(f2like$coeflist[[1]],f2like$coeflist[[2]])
#' plotfit(x,w=Inf,c(param[1]-exp(2*param[2]),param[2]),nclass=30,
#'		main="eq 7 likelihood")
#' f1like=fitadmb(x,w=Inf,likelihood="f1")
#' param=c(f1like$coeflist[[1]],f1like$coeflist[[2]])
#' plotfit(x,w=Inf,c(param[1]-exp(2*param[2]),param[2]),nclass=30,
#'		main="eq 8 likelihood")
#' dev.new()
#' par(mfrow=c(1,2)) 
#' #Mixed effect model
#' x1=simdata(n=2000,w=50,beta_eps=-.5,beta=2,fixed=FALSE,reject=TRUE)
#' x2=simdata(n=1000,w=50,beta_eps=-.5,beta=1,fixed=FALSE,reject=TRUE)
#' df=data.frame(covariate=c(rep(0,length(x1)),rep(1,length(x2))),
#' 		distance=c(x1,x2))
#' fwlike=fitadmb(df,w=50,formula=~covariate,likelihood="f2")
#' param=c(fwlike$coeflist[[1]],fwlike$coeflist[[2]])
#' Nhatwcov=plotfit(df$distance,w=50,
#' 		par=c(param[1]-exp(2*param[3]),param[2],param[3]),
#' 		nclass=30,dm=model.matrix(~covariate,df),main="With covariate")
#' flike=fitadmb(df,w=50,formula=~1,likelihood="f2")
#' param=c(flike$coeflist[[1]],flike$coeflist[[2]])
#' Nhatwocov=plotfit(df$distance,w=50,
#' 		par=c(param[1]-exp(2*param[2]),param[2]),nclass=30,
#' 		main="Without covariate")
#' dev.new()
#' par(mfrow=c(2,2))
#' param=c(fwlike$coeflist[[1]],fwlike$coeflist[[2]])
#' Nhatwcov0=plotfit(df$distance[df$covariate==0],w=50,
#' 		par=c(param[1]-exp(2*param[3]),param[2],param[3]),
#' 		nclass=30,dm=model.matrix(~covariate,df[df$covariate==0,]),
#' 		main="With covariate value=0")
#' Nhatwcov1=plotfit(df$distance[df$covariate==1],w=50,
#' 		par=c(param[1]-exp(2*param[2]),param[2],param[3]),nclass=30,
#' 		dm=model.matrix(~covariate,df[df$covariate==1,]),
#' 		main="With covariate value=1")
#' param=c(flike$coeflist[[1]],flike$coeflist[[2]])
#' Nhatwocov0=plotfit(df$distance[df$covariate==0],w=50,
#' 		par=c(param[1]-exp(2*param[2]),param[2]),nclass=30,
#' 		main="Without covariate value=0")
#' Nhatwocov1=plotfit(df$distance[df$covariate==1],w=50,
#' 		par=c(param[1]-exp(2*param[2]),param[2]),nclass=30,
#' 		main="Without covariate value=1")
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
