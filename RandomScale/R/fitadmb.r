#' Random Scale Detection Function Fitting with admb
#' 
#' Fits a half-normal detection function with random scale for line
#' transect sampling data. Uses one of 3 likelihoods. g,f1,f2
#' 
#' @export 
#' @param x vector of distances or dataframe containing observed distances (distance) and other covariates
#' @param w half-width of strip; if infinite w routine sets w to 2*max(x)
#' @param formula formula for scale function 
#' @param likelihood character string "g","f1","f2"
#' @param extra.args for admb run
#' @param verbose for compile and run
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
fitadmb=function(x,w=Inf,formula=~1,likelihood="f2",extra.args="-est -gh 10",verbose=FALSE)
{
	sdir=system.file(package="RandomScale")
	if(!likelihood%in%c("g","f1","f2"))stop("incorrect likelihood string")
	if(!file.exists("df1b2gh.cpp"))
		file.copy(file.path(sdir,"df1b2gh.cpp"),file.path(getwd(),"df1b2gh.cpp"),overwrite=TRUE)
	if(!file.exists("minfil.cpp"))
		file.copy(file.path(sdir,"minfil.cpp"),file.path(getwd(),"minfil.cpp"),overwrite=TRUE)
	if(!file.exists("xmodelm5.cpp"))
		file.copy(file.path(sdir,"xmodelm5.cpp"),file.path(getwd(),"xmodelm5.cpp"),overwrite=TRUE)
	if(formula!=~1)
		tpl="mixed_hnre_f2"
	else
	    tpl=paste("hnre_",likelihood,sep="")
	con=file(paste(tpl,".dat",sep=""),open="wt")
    dm=NULL
	if(!is.vector(x))
	{
		dm=model.matrix(formula,data=x)
		if(is.null(x$distance))
			stop("missing distance field")
		else
			x=x$distance
	}
	write(length(x),con,append=FALSE)
	if(w==Inf)w=2*max(x)
	write(w,con,append=TRUE)
	write(x,con,ncolumns=1,append=TRUE)
	if(formula!=~1)
	{
		write(ncol(dm),con,append=TRUE)
		write(t(dm),con,ncolumns=ncol(dm),append=TRUE)
	}
	close(con)
	if(!file.exists(paste(tpl,".exe",sep="")))
	{
		clean_admb(tpl)
		if(!file.exists(paste(tpl,".tpl",sep="")))
			file.copy(file.path(sdir,paste(tpl,".tpl",sep="")),file.path(getwd(),paste(tpl,".tpl",sep="")),overwrite=TRUE)
		compile_admb(tpl,re=TRUE,verbose=verbose)	
	}
	if(!file.exists(paste(tpl,".exe",sep="")))stop("problem with tpl compile")		
	run_admb(tpl,extra.args=extra.args,verbose=verbose)
	results=read_admb(tpl,checkterm=FALSE)
	return(results)
}
