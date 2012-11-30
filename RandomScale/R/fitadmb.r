#' Random Scale Detection Function Fitting with admb
#' 
#' Fits a half-normal detection function with random scale for line
#' transect sampling data. Uses one of 3 likelihoods. g,f1,f2
#' 
#' @export 
#' @param x vector of distances or dataframe containing observed distances (distance) and other covariates
#' @param w half-width of strip 
#' @param formula formula for scale function 
#' @param likelihood character string "g","f1","f2"
#' @param extra.args for admb run
#' @param verbose for compile and run
#' @author Jeff Laake
fitadmb=function(x,w,formula=~1,likelihood="f2",extra.args="-est -gh 10",verbose=FALSE)
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
	write(w,con,append=TRUE)
	write(x,con,ncol=1,append=TRUE)
	write(ncol(dm),con,append=TRUE)
	write(t(dm),con,ncol=ncol(dm),append=TRUE)
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
