#' Random Scale Detection Functions for Line Transect Data
#' 
#' The RandomScale package provides R and AD Model Builder(\url{http://admb-project.org})code to fit random effect and mixed effect models for the scale
#' of a half-normal detection function to line transect sampling data. 
#' 
#' The ADMB TPL files are contained in the Inst directory(\url{https://github.com/jlaake/RandomScale/tree/master/RandomScale/inst}).
#' Example AMDB DAT files are contained in the Data directory(\url{https://github.com/jlaake/RandomScale/tree/master/RandomScale/Data}).
#' Version 11 of ADMB is required for the TPL files as distributed but you can modify them by replacing instances of PI with the 
#' value 3.141592654. Executable versions are also provided if you download the binary version for Windows at \url{https://drive.google.com/folder/d/0B77g1ScdUwVeOVJNUVVGS0YtWE0/edit?ddrp=1#}
#' In that case you do not need ADMB but the R2admb package is required because it is used to extract results from the output file. The examples
#' also use the mrds package which is available on CRAN. Both R2admb and mrds should be installed separately.
#'  
#' If you wish to to install ADMB, you need to install a C++ compiler as well as ADMB. We suggest that
#' you install gcc. To use the built-in links in the function prepare_admb, admb should be installed to c:\\admb and the gcc
#' compiler to c:\\mingw. To use different locations you'll need to change function prepare_admb. Unless you use existing 
#' exe's in the package, you'll need to run prepare_admb for each R session. 
#' 
#' For Windows install ADMB-11 (\url{http://admb-project.googlecode.com/files/admb-11-linux-gcc4.6-32bit.zip}) and
#' gcc(\url{http://www.admb-project.org/tools/gcc/gcc452-win32.zip/at_download/file})
#' 
#' For other operating systems see (\url{http://www.admb-project.org/downloads}) and
#' gcc(\url{http://www.admb-project.org/tools/gcc/}). Note that prepare_admb() only works for Windows.
#'  
#' For more information see \url{https://github.com/downloads/jlaake/ADMB-Examples/distance_random_effect.pdf}. Also the file likelihoods.pdf
#' in the RandomScale package directory contains a description of the alternative likelihoods including an incorrect likelihood f1 which was
#' removed from the paper but remains in the code.  
#' @name RandomScale-package
#' @aliases RandomScale
#' @author Jeff Laake
#' @references Fournier, D.A., H.J. Skaug, J. Ancheta, J. Ianelli, A. Magnusson, M.N. Maunder, A. Nielsen, and J. Sibert. 2012. AD Model Builder: using automatic differentiation for statistical inference of highly parameterized complex nonlinear models. Optim. Methods Softw. 27:233-249.
#' @examples
#' # the following runs the examples in Oedekoven et al; 
#' # note that you need to change nreps to 500 in simulation to get the same example
#' example(fitadmb)
#' example(hp)
#' example(simulation)
#' 
NULL


#' Harbor porpoise vessel survey data
#' 
#' Line transect sampling data for habor porpoise from a vessel survey in the inland Washington waters
#' 
#' The data include perpendicular distances to observed harbor porpoise groups and the size of the group.
#' It also includes observer (always 1), detected (always 1) and object (sequential number) that are needed 
#' by mrds package.
#' 
#' @name hp
#' @docType data
#' @format A data frame with 477 observations on the following 5 variables.
#' \describe{ \item{object}{sequential object #}
#'            \item{distance}{perpendicular distance to observed harbor porpoise group}
#'            \item{size}{size of harbor porpoise group}
#'            \item{detected}{always 1}
#'            \item{observer}{always 1}}
#' 
#' @keywords datasets
#' @examples
#'
#' dev.new()
#'par(mfrow=c(1,2))
#'data(hp)
#'# hp=read.delim("hp.txt",header=TRUE)
#'# fit half-normal detection function with group size as a covariate
#'modmcds=ddf(dsmodel=~cds(key="hn",formula=~size),data=hp,method="ds")
#'# show summary of results
#'xx=summary(modmcds)
#'xx
#'# fit mixed effect model with eq (9) and use same truncation as mrds analysis
#'modmixed=fitadmb(hp,formula=~size,w=443.1635)
#'summary(modmixed)
#'# plot results and compute abundances in covered area
#'Nhat.mcds=xx$Nhat
#'Nhat.se.mcds=xx$Nhat.se
#'par(mfrow=c(2,1))
#'plotfit(hp$distance,w=443.1635,dm=model.matrix(~size,hp[,"size",drop=FALSE]),
#'        par=modmixed$coeff[1:3])
#'plot(modmcds,new=FALSE,showpoints=FALSE,pl.den=0,main="")
#'Nhat.mixed=compute_Nhat(par=modmixed$coeff[1:3],x=hp$distance,w=443.1635,
#'                 dm=model.matrix(~size,hp[,"size",drop=FALSE]))
#'Nhat.se.mixed=compute_Nhat.se(par=modmixed$coeff[1:3],vcov=modmixed$vcov[1:3,1:3],
#'                 x=hp$distance,w=443.1635,dm=model.matrix(~size,hp[,"size",drop=FALSE]))
#'Nhat.mcds
#'Nhat.se.mcds
#'Nhat.mixed
#'Nhat.se.mixed

NULL

#' Simulation in Oedekoven et al 
#' 
#' The following is the code used to generate the simulation results in the paper.  To produce the
#' same results set nreps=5 to nreps=500. It will likely take about 10 hours to run all of the simulations
#' depending on your computer. We have set nreps to 5 to reduce time for computation with the example.  
#' 
#' @name simulation
#' @examples
#' #define function for simulations
#'tsims=function(N=4000,df=3,w=40,nreps=100,plot=TRUE,hz=FALSE,
#'              tdf=TRUE,hzscale=1,hzpow=2,sigma=1)
#'{
#'maxw=w
#'func_t=function(x,df) return(dt(x,df)/dt(0,df))
#'TrueN=N/w
#'fN=matrix(NA,ncol=2,nrow=nreps)
#'gN=matrix(NA,ncol=2,nrow=nreps)
#'hrN=matrix(NA,ncol=2,nrow=nreps)
#'fAIC=matrix(NA,ncol=2,nrow=nreps)
#'gAIC=matrix(NA,ncol=2,nrow=nreps)
#'hrAIC=matrix(NA,ncol=2,nrow=nreps)
#'nobs=matrix(NA,ncol=2,nrow=nreps)
#'w=matrix(NA,ncol=2,nrow=nreps)
#'for (i in 1:nreps)
#'{
#'u=runif(N,0,maxw)
#'if(hz)
#'  seen=mrds:::keyfct.hz(u,hzscale,hzpow)>=runif(N,0,1)
#'else
#' if(tdf)
#'     seen=func_t(u,df=df)>=runif(N,0,1)
#'  else
#'     seen=mrds:::keyfct.hn(u,key.scale=exp(hzscale+rnorm(N,0,sigma)))>=runif(N,0,1)
#'xx=u[seen]
#'n=length(xx)
#'data=data.frame(distance=xx,observer=rep(1,n),detected=rep(1,n),object=1:n)
#'for(j in 1:2)
#'{
#'   nobs[i,j]=n
#'   if(j==1)
#'     w[i,j]=max(data$distance)
#'   else
#'	 w[i,j]=min(maxw,2*max(data$distance))
#'# code for hazard rate
#'   if(hz)
#'      mod.hr=ddf(dsmodel=~cds(key="hr"),method="ds",data=data,meta.data=list(width=w[i,j]),
#'              control=list(initial=list(scale=log(hzscale),shape=log(hzpow))))
#'   else
#'      mod.hr=ddf(dsmodel=~cds(key="hr"),method="ds",data=data,
#'                  meta.data=list(width=w[i,j]))
#'# code for g likelihood
#'   mod.admb=fitadmb(data, w = w[i,j], formula = ~1,likelihood="g",keep=TRUE)
#'# code for f2 likelihood
#'   mod.admbf=fitadmb(data, w = w[i,j], formula = ~1,likelihood="f2",keep=TRUE)
#'   hrN[i,j]=mod.hr$Nhat/w[i,j]
#'   hrAIC[i,j]=mod.hr$criterion
#'   gN[i,j]=compute_Nhat(mod.admb$coefficients[1:2],x=data$distance,
#'                               w=w[i,j],adjust=F)/w[i,j]
#'   gAIC[i,j]=-2*mod.admb$loglik+4
#'   fN[i,j]=compute_Nhat(mod.admbf$coefficients[1:2],x=data$distance,
#'                               w=w[i,j],adjust=T)/w[i,j]
#'   fAIC[i,j]=-2*mod.admbf$loglik+4
#'   nclass=ceiling(sqrt(n)*2)
#'   if(plot)
#'   {
#'   par(mfrow=c(3,1))
#'   param=c(mod.admb$coefficients[1],mod.admb$coefficients[2])
#'   plotfit(data$distance,w=w[i,j],param,nclass=nclass,adjust=F)
#'   if(hz)
#'      lines((0:(ceiling(w[i,j]*100)))/100,
#'            mrds:::keyfct.hz((0:(ceiling(w[i,j]*100)))/100,hzscale,hzpow),lty=2)
#'   else
#'      lines((0:(ceiling(w[i,j]*100)))/100,
#'            func_t((0:(ceiling(w[i,j]*100)))/100,df=df),lty=2) 
#'   param=c(mod.admbf$coefficients[1],mod.admbf$coefficients[2])
#'   plotfit(data$distance,w=w[i,j],c(param[1]-exp(2*param[2]),param[2]),
#'            nclass=nclass,adjust=F)
#'   if(hz)
#'      lines((0:(ceiling(w[i,j]*100)))/100,
#'            mrds:::keyfct.hz((0:(ceiling(w[i,j]*100)))/100,hzscale,hzpow),lty=2)
#'   else
#'      lines((0:(ceiling(w[i,j]*100)))/100,
#'            func_t((0:(ceiling(w[i,j]*100)))/100,df=df),lty=2)
#'   plot(mod.hr,nc=nclass,showpoints=FALSE)
#'   if(hz)
#'      lines((0:(ceiling(w[i,j]*100)))/100,
#'            mrds:::keyfct.hz((0:(ceiling(w[i,j]*100)))/100,hzscale,hzpow),lty=2)
#'   else
#'      lines((0:(ceiling(w[i,j]*100)))/100,
#'            func_t((0:(ceiling(w[i,j]*100)))/100,df=df),lty=2)
#'   }
#'}
#'}
#' return(list(TrueN=TrueN,fN=fN,gN=gN,hrN=hrN,fAIC=fAIC,
#'              gAIC=gAIC,hrAIC=hrAIC,nobs=nobs,w=w))
#'}
#'# define function to summarize results
#'sim_summary=function(z)
#'{ 
#'   res=NULL
#'   TrueN=z$TrueN
#'   for (i in 1:2)
#'   {
#'   x=lapply(z[2:9],function(x)return(x[,i]))
#'   AvgN=with(x,c(fN[fAIC<hrAIC],hrN[fAIC>=hrAIC]))
#'   prb_avg=100*(AvgN-TrueN)/TrueN
#'   res=rbind(res,c(w=mean(x$w),nobs=mean(x$nobs),prb_f=mean(100*(x$fN-TrueN)/TrueN),
#'       prb_hr=mean(100*(x$hrN-TrueN)/TrueN),prb_avg=mean(prb_avg),
#'       se_prb=sqrt(var(prb_avg)/length(x$w)),
#'   prop_f_better_than_hr=with(x,mean(fAIC<hrAIC)),
#'   prop_f_better_than_g=with(x,mean(fAIC<gAIC)),
#'   prop_g_better_than_f=with(x,mean(fAIC>gAIC)),
#'   rmse_hr =with(x,100*sqrt((TrueN-mean(x$hrN))^2+var(x$hrN)))/TrueN,
#'   rmse_f = with(x,100*sqrt((TrueN-mean(x$fN))^2+var(x$fN)))/TrueN,
#'   rmse_avg =with(x,100*sqrt((TrueN-mean(AvgN))^2+var(AvgN)))/TrueN))
#'   }
#'   return(res)
#'}
#' # perform simulations
#' nreps=5
#' set.seed(93851)
#'t.df3.n180=tsims(plot=F,nreps=nreps)
#'t.df5.n180=tsims(df=5,plot=F,nreps=nreps)
#'t.df10.n180=tsims(df=10,plot=F,nreps=nreps)
#'t.df3.n90=tsims(N=2000,plot=F,nreps=nreps)
#'t.df5.n90=tsims(N=2000,df=5,plot=F,nreps=nreps)
#'t.df10.n90=tsims(N=2000,df=10,plot=F,nreps=nreps)
#'hn.n180=tsims(N=8000,td=FALSE,hzscale=-.5,sigma=0.5,plot=F,nreps=nreps)
#'hn.n90=tsims(N=4000,td=FALSE,hzscale=-.5,sigma=0.5,plot=F,nreps=nreps)
#'hr.n180=tsims(N=7000,hz=TRUE,hzscale=.7,hzpow=2.5,plot=F,nreps=nreps)
#'hr.n90=tsims(N=3500,hz=TRUE,hzscale=.7,hzpow=2.5,plot=F,nreps=nreps)
#'# pull together results into a table which can be used by xtable for TeX file
#'sim_results=rbind(sim_summary(t.df3.n180),sim_summary(t.df5.n180),
#'                  sim_summary(t.df10.n180),sim_summary(t.df3.n90),
#'                  sim_summary(t.df5.n90),sim_summary(t.df10.n90),
#'                  sim_summary(hn.n180),sim_summary(hn.n90),sim_summary(hr.n180),
#'                  sim_summary(hr.n90))
#'sim_results=cbind(c("t(df=3)","","t(df=5)","","t(df=10)","","t(df=3)","","t(df=5)","",
#'                  "t(df=10)","","hn","","hn","","hz","","hz",""),
#'                   apply(sim_results,2,function(x) sprintf("%0.2f",x)))
#'colnames(sim_results)[1]="Function"
#'colnames(sim_results)[2]="$\\bar{w}\\;\\;$"
#'colnames(sim_results)[3]="$\\bar{n}\\;\\;\\;$"
#'colnames(sim_results)[4]="$PRB_F$"
#'colnames(sim_results)[5]="$PRB_{HR}$"
#'colnames(sim_results)[6]="$PRB_{AVG}$"
#'colnames(sim_results)[7]="$se(PRB_{AVG})$"
#'colnames(sim_results)[8]="$AIC_F<AIC_{HR}$"
#'colnames(sim_results)[9]="$AIC_F<AIC_G$"
#'colnames(sim_results)[10]="$AIC_G<AIC_F$"
#'colnames(sim_results)[11]="$RMSE_F$"
#'colnames(sim_results)[12]="$RMSE_{HR}$"
#'colnames(sim_results)[13]="$RMSE_{AVG}$"
#' #commented out code that produces table in paper
#' #library(xtable)
#' #print(xtable(sim_results,caption="Percent relative bias (PRB) and root mean square error (RMSE) as proportion of true abundance for random scale half-normal and hazard rate detection function models for distances generated from t-distribution, random scale half-normal and hazard rate detection functions. Each value is the summary for 100 replicate simulations. The subscripts F, G and HR refer to the likelihoods eq (6), eq (9) and the hazard rate. AVG subscript represents the values in which the estimate was generated from the model that had the lowest AIC for each replicate.",
#' #label="simresults",align=c("c","c",rep("r",12))),
#' #caption.placement="top",latex.environments="center",size="scriptsize",include.rownames=FALSE,
#' #sanitize.colnames.function = function(x){x}, sanitize.text.function=function(x){x})
NULL

