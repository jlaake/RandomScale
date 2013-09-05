### R code from vignette source 'RandomScalePaper.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: RandomScalePaper.Rnw:19-20
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: RandomScalePaper.Rnw:69-75
###################################################
library(RandomScale)
library(R2admb)
#**************************************************
# NOTE: ADMB and MINGW C++ must be installed and use
# these directories for the code to work or you
# you need to change the path specifications below
#**************************************************
Sys.setenv(PATH = paste("c:/admb/bin;c:admb/utilities;c:/MinGW/bin;", 
            Sys.getenv("PATH"), sep = ";"))
Sys.setenv(ADMB_HOME = "c:/admb")
set.seed(123)


###################################################
### code chunk number 3: RandomScalePaper.Rnw:528-530
###################################################
png("RcodeFitting.png")
par(mfrow=c(1,3)) 


###################################################
### code chunk number 4: RandomScalePaper.Rnw:533-545
###################################################
# simulate data
x=simdata(n=500,w=Inf,beta=2,beta_eps=-.5)
# fit data with g likelihood
results_random=fitdata(x,w=Inf)
plotfit(x,w=max(x),results_random$model$par,nclass=30,
        main="eq 6 likelihood")
# fit data with incorrect likelihood
results_random_wrong=fitdata(x,w=Inf,wrong=TRUE)
plotfit(x,w=max(x),results_random_wrong$model$par,nclass=30,
        main="eq 11 likelihood",adjust=FALSE)
plotfit(x,w=max(x),results_random_wrong$model$par,nclass=30,
       main="eq 11 likelihood\nadjusted beta")


###################################################
### code chunk number 5: RandomScalePaper.Rnw:548-549
###################################################
dev.off()


###################################################
### code chunk number 6: RandomScalePaper.Rnw:565-567
###################################################
png("ADMBcodeFitting.png")
par(mfrow=c(1,3)) 


###################################################
### code chunk number 7: RandomScalePaper.Rnw:570-579
###################################################
glike=fitadmb(x,w=Inf,likelihood="g",verbose=FALSE)
plotfit(x,w=Inf, glike$coefficients[1:2],nclass=30,
                 main="eq 6 likelihood",adjust=FALSE)
f2like=fitadmb(x,w=Inf,likelihood="f2",verbose=FALSE)
plotfit(x,w=Inf,f2like$coefficients[1:2],nclass=30,
                 main="eq 9 likelihood\nadjusted beta")
f1like=fitadmb(x,w=Inf,likelihood="f1",verbose=FALSE)
plotfit(x,w=Inf,f1like$coefficients[1:2],nclass=30,
                 main="eq 11 likelihood\nadjusted beta")


###################################################
### code chunk number 8: RandomScalePaper.Rnw:583-584
###################################################
dev.off()


###################################################
### code chunk number 9: RandomScalePaper.Rnw:609-611
###################################################
png("ADMBWithandWithoutCovariate.png")
par(mfrow=c(1,2))


###################################################
### code chunk number 10: RandomScalePaper.Rnw:614-638
###################################################
# simulate data
x1=simdata(n=2000,w=50,beta_eps=-.5,beta=2,
          fixed=FALSE,reject=TRUE)
x2=simdata(n=1000,w=50,beta_eps=-.5,beta=1,
          fixed=FALSE,reject=TRUE)
df=data.frame(covariate=c(rep(0,length(x1)),
          rep(1,length(x2))),distance=c(x1,x2))
# fit data with covariate
fwlike=fitadmb(df,w=50,formula=~covariate,
               likelihood="f2",verbose=FALSE)
param=fwlike$coefficients[1:3]
Nhatwcov=plotfit(df$distance,w=50,par=param,nclass=30,
                 dm=model.matrix(~covariate,df),
                 main="With covariate")
Nhatwcov.se=compute_Nhat.se(par=param,fwlike$vcov[1:3,1:3],
                 x=df,w=50,dm=model.matrix(~covariate,df))
# fit data without covariate
flike=fitadmb(df,w=50,formula=~1,
              likelihood="f2",verbose=FALSE)
param=flike$coefficients[1:2]
Nhatwocov=plotfit(df$distance,w=50,par=param,nclass=30,
                 main="Without covariate")
Nhatwocov.se=compute_Nhat.se(par=param,flike$vcov[1:2,1:2],
                 x=df,w=50,dm=model.matrix(~1,df))


###################################################
### code chunk number 11: RandomScalePaper.Rnw:642-643
###################################################
dev.off()


###################################################
### code chunk number 12: RandomScalePaper.Rnw:658-673
###################################################
png("ADMBByCovariate.png")
par(mfrow=c(2,2))
param=fwlike$coefficients[1:3]
Nhatwcov0=plotfit(df$distance[df$covariate==0],w=50,par=param,
		nclass=30,dm=model.matrix(~covariate,df[df$covariate==0,]),
		main="With covariate value=0")
Nhatwcov0.se=compute_Nhat.se(par=param,fwlike$vcov[1:3,1:3],
		x=df[df$covariate==0,],w=50,
		dm=model.matrix(~covariate,df[df$covariate==0,]))
Nhatwcov1=plotfit(df$distance[df$covariate==1],w=50,par=param,
		nclass=30,dm=model.matrix(~covariate,df[df$covariate==1,]),
		main="With covariate value=1")
Nhatwcov1.se=compute_Nhat.se(par=param,fwlike$vcov[1:3,1:3],
		x=df[df$covariate==1,],w=50,
		dm=model.matrix(~covariate,df[df$covariate==1,]))
param=flike$coefficients[1:2]
Nhatwocov0=plotfit(df$distance[df$covariate==0],w=50,par=param,
		nclass=30, main="Without covariate value=0")
Nhatwocov0.se=compute_Nhat.se(par=param,flike$vcov[1:2,1:2],
		x=df[df$covariate==0,],w=50,
		dm=model.matrix(~1,df[df$covariate==0,]))
Nhatwocov1=plotfit(df$distance[df$covariate==1],w=50,par=param,
		nclass=30,main="Without covariate value=1")
Nhatwocov1.se=compute_Nhat.se(par=param,flike$vcov[1:2,1:2],
		x=df[df$covariate==1,],w=50,
		dm=model.matrix(~1,df[df$covariate==1,]))
dev.off()


