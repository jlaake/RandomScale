Random Scale
========

The `RandomScale` package provides R code to fit a random scale half-normal detection function to line transect sampling data. The package also 
contains R code to use [AD Model Builder](http://admb-project.org) with the R package [R2admb](https://github.com/bbolker/R2admb) to fit mixed-effects models with known covariate values and a random component. 
The ADMB TPL files are contained in the [Inst directory](https://github.com/jlaake/RandomScale/tree/master/RandomScale/inst).
Example AMDB DAT files are contained in the [Inst/Data directory](https://github.com/jlaake/RandomScale/tree/master/RandomScale/inst/Data). 
Version 11 of ADMB is required for the TPL files as distributed but you can modify them by replacing instances of PI with the value 3.141592654.

To install ADMB, you need to install a C++ compiler as well as ADMB. We suggest that
you install gcc. To use the built-in links, [ADMB-11](http://admb-project.googlecode.com/files/admb-11-linux-gcc4.6-32bit.zip) should be installed 
to c:\admb and [gcc compiler](http://www.admb-project.org/tools/gcc/gcc452-win32.zip/at_download/file) should be installed to c:\mingw. 
To use different locations you'll need to change function prepare_admb which is run when the package is attached.

For other operating systems see (http://www.admb-project.org/downloads) and
(http://www.admb-project.org/tools/gcc/). Note that prepare_admb() only works for Windows.

Download [Windows package binary](https://docs.google.com/folder/d/0B77g1ScdUwVeOVJNUVVGS0YtWE0/edit). From link, browse to RandomScale and then click on
the version of the package you want. You should see a listing of the package contents as files.  Select File/Download. 
To install in R, from the R menu, use Packages\Install from Local Zip file and browse to location of downloaded zip. 
