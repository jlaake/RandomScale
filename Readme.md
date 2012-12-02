Random Scale
========

The `RandomScale` package provides R code to fit a random scale half-normal detection function to line transect sampling data. The package also 
contains R code to use [AD Model Builder](http://admb-project.org) with the R package [R2admb](https://github.com/bbolker/R2admb) to fit mixed-effects models with known covariate values and a random component. 
The ADMB TPL files are contained in the [Inst directory](https://github.com/jlaake/RandomScale/tree/master/RandomScale/inst).
Example AMDB DAT files are contained in the [Data directory](https://github.com/jlaake/RandomScale/tree/master/RandomScale/Data). 
Version 11 of ADMB is required for the TPL files as distributed but you can modify them by replacing instances of PI with the value 3.141592654.

To install ADMB, you need to install a C++ compiler as well as ADMB. We suggest that
you install gcc. To use the built-in links, admb should be installed to c:\admb and the gcc
compiler to c:\mingw. To use different locations you'll need to change function prepare_admb
which is run when the package is attached.

For Windows install [ADMB-11](http://admb-project.googlecode.com/files/admb-11-linux-gcc4.6-32bit.zip) and
[gcc](http://www.admb-project.org/tools/gcc/gcc452-win32.zip/at_download/file)

For other operating systems see (http://www.admb-project.org/downloads) and
[gcc](http://www.admb-project.org/tools/gcc/). Note that prepare_admb() only works for Windows.
 