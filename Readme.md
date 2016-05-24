Random Scale
========

The `RandomScale` package provides R code to fit a random scale half-normal detection function to line transect sampling data. The package also 
contains R code to use [AD Model Builder](http://admb-project.org) with the R package [R2admb](https://github.com/bbolker/R2admb) to fit mixed-effects models with known covariate values and a random component. 
The ADMB TPL files are contained in the [Inst directory](https://github.com/jlaake/RandomScale/tree/master/RandomScale/inst).
Example AMDB DAT files are contained in the [Inst/Data directory](https://github.com/jlaake/RandomScale/tree/master/RandomScale/inst/Data). 
Version 11 of ADMB is required for the TPL files as distributed but you can modify them by replacing instances of PI with the value 3.141592654.

Executable versions are also provided if you download the [Windows binary](https://drive.google.com/folder/d/0B77g1ScdUwVeOVJNUVVGS0YtWE0/edit?ddrp=1#).
In that case you do not need ADMB but the R2admb package is required because it is used to extract results from the output file. The R2admb package should be installed separately.

Examples from Oedekoven et al can be found in the package. To run the hp and simulation example, you'll need to install the mrds package from CRAN.  In turn it needs several other packages to be installed manually. See ?RandomScale for links to the examples.

To install ADMB, you need to install a C++ compiler as well as ADMB. We suggest that
you install gcc. To use the built-in links, [ADMB-11](http://admb-project.googlecode.com/files/admb-11-linux-gcc4.6-32bit.zip) should be installed 
to c:\admb and [gcc compiler](http://www.admb-project.org/tools/gcc/gcc452-win32.zip/at_download/file) should be installed to c:\mingw. 
To use different locations you'll need to change function prepare_admb. prepare_admb() must be run prior to using ADMB without downloaded executables.

For other operating systems see (http://www.admb-project.org/downloads) and
(http://www.admb-project.org/tools/gcc/). Note that prepare_admb() only works for Windows.

Download [Windows package binary](https://docs.google.com/folder/d/0B77g1ScdUwVeOVJNUVVGS0YtWE0/edit). From link, browse to RandomScale and then click on
the version of the package you want. You should see a listing of the package contents as files.  Select File/Download. 
To install in R, from the R menu, use Packages\Install from Local Zip file and browse to location of downloaded zip. 

This repository is a scientific product and is not official communication of the Alaska Fisheries Science Center, the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. AFSC-MML has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
