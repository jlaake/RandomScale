#' Random Scale Detection Functions for Line Transect Data
#' 
#' The RandomScale package provides R and AD Model Builder(\url{http://admb-project.org})code to fit random effect and mixed effect models for the scale
#' of a half-normal detection function to line transect sampling data. 
#' 
#' The ADMB TPL files are contained in the Inst directory(\url{https://github.com/jlaake/RandomScale/tree/master/RandomScale/inst}).
#' Example AMDB DAT files are contained in the Data directory(\url{https://github.com/jlaake/RandomScale/tree/master/RandomScale/Data}).
#' Version 11 of ADMB is required for the TPL files as distributed but you can modify them by replacing instances of PI with the 
#' value 3.141592654
#'  
#' To install ADMB, you need to install a C++ compiler as well as ADMB. We suggest that
#' you install gcc. To use the built-in links, admb should be installed to c:\admb and the gcc
#' compiler to c:\mingw. To use different locations you'll need to change function prepare_admb in the package
#' which is run when the package is attached.
#' 
#' For Windows install ADMB-11 (\url{http://admb-project.googlecode.com/files/admb-11-linux-gcc4.6-32bit.zip}) and
#' gcc(\url{http://www.admb-project.org/tools/gcc/gcc452-win32.zip/at_download/file})
#' 
#' For other operating systems see (\url{http://www.admb-project.org/downloads}) and
#' gcc(\url{http://www.admb-project.org/tools/gcc/}). Note that prepare_admb() only works for Windows.
#'  
#' For more information see \url{https://github.com/downloads/jlaake/ADMB-Examples/distance_random_effect.pdf}
#' @name RandomScale-package
#' @aliases RandomScale
#' @author Jeff Laake
#' @references Fournier, D.A., H.J. Skaug, J. Ancheta, J. Ianelli, A. Magnusson, M.N. Maunder, A. Nielsen, and J. Sibert. 2012. AD Model Builder: using automatic differentiation for statistical inference of highly parameterized complex nonlinear models. Optim. Methods Softw. 27:233-249.
NULL

