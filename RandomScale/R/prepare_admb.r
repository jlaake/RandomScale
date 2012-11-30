#' Prepare to use ADMB
#' 
#' Sets environment variables for admb if Windows
#' 
#' @export
#' @return NULL
#' 
prepare_admb=function()
{
	if(R.Version()$os=="mingw32")
	{
		Sys.setenv(PATH = paste("c:/admb/bin;c:admb/utilities;c:/MinGW/bin;",Sys.getenv("PATH"),sep=";"))
		Sys.setenv(ADMB_HOME = "c:/admb")
	}
	invisible()
}
