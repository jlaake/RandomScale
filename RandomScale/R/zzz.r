print.RandomScale.version <- function()
{ library(help=RandomScale)$info[[1]] -> version
	version <- version[pmatch("Version",version)]
	if(!is.null(version))
	{
		um <- strsplit(version," ")[[1]]
		version <- um[nchar(um)>0][2]
	}
	hello <- paste("This is RandomScale ",version,"\n",sep="")
	packageStartupMessage(hello)
}



.onAttach <- function(...) { 
	print.RandomScale.version()
	#prepare_admb()
}
