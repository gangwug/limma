#  ZZZ.R

#.First.lib <- function(libname, pkgname, where)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	23 Jan 2004. Last revised 23 Oct 2004.
#{
#	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) {
#		winMenuAddItem("Vignettes","limma","limmaUsersGuide()")
#	}
#}

.onLoad <- function(libname, pkgname)
#	Load normexp C code
#	Gordon Smyth
#	4 Jan 2009.
{
	library.dynam("limma", pkgname, libname)
}

.onUnload <- function(libpath)
#	Unload normexp C code
#	Gordon Smyth
#	15 Jan 2009.
{
	library.dynam.unload("limma", libpath)
}

.onAttach <- function(libname, pkgname)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	23 Jan 2004. Last revised 4 Jan 2009.
{
	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) winMenuAddItem("Vignettes","limma","limmaUsersGuide()")
}
