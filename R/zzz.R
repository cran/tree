.First.lib <- function(lib, pkg) library.dynam("tree", pkg, lib)

if(!exists("na.pass", envir=NULL)) na.pass <- function(x) x

