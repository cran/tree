.First.lib <- function(lib, pkg) library.dynam("tree", pkg, lib)

if(version$major == "0" && version$minor < "0.63")
  labels <- function(object, ...) UseMethod("labels")

