% This file is part of the 'tree' package for R
% Copyright (C) B. D. Ripley 1998-2007
% It is distributed under the GPL version 2 or later

\name{deviance.tree}
\alias{deviance.tree}
\alias{deviance.singlenode}
\title{Extract Deviance from a Tree Object}
\description{
  Extract deviance from a tree object.
}
\usage{
\method{deviance}{tree}(object, detail = FALSE, \dots)
}
\arguments{
 \item{object}{an object of calls \code{"tree"}}
 \item{detail}{logical. If true, returns a vector of deviance
   contributions from each node.}
 \item{\dots}{arguments to be passed to or from other methods.}
}
\value{
  The overall deviance, or a vector of contributions from the cases
  at each node. The overall deviance is the sum over leaves in the
  latter case.
}
\keyword{tree}
