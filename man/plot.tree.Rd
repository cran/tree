% This file is part of the 'tree' package for R
% Copyright (C) B. D. Ripley 1998-2007
% It is distributed under the GPL version 2 or 3

\name{plot.tree}
\alias{plot.tree}
\title{Plot a Tree Object}
\description{
  Plot a tree object on the current graphical device
}
\usage{
\method{plot}{tree}(x, y = NULL, type = c("proportional", "uniform"), ...)
}
\arguments{
 \item{x}{an object of class \code{"tree"}.}
 \item{y}{ignored. Used for positional matching of \code{type}.}
 \item{type}{character string. If this partially matches
   \code{"uniform"}, the branches are of uniform
   length. Otherwise they are proportional to the decrease in impurity.}
 \item{\dots}{graphical parameters.}
}
\value{
  An (invisible) list with components \code{x} and \code{y}
  giving the coordinates of the tree nodes.

  As a side effect, the value of \code{type == "uniform"} is stored in
  the variable \code{.Tree.unif.?} in the global environment, where \code{?}
  is the device number.
}
\author{B. D. Ripley}


\seealso{\code{\link{tree}}}
\keyword{tree}
\keyword{hplot}
