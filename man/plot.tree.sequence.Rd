% This file is part of the 'tree' package for R
% Copyright (C) B. D. Ripley 1998-2007
% It is distributed under the GPL version 2 or 3

\name{plot.tree.sequence}
\alias{plot.tree.sequence}
\title{
  Plot a Tree Sequence
}
\description{
  Allows the user to plot a tree sequence.
}
\usage{
\method{plot}{tree.sequence}(x, \dots, type = "l", ylim = range(x$dev),
    order = c("increasing", "decreasing"))
}
\arguments{
  \item{x}{object of class \code{tree.sequence}.
    This is assumed to be the result of some function that produces
    an object with the same named components (\code{size},
    \code{deviance}, \code{k}) as that returned by \code{\link{prune.tree}}.}
  \item{order}{of \code{size} on the plot.  Use \code{"decreasing"} for
    the natural ordering of \code{k} and the amount of pruning.  Only
    the first character is needed.}
  \item{type, ylim, \dots}{graphical parameters.}
}
\section{Side Effects}{
  Plots deviance or number of misclassifications (or total loss)
  versus size for a sequence of trees.
}
\details{
  This function is a method for the generic function
  \code{plot()} for class \code{tree.sequence}.
  It can be invoked by calling \code{plot(x)} for an
  object \code{x} of the appropriate class, or directly by
  calling \code{plot.tree.sequence(x)} regardless of the
  class of the object.
}
\examples{
data(cpus, package="MASS")
cpus.ltr <- tree(log(perf) ~ syct + mmin + mmax + cach + chmin + chmax,
                 data = cpus)
plot(prune.tree(cpus.ltr))
}
\keyword{tree}
