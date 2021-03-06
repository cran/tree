% This file is part of the 'tree' package for R
% Copyright (C) B. D. Ripley 1998-2016
% It is distributed under the GPL version 2 or 3

\name{tree.control}
\title{Select Parameters for Tree}
\usage{
tree.control(nobs, mincut = 5, minsize = 10, mindev = 0.01)
}
\alias{tree.control}
\arguments{
  \item{nobs}{The number of observations in the training set.}
  \item{mincut}{The minimum number of observations to include in either
    child node.  This is a weighted quantity; the observational weights are used
    to compute the \sQuote{number}. The default is 5.}
  \item{minsize}{The smallest allowed node size: a weighted quantity.  The
    default is 10.}
  \item{mindev}{The within-node deviance must be at least this times that
    of the root node for the node to be split.}
}
\description{
  A utility function for use with the \code{control} argument of \code{tree}.
}
\details{ 
  This function produces default values of \code{mincut} and
  \code{minsize}, and ensures that \code{mincut} is at most half
  \code{minsize}.

  To produce a tree that fits the data perfectly, set \code{mindev = 0}
  and \code{minsize = 2}, if the limit on tree depth allows such a tree.
}
\value{
  A list:
  \item{mincut}{The maximum of the input or default \code{mincut} and 1}
  \item{minsize}{The maximum of the input or default \code{minsize} and 2.}
  \item{nmax}{A estimate of the maximum number of nodes that might be grown.}
  \item{nobs}{The input \code{nobs}.}
}
\note{
  The interpretation of \code{mindev} given here is that of Chambers and
  Hastie (1992, p. 415), and apparently not what is actually implemented
  in S.  It seems S uses an absolute bound.
}
\author{B. D. Ripley}

\references{
  Chambers, J. M. and Hastie, T. J. (1992)
  \emph{Statistical Models in S}. Wadsworth & Brooks/Cole.
}


\seealso{\code{\link{tree}}}

\keyword{tree}
