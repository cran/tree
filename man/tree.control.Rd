\name{tree.control}
\title{Select Parameters for Tree}
\usage{
tree.control(nobs, mincut=5, minsize=10, mindev=0.01)
}
\alias{tree.control}
\arguments{
 \item{nobs}{The number of observations in the training set.}
 \item{mincut}{The minimum number of observations to include in either
   child node. This is a weighted quantity; the observational weights are used
   to compute the `number'. The default is 5.}
 \item{minsize}{The smallest allowed node node: a weighted quantity. The
   default is 10.}
 \item{mindev}{The within-node deviance must be at least this times that
 of the root node for the node to be split.}
}
\description{
  The produces default values of \code{mincut} and \code{minsize}, and
  ensures that \code{mincut} is at least half \code{minsize}.

  To produce a tree that perfectly fits the data, set \code{mindev=-1}
  and \code{minsize = 2}.
}
\value{
  A list
  \item{mincut}{The maximum of the input or default mincut and 1}
  \item{minsize}{The maximum of the input or default minsize and 2.}
  \item{nmax}{A estimate of the maximum number of nodes that might be grown.}
  \item{nobs}{The input \code{nobs}.}
}
\author{B.D. Ripley}

\seealso{\code{\link{tree}}}

\keyword{tree}
