\name{deviance.tree}
\alias{deviance.tree}
\title{Extract Deviance from a Tree Object}
\description{
  Extract deviance from a tree object.
}
\usage{
deviance.tree(object, detail = FALSE)
}
\arguments{
 \item{object}{an object of calls \code{"tree"}}
 \item{detail}{logical. If true, returns a vector of deviance
   contributions from each node.}
}
\value{
  The overall deviance, or a vector of contributions from each the cases
  at each node. The overall deviance is the sum over leaves in the
  latter case.
}
\keyword{tree}
