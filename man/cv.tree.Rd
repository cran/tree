\name{cv.tree}
\title{Cross-validation for Choosing Tree Complexity}
\usage{
cv.tree(object, rand, FUN = prune.tree, K = 10, ...)
}
\alias{cv.tree}
\arguments{
  \item{object}{An object of class \code{"tree"}.}
  \item{rand}{Optionally an integer vector of the length the number of
    cases used to create \code{object}, assigning the cases to different
    groups for cross-validation.}
  \item{FUN}{The function to do the pruning.}
  \item{K}{The number of folds of the cross-validation.}
  \item{\dots}{Additional arguments to \code{FUN}.}
}
\description{
  Runs a K-fold cross-validation experiment to find the deviance or
  number of misclassifications as a function of the cost-complexity
  parameter \code{k}.
}
\value{
  A copy of \code{FUN} applied to \code{object}, with component
  \code{dev} replaced by the cross-validated results from the
  sum of the \code{dev} components of each fit.
}
\author{B. D. Ripley}

\seealso{\code{\link{tree}}, \code{\link{prune.tree}}}

\examples{
library(MASS)
data(cpus)
cpus.ltr <- tree(log10(perf) ~ syct + mmin + mmax + cach
     + chmin + chmax, data=cpus)
cv.tree(cpus.ltr, , prune.tree)
}
\keyword{tree}
